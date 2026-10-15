#!/usr/bin/env python3
"""
Detect synteny breakpoints from AGORA ancestral genome reconstruction output.

Sorts genes by each of three genomes (ANC, Q1, Q2) and detects where the
other genomes change chromosome — indicating a breakpoint event.

Classification logic
--------------------
Sorted by ANC (per ancestral chromosome):
  - Q1 and Q2 both change chr  → shared_fission
  - only Q1 changes            → fission_Q1
  - only Q2 changes            → fission_Q2

Sorted by Q1 / Q2 (per focal chromosome):
  - ANC and non-focal both change → fusion_<focal>
  - only ANC changes, same gene order in non-focal → shared_fusion
  - only ANC changes, opposite gene order in non-focal → fusion_Q1 + fusion_Q2
  - only non-focal changes        → fission_<non-focal>

Usage
-----
    python detect_breakpoints_agora_v2.py -i input.bz2 --q1 speciesA --q2 speciesB -o all.tsv --bp-out bp_only.tsv
"""

import argparse
import bz2
import logging
import os
import pathlib
import sys

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

LARGE_INT = 10**12
MISSING = {"", ".", "NA", None}


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def open_maybe_bz2(path: str):
    """Open plain-text or bz2-compressed files transparently."""
    if not os.path.isfile(path):
        sys.exit(f"ERROR: input file not found: {path}")
    return bz2.open(path, "rt") if path.endswith(".bz2") else open(path, "r")


def parse_token(tok: str) -> tuple:
    """
    Parse an AGORA gene token.
    Format: <species>.<geneID>_…_SUPER_<chr>_<start>-<end>
    Returns (species, gene_id, chr, start, end).  gene_id is the part
    before the first underscore after the species dot.
    """
    species, rest = tok.split(".", 1) if "." in tok else ("", tok)

    if "_SUPER_" not in rest:
        return species, "", "", "", ""

    gene_part, tail = rest.split("_SUPER_", 1)
    gene_id = gene_part.split("_", 1)[0] if gene_part else ""

    parts = tail.split("_", 1)
    chr_ = parts[0]
    start, end = "", ""
    if len(parts) > 1 and "-" in parts[1]:
        start, end = parts[1].split("-", 1)

    return species, gene_id, chr_, start, end


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_input(path: str, q1: str, q2: str) -> pd.DataFrame:
    """Read AGORA file and return a DataFrame with one row per gene."""
    records = []
    n_skipped = 0

    with open_maybe_bz2(path) as fh:
        for lineno, line in enumerate(fh, 1):
            fields = line.rstrip("\n").split()
            if not fields:
                n_skipped += 1
                continue

            anc_chr = fields[0]
            anc_start = fields[1] if len(fields) > 1 else ""
            anc_end = fields[2] if len(fields) > 2 else ""

            q1_chr = q1_start = q1_end = ""
            q2_chr = q2_start = q2_end = ""
            gene_ids = set()

            for tok in fields[5:]:
                species, gid, chr_, start, end = parse_token(tok)
                if gid:
                    gene_ids.add(gid)
                if species == q1 and not q1_chr:
                    q1_chr, q1_start, q1_end = chr_, start, end
                elif species == q2 and not q2_chr:
                    q2_chr, q2_start, q2_end = chr_, start, end

            gene_id = ";".join(sorted(gene_ids)) if gene_ids else ""

            records.append(
                {
                    "lineno": lineno,
                    "anc_chr": anc_chr,
                    "anc_start": anc_start,
                    "anc_end": anc_end,
                    "gene_id": gene_id,
                    "q1_chr": q1_chr,
                    "q1_start": q1_start,
                    "q1_end": q1_end,
                    "q2_chr": q2_chr,
                    "q2_start": q2_start,
                    "q2_end": q2_end,
                    "raw": line.rstrip("\n"),
                }
            )

    if n_skipped:
        log.warning("Skipped %d empty lines", n_skipped)

    df = pd.DataFrame(records)

    if df.empty:
        sys.exit("ERROR: no records parsed — check the input file and species names")

    # replace empty strings with NA for cleaner handling
    for col in ("anc_chr", "q1_chr", "q2_chr"):
        df[col] = df[col].replace("", pd.NA)

    # convert coordinate columns to nullable integers
    for col in ("anc_start", "anc_end", "q1_start", "q1_end", "q2_start", "q2_end"):
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    # remove rows where multiple gene IDs mapped (different BUSCO genes
    # for Q1 vs Q2 in the same ancestral locus)
    multi_gene = df["gene_id"].str.contains(";", na=False)
    n_multi = multi_gene.sum()
    if n_multi:
        log.info("Removed %d rows with multiple gene IDs", n_multi)
        df = df.loc[~multi_gene].reset_index(drop=True)

    # initialise breakpoint columns
    df["bp_ids"] = [[] for _ in range(len(df))]
    df["bp_class"] = [[] for _ in range(len(df))]

    n_q1 = df["q1_chr"].notna().sum()
    n_q2 = df["q2_chr"].notna().sum()
    log.info(
        "Parsed %d genes (%d with Q1 coords, %d with Q2 coords)",
        len(df), n_q1, n_q2,
    )

    return df


# ---------------------------------------------------------------------------
# Overlap filtering
# ---------------------------------------------------------------------------

def remove_overlapping_genes(df: pd.DataFrame) -> pd.DataFrame:
    """Remove rows whose gene coordinates overlap any other gene on the
    same chromosome within Q1 or Q2 genomes.

    ANC coordinates are ordinal positions (not genomic bp) so they are
    not checked for overlap.

    Two genes overlap when they share at least one base on the same
    chromosome: not (a_end < b_start or b_end < a_start).

    Both members of every overlapping pair are removed.
    """
    bad = set()  # indices to drop

    for prefix in ("q1", "q2"):
        chr_col = f"{prefix}_chr"
        start_col = f"{prefix}_start"
        end_col = f"{prefix}_end"

        # work only with rows that have coords for this genome
        mask = df[chr_col].notna() & df[start_col].notna() & df[end_col].notna()
        sub = df.loc[mask, [chr_col, start_col, end_col]].copy()
        if sub.empty:
            continue

        sub[start_col] = sub[start_col].astype(int)
        sub[end_col] = sub[end_col].astype(int)

        # sort by (chr, start) then sweep
        sub = sub.sort_values([chr_col, start_col])
        indices = sub.index.to_numpy()
        chrs = sub[chr_col].to_numpy()
        starts = sub[start_col].to_numpy()
        ends = sub[end_col].to_numpy()

        n = len(indices)
        for i in range(n):
            # compare with following genes on the same chromosome
            for j in range(i + 1, n):
                if chrs[j] != chrs[i]:
                    break  # different chromosome, no more overlaps
                if starts[j] > ends[i]:
                    break  # sorted by start: no further overlap possible
                # overlap exists
                bad.add(indices[i])
                bad.add(indices[j])

    if bad:
        log.info(
            "Removed %d rows with overlapping gene coordinates", len(bad)
        )
    return df.drop(index=bad).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Breakpoint detection
# ---------------------------------------------------------------------------

def _chr_val(series_val):
    """Return None for missing/empty chromosome values, else the value."""
    if pd.isna(series_val):
        return None
    v = str(series_val).strip()
    return None if v in MISSING else v


def _sort_key(df: pd.DataFrame, chr_col: str, start_col: str) -> pd.DataFrame:
    """Return df sorted by (chr, start) with missing values last."""
    tmp = df.copy()
    tmp["_sort_chr"] = tmp[chr_col].fillna("ZZZ")
    tmp["_sort_start"] = tmp[start_col].fillna(LARGE_INT)
    tmp = tmp.sort_values(["_sort_chr", "_sort_start"]).drop(
        columns=["_sort_chr", "_sort_start"]
    )
    return tmp


def detect_breakpoints(df: pd.DataFrame, focal: str, counter: int = 1) -> int:
    """
    Sort by `focal` genome, walk consecutive rows within each focal chromosome,
    and emit breakpoints where the other genomes break contiguity.

    A contiguity break is either:
      - a chromosome change in a non-focal genome, OR
      - a direction reversal (ascending ↔ descending start positions)
        within the same chromosome in a non-focal genome.

    Modifies df in-place (bp_ids, bp_class).  Returns next counter value.
    """
    if focal == "ANC":
        chr_col, start_col = "anc_chr", "anc_start"
        other_cols = {"Q1": "q1_chr", "Q2": "q2_chr"}
    elif focal == "Q1":
        chr_col, start_col = "q1_chr", "q1_start"
        other_cols = {"ANC": "anc_chr", "Q2": "q2_chr"}
    elif focal == "Q2":
        chr_col, start_col = "q2_chr", "q2_start"
        other_cols = {"ANC": "anc_chr", "Q1": "q1_chr"}
    else:
        sys.exit(f"ERROR: unknown focal genome '{focal}' (expected ANC, Q1, Q2)")

    # only consider rows with all three genomes present
    complete = (
        df["anc_chr"].notna() & df["q1_chr"].notna() & df["q2_chr"].notna()
    )
    idx_complete = df.index[complete]

    # Block tracking: consecutive genes without a break share a block_id
    block_col = f"block_{focal.lower()}"
    df[block_col] = pd.NA

    if idx_complete.empty:
        log.warning("No complete rows — skipping %s pass", focal)
        return counter

    sorted_idx = (
        df.loc[idx_complete]
        .pipe(_sort_key, chr_col, start_col)
        .index
    )

    # Build start-column lookup for direction tracking
    other_keys = list(other_cols.keys())   # e.g. ['Q1','Q2'] when focal=ANC
    other_chr_cols = list(other_cols.values())
    other_start_cols = [
        c.replace("_chr", "_start") for c in other_chr_cols
    ]

    # Per-genome traversal direction: +1 ascending, -1 descending, None unknown.
    # Reset when the focal chromosome changes.
    prev_dir = [None] * len(other_keys)
    prev_i = None
    block_id = 1

    for curr_i in sorted_idx:
        if prev_i is None:
            df.at[curr_i, block_col] = block_id
            prev_i = curr_i
            continue

        # new focal chromosome → reset direction tracking, new block
        if _chr_val(df.at[prev_i, chr_col]) != _chr_val(df.at[curr_i, chr_col]):
            block_id += 1
            df.at[curr_i, block_col] = block_id
            prev_dir = [None] * len(other_keys)
            prev_i = curr_i
            continue

        prev_vals = [_chr_val(df.at[prev_i, c]) for c in other_chr_cols]
        curr_vals = [_chr_val(df.at[curr_i, c]) for c in other_chr_cols]

        # skip if any of the four values is missing (missing ≠ breakpoint)
        if None in prev_vals or None in curr_vals:
            df.at[curr_i, block_col] = block_id
            prev_i = curr_i
            continue

        # Detect breaks: chromosome change OR direction reversal within
        # the same chromosome.  Track reason per genome.
        # reason: None = no break, 'chr' = chr change, 'dir' = direction reversal
        changed = [False] * len(other_keys)
        change_reason = [None] * len(other_keys)
        for k in range(len(other_keys)):
            if prev_vals[k] != curr_vals[k]:
                # chromosome change → definite break
                changed[k] = True
                change_reason[k] = 'chr'
                prev_dir[k] = None  # reset direction
            else:
                # same chromosome: check traversal direction
                p_start = df.at[prev_i, other_start_cols[k]]
                c_start = df.at[curr_i, other_start_cols[k]]
                if pd.notna(p_start) and pd.notna(c_start):
                    delta = int(c_start) - int(p_start)
                    if delta == 0:
                        pass  # identical position – no change
                    else:
                        curr_sign = 1 if delta > 0 else -1
                        if prev_dir[k] is not None and curr_sign != prev_dir[k]:
                            # direction reversal within the same chr
                            changed[k] = True
                            change_reason[k] = 'dir'
                        prev_dir[k] = curr_sign

        # For Q1/Q2 passes: when only ANC changes, check whether the
        # non-focal query has genes in the same order as the focal query.
        # Same order → shared_fusion (single event in common ancestor).
        # Opposite order → independent fusions in both queries.
        nonf_same_direction = None
        if focal in ("Q1", "Q2"):
            nonf_key = other_keys[1]  # 'Q1' or 'Q2'
            nonf_start_col = f"{nonf_key.lower()}_start"
            p_nonf_start = df.at[prev_i, nonf_start_col]
            c_nonf_start = df.at[curr_i, nonf_start_col]
            if pd.notna(p_nonf_start) and pd.notna(c_nonf_start):
                nonf_same_direction = (int(p_nonf_start) < int(c_nonf_start))

        # classify
        labels = _classify(focal, other_keys, changed, change_reason, nonf_same_direction)
        if not labels:
            df.at[curr_i, block_col] = block_id
            prev_i = curr_i
            continue

        for label in labels:
            bp_id = f"{focal}_BP{counter:06d}"
            counter += 1

            df.at[prev_i, "bp_ids"].append(bp_id)
            df.at[curr_i, "bp_ids"].append(bp_id)
            df.at[prev_i, "bp_class"].append(label)
            df.at[curr_i, "bp_class"].append(label)

        # start a new block after a breakpoint
        block_id += 1
        df.at[curr_i, block_col] = block_id

        prev_i = curr_i

    return counter


def _classify(
    focal: str,
    other_keys: list,
    changed: list,
    change_reason: list,
    nonf_same_direction: bool | None = None,
) -> list[str]:
    """Return list of breakpoint labels (empty if no breakpoint).

    Parameters
    ----------
    focal : str
        The genome used for sorting ('ANC', 'Q1', or 'Q2').
    other_keys : list[str]
        The two other genome labels, e.g. ['Q1','Q2'] when focal='ANC'.
    changed : list[bool]
        Whether each of the other genomes broke contiguity.
    change_reason : list[str | None]
        Per-genome reason: 'chr' for chromosome change, 'dir' for
        direction reversal only, None if no break.
    nonf_same_direction : bool | None
        When focal is Q1/Q2 and only ANC changed: True if the non-focal
        query has the two genes in the same order as the focal query
        (→ shared_fusion), False if reversed (→ independent fusions in
        both queries).  None if coordinates were missing.
    """
    n_changed = sum(changed)
    if n_changed == 0:
        return []

    if focal == "ANC":
        # other_keys = ['Q1','Q2']
        if all(changed):
            # if both are direction reversals only → shared_inversion
            if change_reason[0] == 'dir' and change_reason[1] == 'dir':
                return ["shared_inversion"]
            # if both are chr changes → shared_fission
            if change_reason[0] == 'chr' and change_reason[1] == 'chr':
                return ["shared_fission"]
            # mixed: one chr change, one inversion
            labels = []
            for k in range(2):
                if change_reason[k] == 'dir':
                    labels.append(f"inversion_{other_keys[k]}")
                else:
                    labels.append(f"fission_{other_keys[k]}")
            return labels
        # only one changed
        idx = 0 if changed[0] else 1
        if change_reason[idx] == 'dir':
            return [f"inversion_{other_keys[idx]}"]
        return [f"fission_{other_keys[idx]}"]

    # focal is Q1 or Q2; other_keys = ['ANC', '<non-focal>']
    anc_changed = changed[0]
    nonf_changed = changed[1]
    anc_reason = change_reason[0]
    nonf_reason = change_reason[1]
    nonf = other_keys[1]

    if anc_changed and nonf_changed:
        # both break: if non-focal is only a direction reversal
        if nonf_reason == 'dir' and anc_reason == 'chr':
            return [f"fusion_{focal}", f"inversion_{nonf}"]
        if nonf_reason == 'chr' and anc_reason == 'dir':
            return [f"inversion_{focal}", f"fission_{nonf}"]
        if nonf_reason == 'dir' and anc_reason == 'dir':
            return [f"inversion_{focal}", f"inversion_{nonf}"]
        # both chr changes
        return [f"fusion_{focal}"]

    if anc_changed and not nonf_changed:
        if anc_reason == 'dir':
            return [f"inversion_{focal}"]
        # chr change in ANC only → check directionality for shared vs independent fusion
        if nonf_same_direction is None:
            return ["missing_data"]
        if nonf_same_direction:
            return ["shared_fusion"]
        else:
            return ["fusion_Q1", "fusion_Q2"]

    if nonf_changed and not anc_changed:
        if nonf_reason == 'dir':
            return [f"inversion_{nonf}"]
        return [f"fission_{nonf}"]

    return []


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def build_output(df: pd.DataFrame, q1: str, q2: str) -> pd.DataFrame:
    """Prepare a clean output DataFrame."""
    out = df[
        [
            "anc_chr", "anc_start", "anc_end", "gene_id",
            "q1_chr", "q1_start", "q1_end",
            "q2_chr", "q2_start", "q2_end",
            "block_anc", "block_q1", "block_q2",
            "bp_ids", "bp_class",
        ]
    ].copy()

    out.insert(7, "q1_species", q1)
    out.insert(11, "q2_species", q2)


    # collapse list columns to semicolon-separated strings
    out["bp_ids"] = out["bp_ids"].apply(lambda x: ";".join(x) if x else "NA")
    out["bp_class"] = out["bp_class"].apply(
        lambda x: ";".join(sorted(set(x))) if x else "NA"
    )

    # fill remaining NAs
    out = out.fillna("NA")

    return out


def load_chr_sizes(bed_dir: str, species: str) -> dict[str, int]:
    """Load chromosome sizes from a BED file (chr  0  size).

    The SUPER_ prefix is stripped so chromosome names match the data.
    Returns {chr_name: size}.
    """
    path = pathlib.Path(bed_dir) / f"{species}.bed"
    sizes: dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            name = parts[0].replace("SUPER_", "")
            sizes[name] = int(parts[2])
    log.info("Loaded %d chromosome sizes from %s", len(sizes), path)
    return sizes


def build_event_table(
    df: pd.DataFrame,
    q1: str,
    q2: str,
    chr_sizes: dict[str, dict[str, int]],
) -> pd.DataFrame:
    """Build a per-event breakpoint table with coordinates in all three genomes.

    For each breakpoint event (bp_id), the two flanking genes that share that
    bp_id define the breakpoint region in every genome:

    - Same chromosome: the gap between the two genes (closest endpoints).
    - Different chromosomes (fission): first check if the gene has co-located
      events (other bp_ids) with a partner on the SAME chromosome — if so, use
      the gap to that partner gene.  Otherwise, two rows from the gene to the
      chromosome terminus (using actual chromosome sizes from BED files).

    Parameters
    ----------
    chr_sizes : dict mapping genome prefix ("q1"/"q2") to {chr: length}.
    """
    chr_len: dict[str, dict[str, int]] = {"q1": {}, "q2": {}}
    chr_len.update(chr_sizes)

    # Map bp_id → {class, rows: [idx1, idx2]}
    event_map: dict[str, dict] = {}
    for idx in df.index:
        for bp_id, bp_cl in zip(df.at[idx, "bp_ids"], df.at[idx, "bp_class"]):
            ev = event_map.setdefault(bp_id, {"class": bp_cl, "rows": []})
            if idx not in ev["rows"]:
                ev["rows"].append(idx)

    # Reverse map: row_idx → set of bp_ids for that gene
    gene_bpids: dict[int, set] = {}
    for idx in df.index:
        bps = df.at[idx, "bp_ids"]
        if bps:
            gene_bpids[idx] = set(bps)

    GENOME_COLS = [
        ("anc", "anc_chr", "anc_start", "anc_end"),
        ("q1",  "q1_chr",  "q1_start",  "q1_end"),
        ("q2",  "q2_chr",  "q2_start",  "q2_end"),
    ]

    def _colocated_gap(
        row_idx: int, bp_id: str,
        chr_col: str, start_col: str, end_col: str,
        this_chr: str, this_start: int, this_end: int,
    ) -> tuple | None:
        """Check if *row_idx* has another bp_id whose partner gene is on
        *this_chr* in the same genome.  If so, return (chr, gap_start, gap_end).
        """
        other_bpids = gene_bpids.get(row_idx, set()) - {bp_id}
        for other_bp in other_bpids:
            other_info = event_map.get(other_bp)
            if not other_info:
                continue
            for other_row in other_info["rows"]:
                if other_row == row_idx:
                    continue
                oc = _chr_val(df.at[other_row, chr_col])
                if oc != this_chr:
                    continue
                os_ = df.at[other_row, start_col]
                oe_ = df.at[other_row, end_col]
                if pd.isna(os_) or pd.isna(oe_):
                    continue
                os_, oe_ = int(os_), int(oe_)
                bp_s = min(this_end, oe_)
                bp_e = max(this_start, os_)
                return (this_chr, bp_s, bp_e)
        return None

    records = []
    for bp_id in sorted(event_map):
        info = event_map[bp_id]
        bp_cls = info["class"]
        rows = info["rows"]
        if len(rows) != 2:
            log.warning("Skipping %s: expected 2 gene rows, got %d", bp_id, len(rows))
            continue

        g1, g2 = df.loc[rows[0]], df.loc[rows[1]]

        # For each genome compute breakpoint location(s)
        genome_locs: dict[str, list[tuple]] = {}
        has_multi = False

        for gname, chr_col, start_col, end_col in GENOME_COLS:
            c1, c2 = _chr_val(g1[chr_col]), _chr_val(g2[chr_col])
            s1, e1 = g1[start_col], g1[end_col]
            s2, e2 = g2[start_col], g2[end_col]

            if c1 is None or c2 is None or pd.isna(s1) or pd.isna(e1) or pd.isna(s2) or pd.isna(e2):
                genome_locs[gname] = [("NA", "NA", "NA")]
                continue

            s1, e1, s2, e2 = int(s1), int(e1), int(s2), int(e2)

            if c1 == c2:
                # Same chromosome → gap between closest endpoints
                bp_s = min(e1, e2)
                bp_e = max(s1, s2)
                genome_locs[gname] = [(c1, bp_s, bp_e)]
            else:
                # Different chromosomes → check co-located events first
                has_multi = True
                if gname == "anc":
                    genome_locs[gname] = [(c1, s1, e1), (c2, s2, e2)]
                else:
                    locs = []
                    for c, s, e, ridx in [
                        (c1, s1, e1, rows[0]),
                        (c2, s2, e2, rows[1]),
                    ]:
                        # Try co-located event first
                        gap = _colocated_gap(
                            ridx, bp_id, chr_col, start_col, end_col, c, s, e,
                        )
                        if gap:
                            locs.append(gap)
                        else:
                            # Fall back to terminus
                            cl = chr_len[gname].get(c, e)
                            if s <= (cl - e):
                                locs.append((c, 0, s))
                            else:
                                locs.append((c, e, cl))
                    genome_locs[gname] = locs

        if not has_multi:
            rec = {"bp_id": bp_id, "bp_class": bp_cls}
            for gname in ("anc", "q1", "q2"):
                c, s, e = genome_locs[gname][0]
                rec[f"{gname}_chr"] = c
                rec[f"{gname}_bp_start"] = s
                rec[f"{gname}_bp_end"] = e
            records.append(rec)
        else:
            for i in range(2):
                rec = {"bp_id": bp_id, "bp_class": bp_cls}
                for gname in ("anc", "q1", "q2"):
                    locs = genome_locs[gname]
                    loc = locs[i] if i < len(locs) else locs[0]
                    c, s, e = loc
                    rec[f"{gname}_chr"] = c
                    rec[f"{gname}_bp_start"] = s
                    rec[f"{gname}_bp_end"] = e
                records.append(rec)

    result = pd.DataFrame(records)
    if result.empty:
        return result

    # Insert species name columns
    q1_pos = list(result.columns).index("q1_chr")
    q2_pos = list(result.columns).index("q2_chr")
    result.insert(q1_pos, "q1_species", q1)
    result.insert(q2_pos + 1, "q2_species", q2)

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Detect synteny breakpoints from AGORA output."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="AGORA ancestral genome file (.txt or .bz2)",
    )
    parser.add_argument("--q1", required=True, help="Species name for query 1")
    parser.add_argument("--q2", required=True, help="Species name for query 2")
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV with all genes, summary of events and bedfiles for each genome",
    )
    parser.add_argument(
        "--bed-dir", required=True,
        help="Directory with <species>.bed files (chr  0  size) for chromosome lengths",
    )
    args = parser.parse_args()

    # parse
    df = parse_input(args.input, args.q1, args.q2)

    # remove genes that overlap any other gene in any genome
    df = remove_overlapping_genes(df)

    if df.empty:
        sys.exit("ERROR: no rows remain after removing overlapping genes")

    # filter: keep only rows where all three genomes have chromosome and coordinates
    n_before = len(df)
    complete = (
        df["anc_chr"].notna()
        & df["q1_chr"].notna()
        & df["q2_chr"].notna()
        & df["anc_start"].notna()
        & df["q1_start"].notna()
        & df["q2_start"].notna()
    )
    df = df.loc[complete].reset_index(drop=True)
    n_after = len(df)
    log.info(
        "Filtered incomplete rows: %d → %d (%d removed)",
        n_before, n_after, n_before - n_after,
    )

    if df.empty:
        sys.exit("ERROR: no complete rows remain after filtering — check species names")

    # re-initialise bp columns after filtering
    df["bp_ids"] = [[] for _ in range(len(df))]
    df["bp_class"] = [[] for _ in range(len(df))]

    # detect breakpoints in three passes
    counter = 1
    for focal in ("ANC", "Q1", "Q2"):
        counter = detect_breakpoints(df, focal, counter)

    total_bp = (df["bp_ids"].apply(len) > 0).sum()
    log.info("Total breakpoint events: %d (rows involved: %d)", counter - 1, total_bp)


    # build output
    out = build_output(df, args.q1, args.q2)

    # if directory is specified, create it if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # write full table
    out.to_csv(args.output, sep="\t", index=False)
    log.info("Wrote all genes      → %s  (%d rows)", args.output, len(out))

    # load chromosome sizes from BED files
    chr_sizes = {
        "q1": load_chr_sizes(args.bed_dir, args.q1),
        "q2": load_chr_sizes(args.bed_dir, args.q2),
    }

    # build and write breakpoint event table + per-genome BED files
    events = build_event_table(df, args.q1, args.q2, chr_sizes)
    base, ext = os.path.splitext(args.output)
    events_path = f"{base}.events{ext}"
    events.to_csv(events_path, sep="\t", index=False)
    log.info("Wrote event table    → %s  (%d rows)", events_path, len(events))

    for gname, species in [("q1", args.q1), ("q2", args.q2)]:
        bp_chr = f"{gname}_chr"
        bp_start = f"{gname}_bp_start"
        bp_end = f"{gname}_bp_end"
        bed = events[[bp_chr, bp_start, bp_end,"bp_id","bp_class"]].copy()
        bed = bed[bed[bp_chr] != "NA"]
        bed = bed.drop_duplicates()
        bed[bp_start] = bed[bp_start].apply(
            lambda x: max(0, int(x) - 1) if x != "NA" else 0
        )
        bed = bed.sort_values([bp_chr, bp_start]).reset_index(drop=True)
        bed_path = f"{base}.{species}_breakpoints.bed"
        bed.to_csv(bed_path, sep="\t", index=False, header=False)
        median_len = int(np.median(bed[bp_end].astype(int) - bed[bp_start].astype(int))) if not bed.empty else 0
        max_len = int(np.max(bed[bp_end].astype(int) - bed[bp_start].astype(int))) if not bed.empty else 0
        min_len = int(np.min(bed[bp_end].astype(int) - bed[bp_start].astype(int))) if not bed.empty else 0
        log.info(
            "Wrote %s breakpoint BED    → %s  (%d rows)  [median length: %d bp, min length: %d bp, max length: %d bp]",
            gname.upper(), bed_path, len(bed), median_len, min_len, max_len,
        )


if __name__ == "__main__":
    main()
