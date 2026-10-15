#!/usr/bin/env python3
"""
Refine AGORA breakpoints using cactus conserved alignment blocks.

For each breakpoint region in the per-genome BED files produced by
detect_breakpoints_agora.py, narrow the interval to the tightest gap
between flanking conserved alignment blocks from cactus.

Inputs
------
  --input-bed-q1       events BED for Q1  (chr  start  end  bp_id  bp_class)
  --input-bed-q2       events BED for Q2  (chr  start  end  bp_id  bp_class)
  --blocks-q1          conserved_blocks_by_query1.tsv  (sorted by Q1 coords)
  --blocks-q2          conserved_blocks_by_query2.tsv  (sorted by Q2 coords)

Outputs
-------
  --bed-q1             BED file of refined breakpoints in Q1 coordinates
  --bed-q2             BED file of refined breakpoints in Q2 coordinates
  
  Additional outputs (same format as --bed-q1/2 but with _anc_breakpoints.bed suffix) are the refined breakpoints before merging/filtering by genome.
  
Usage
-----
    python refine_breakpoints.py \\
        --input-bed-q1 gene_alignment_q1_events.bed \\
        --input-bed-q2 gene_alignment_q2_events.bed \\
        --blocks-q1 conserved_blocks_by_query1.tsv \\
        --blocks-q2 conserved_blocks_by_query2.tsv \\
        --bed-q1 refined_q1.bed \\
        --bed-q2 refined_q2.bed
"""

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Load input BED files (from detect_breakpoints_agora.py)
# ---------------------------------------------------------------------------

def load_input_bed(path: str) -> pd.DataFrame:
    """Load a headerless BED file (chr  start  end  bp_id  bp_class)."""
    if not os.path.isfile(path):
        sys.exit(f"ERROR: input BED not found: {path}")

    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["chr", "start", "end", "bp_id", "bp_class"],
        dtype={"chr": str, "bp_id": str},
    )
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    log.info("Loaded %d breakpoint regions from %s", len(df), path)
    return df


# ---------------------------------------------------------------------------
# Load conserved blocks
# ---------------------------------------------------------------------------

def load_blocks(path: str, genome: str) -> pd.DataFrame:
    """Load conserved block TSV and return a DataFrame sorted by the
    indicated genome's coordinates.

    The SUPER_ prefix is stripped from chromosome names so they match
    the BED files produced by detect_breakpoints_agora.py.

    Parameters
    ----------
    path : str
        Path to conserved_blocks_by_query{1,2}.tsv.
    genome : str
        'q1' or 'q2' — which genome's coordinates to use as reference.

    Returns a DataFrame with columns:
        ref_chr, ref_start, ref_end, other_chr, other_start, other_end,
        strand, block_id
    """
    if not os.path.isfile(path):
        sys.exit(f"ERROR: blocks file not found: {path}")

    df = pd.read_csv(path, sep="\t")

    # Standardise column names
    if genome == "q1":
        df = df.rename(columns={
            "query1_chr": "ref_chr", "query1_start": "ref_start", "query1_end": "ref_end",
            "query2_chr": "other_chr", "query2_start": "other_start", "query2_end": "other_end",
        })
    elif genome == "q2":
        df = df.rename(columns={
            "query2_chr": "ref_chr", "query2_start": "ref_start", "query2_end": "ref_end",
            "query1_chr": "other_chr", "query1_start": "other_start", "query1_end": "other_end",
        })
    else:
        sys.exit(f"ERROR: genome must be 'q1' or 'q2', got '{genome}'")

    # Strip SUPER_ prefix from chromosome names
    df["ref_chr"] = df["ref_chr"].astype(str).str.replace("SUPER_", "", regex=False)
    df["other_chr"] = df["other_chr"].astype(str).str.replace("SUPER_", "", regex=False)

    for col in ("ref_start", "ref_end", "other_start", "other_end"):
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

    # ensure start <= end within each row
    swap = df["ref_start"] > df["ref_end"]
    df.loc[swap, ["ref_start", "ref_end"]] = (
        df.loc[swap, ["ref_end", "ref_start"]].values
    )

    df = df.sort_values(["ref_chr", "ref_start"]).reset_index(drop=True)

    log.info("Loaded %d conserved blocks from %s (ref = %s)", len(df), path, genome.upper())
    return df


# ---------------------------------------------------------------------------
# Build an interval index for fast block lookup
# ---------------------------------------------------------------------------

def build_block_index(blocks: pd.DataFrame) -> dict:
    """Return {chr: (starts_array, ends_array)} sorted by start."""
    idx = {}
    for chrom, grp in blocks.groupby("ref_chr"):
        starts = grp["ref_start"].to_numpy(dtype=np.int64)
        ends = grp["ref_end"].to_numpy(dtype=np.int64)
        idx[str(chrom)] = (starts, ends)
    return idx


def find_narrowest_gap(block_idx: dict, chrom: str, bp_start: int, bp_end: int):
    """Find the tightest gap between conserved blocks that overlaps the
    breakpoint region [bp_start, bp_end].

    Strategy:
      1. Find all conserved-block gaps that overlap the breakpoint region.
      2. For each overlapping gap, clip it to the breakpoint region.
      3. Return the intersection: the narrowest interval that is within both
         the breakpoint region and the gap between flanking conserved blocks.

    If no gap overlaps the region (breakpoint is entirely inside a conserved
    block), returns None.
    """
    if chrom not in block_idx:
        return None

    starts, ends = block_idx[chrom]
    if len(starts) == 0:
        return None

    # Find blocks that might interact with our region.
    # The first block whose end > bp_start and the last whose start < bp_end.
    i_left = np.searchsorted(ends, bp_start, side="right")  # first block ending after bp_start
    i_right = np.searchsorted(starts, bp_end, side="left")  # first block starting at/after bp_end

    # Collect candidate gaps: between consecutive blocks, plus the edges
    # (before first block, after last block).
    best = None

    # Candidate gap indices: gap_k is the space between block k-1 and block k
    # We also consider gap before block 0 (gap_k = 0 with left edge = 0)
    # and gap after the last block (gap_k = len(starts) with right edge = inf).
    for k in range(max(0, i_left), min(len(starts), i_right) + 1):
        gap_left = ends[k - 1] if k > 0 else 0
        gap_right = starts[k] if k < len(starts) else None

        if gap_right is not None and gap_left >= gap_right:
            continue  # no gap (blocks overlap or are adjacent)

        # Clip gap to breakpoint region
        clipped_start = max(gap_left, bp_start)
        clipped_end = min(gap_right, bp_end) if gap_right is not None else bp_end

        if clipped_start > clipped_end:
            continue  # no overlap

        width = clipped_end - clipped_start
        if best is None or width < best[2]:
            best = (clipped_start, clipped_end, width)

    if best is None:
        return None
    return (best[0], best[1])


# ---------------------------------------------------------------------------
# Refine breakpoints
# ---------------------------------------------------------------------------

def refine_breakpoints(bp_df: pd.DataFrame, block_idx: dict) -> pd.DataFrame:
    """For each breakpoint row in the input BED, narrow the breakpoint region
    using the conserved-block gaps.

    Parameters
    ----------
    bp_df : DataFrame
        Breakpoint BED rows (chr, start, end, bp_class).
    block_idx : dict
        Block index from build_block_index().

    Returns a DataFrame with BED columns:
        chr, start, end, bp_class
    """
    records = []

    for _, row in bp_df.iterrows():
        chrom = str(row["chr"])
        bp_id = str(row["bp_id"])
        bp_class = str(row["bp_class"])

        try:
            bp_start = int(row["start"])
            bp_end = int(row["end"])
        except (ValueError, TypeError):
            continue

        if bp_start > bp_end:
            bp_start, bp_end = bp_end, bp_start

        gap = find_narrowest_gap(block_idx, chrom, bp_start, bp_end)

        if gap is not None:
            refined_start, refined_end = int(gap[0]), int(gap[1])
        else:
            # Breakpoint region is entirely inside a conserved block
            # — keep original coordinates, cannot narrow further
            refined_start, refined_end = bp_start, bp_end

        if refined_start > refined_end:
            refined_start, refined_end = refined_end, refined_start

        records.append({
            "chr": chrom,
            "start": refined_start,
            "end": refined_end,
            "bp_id": bp_id,
            "bp_class": bp_class,
        })

    if not records:
        log.warning("No breakpoints to refine")
        return pd.DataFrame(columns=["chr", "start", "end", "bp_id", "bp_class"])

    out = pd.DataFrame(records)
    out = out.sort_values(["chr", "start"]).reset_index(drop=True)
    return out


# ---------------------------------------------------------------------------
# Write BED
# ---------------------------------------------------------------------------

def merge_and_filter(df: pd.DataFrame, query: str) -> pd.DataFrame:
    """Merge rows sharing chr/start/end/suffix and filter to relevant genome.

    For each bp_class value like ``fusion_Q1`` the suffix is ``Q1``;
    ``shared_fusion`` has no genome suffix and goes to both outputs.

    Steps
    -----
    1. Parse suffix (Q1/Q2/shared) from bp_class.
    2. Keep only rows whose suffix matches *query* or is ``shared``.
    3. Strip the suffix from the type (``fusion_Q1`` → ``fusion``).
    4. Group by chr/start/end and join types with ``;`` (deduped).
    """
    if df.empty:
        return df

    records = []
    for _, row in df.iterrows():
        bp = str(row["bp_class"])

        if bp.startswith("shared_"):
            suffix = "shared"
            base_type = bp                       # keep as-is
        elif bp.endswith("_Q1") or bp.endswith("_Q2"):
            suffix = bp.rsplit("_", 1)[1]        # Q1 or Q2
            base_type = bp.rsplit("_", 1)[0]     # e.g. fusion
        else:
            suffix = "shared"                    # treat unknown as shared
            base_type = bp

        # Filter: keep Q1+shared for q1 output, Q2+shared for q2 output
        if query == "q1" and suffix not in ("Q1", "shared"):
            continue
        if query == "q2" and suffix not in ("Q2", "shared"):
            continue

        records.append({
            "chr": row["chr"],
            "start": int(row["start"]),
            "end": int(row["end"]),
            "bp_id": str(row["bp_id"]),
            "bp_type": base_type,
        })

    if not records:
        return pd.DataFrame(columns=["chr", "start", "end", "bp_id", "bp_class"])

    tmp = pd.DataFrame(records)

    # Group by position, merge types and bp_ids with ;
    merged = (
        tmp.groupby(["chr", "start", "end"], sort=False)
        .agg(
            bp_id=("bp_id", lambda x: ";".join(sorted(set(x)))),
            bp_class=("bp_type", lambda x: ";".join(sorted(set(x)))),
        )
        .reset_index()
    )
    merged = merged.sort_values(["chr", "start"]).reset_index(drop=True)

    log.info("  After merge+filter (%s): %d rows", query.upper(), len(merged))
    log.info("  Breakpoint type counts:\n%s", merged["bp_class"].value_counts().to_string())
    log.info("  Median breakpoint size: %d bp", int(np.median(merged["end"] - merged["start"])))
    log.info("  Min breakpoint size: %d bp", int(np.min(merged["end"] - merged["start"])))
    log.info("  Max breakpoint size: %d bp", int(np.max(merged["end"] - merged["start"])))  
    return merged


def write_bed(df: pd.DataFrame, path: str):
    """Write a BED-like file (tab-separated, no header)."""
    df.to_csv(path, sep="\t", index=False, header=False)
    log.info("Wrote %d breakpoint regions → %s", len(df), path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Refine AGORA breakpoints using cactus conserved alignment blocks."
    )
    parser.add_argument(
        "--input-bed-q1", required=True,
        help="Events BED for Q1 from detect_breakpoints_agora.py",
    )
    parser.add_argument(
        "--input-bed-q2", required=True,
        help="Events BED for Q2 from detect_breakpoints_agora.py",
    )
    parser.add_argument(
        "--blocks-q1", required=True,
        help="conserved_blocks_by_query1.tsv (sorted by Q1 coordinates)",
    )
    parser.add_argument(
        "--blocks-q2", required=True,
        help="conserved_blocks_by_query2.tsv (sorted by Q2 coordinates)",
    )
    parser.add_argument(
        "--bed-q1", default="refined_q1.bed",
        help="Output BED for Q1 genome (default: refined_q1.bed)",
    )
    parser.add_argument(
        "--bed-q2", default="refined_q2.bed",
        help="Output BED for Q2 genome (default: refined_q2.bed)",
    )
    args = parser.parse_args()

    # load conserved blocks & build indices
    blocks_q1 = load_blocks(args.blocks_q1, "q1")
    blocks_q2 = load_blocks(args.blocks_q2, "q2")
    idx_q1 = build_block_index(blocks_q1)
    idx_q2 = build_block_index(blocks_q2)

    # load and refine Q1
    bp_q1 = load_input_bed(args.input_bed_q1)
    refined_q1 = refine_breakpoints(bp_q1, idx_q1)
    q1_anc_path = args.bed_q1.replace(".bed", "_anc_breakpoints.bed")
    write_bed(refined_q1, q1_anc_path)
    refined_q1 = merge_and_filter(refined_q1, "q1")
    write_bed(refined_q1, args.bed_q1)

    # load and refine Q2
    bp_q2 = load_input_bed(args.input_bed_q2)
    refined_q2 = refine_breakpoints(bp_q2, idx_q2)
    q2_anc_path = args.bed_q2.replace(".bed", "_anc_breakpoints.bed")
    write_bed(refined_q2, q2_anc_path)
    refined_q2 = merge_and_filter(refined_q2, "q2")
    write_bed(refined_q2, args.bed_q2)

    log.info("Done.")


if __name__ == "__main__":
    main()

