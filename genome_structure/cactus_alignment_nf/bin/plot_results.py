#!/usr/bin/env python3
"""
Plot conserved blocks and classified breakpoints for two genomes.

Creates a two-panel figure (one per genome) showing:
- Conserved synteny blocks coloured by block_id
- Classified breakpoints (FUSION, FISSION, INVERSION) as vertical markers
  coloured by type.  UNKNOWN breakpoints are excluded.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import sys
import numpy as np
import re

# ── Breakpoint type colours ──────────────────────────────────────────────────
BP_COLORS = {
    "FUSION":    "#e41a1c",   # red
    "FISSION":   "#377eb8",   # blue
    "INVERSION": "#4daf4a",   # green
}

# ── Helpers ──────────────────────────────────────────────────────────────────

def _chr_sort_key(name: str):
    """Natural sort: autosomes first (by number), then sex chromosomes."""
    s = name.replace("SUPER_", "")
    m = re.match(r"^(\d+)$", s)
    if m:
        return (0, int(m.group(1)), "")
    m = re.match(r"^([A-Za-z]+)(\d*)$", s)
    if m:
        return (1, int(m.group(2)) if m.group(2) else 0, m.group(1))
    return (2, 0, s)


# ── Loading functions ────────────────────────────────────────────────────────

def load_conserved_blocks(blocks_dir: str) -> dict[str, pd.DataFrame]:
    """Return {query_key: DataFrame} from conserved_blocks_by_query{1,2}.tsv."""
    blocks = {}
    for f in sorted(Path(blocks_dir).glob("conserved_blocks_by_*.tsv")):
        df = pd.read_csv(f, sep="\t")
        key = f.stem.replace("conserved_blocks_by_", "")   # "query1" / "query2"
        blocks[key] = df
        print(f"  Loaded {f.name}: {len(df)} blocks", file=sys.stderr)
    return blocks


def pick_blocks_df(blocks: dict[str, pd.DataFrame], preferred_key: str) -> pd.DataFrame | None:
    """Return the preferred conserved-block table, or fall back to any available one."""
    if preferred_key in blocks:
        return blocks[preferred_key]

    if not blocks:
        return None

    fallback_key = next(iter(blocks))
    print(
        f"  Missing conserved_blocks_by_{preferred_key}.tsv; using {fallback_key} coordinates instead",
        file=sys.stderr,
    )
    return blocks[fallback_key]


def load_classified_breakpoints(bp_path: str) -> pd.DataFrame:
    """Load one *_breakpoints_classified.tsv, drop UNKNOWNs."""
    df = pd.read_csv(bp_path, sep="\t")
    n_total = len(df)
    df = df[df["type"] != "UNKNOWN"].copy()
    print(f"  Loaded {Path(bp_path).name}: {n_total} total, "
          f"{len(df)} after removing UNKNOWN", file=sys.stderr)
    return df


# ── Plot function ────────────────────────────────────────────────────────────

def plot_genome(ax, blocks_df, chr_col, start_col, end_col,
                bp_df, block_colors, title):
    """
    Draw one genome panel: horizontal bars for conserved blocks,
    vertical lines for classified breakpoints.
    """
    if blocks_df is None or len(blocks_df) == 0:
        ax.text(0.5, 0.5, f"No data", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title(title)
        return

    chromosomes = sorted(blocks_df[chr_col].unique(), key=_chr_sort_key)

    y_pos = 0
    y_labels, y_ticks = [], []

    for chromosome in chromosomes:
        sel = blocks_df[blocks_df[chr_col] == chromosome]
        if sel.empty:
            continue

        y = y_pos
        y_pos += 1
        y_labels.append(chromosome.replace("SUPER_", ""))
        y_ticks.append(y)

        chr_min = sel[[start_col, end_col]].min(axis=1).min()
        chr_max = sel[[start_col, end_col]].max(axis=1).max()

        # Backbone
        ax.plot([chr_min, chr_max], [y, y], "k-", lw=0.8, alpha=0.25)

        # Blocks
        for _, row in sel.iterrows():
            bid = row["block_id"]
            xs = min(row[start_col], row[end_col])
            xe = max(row[start_col], row[end_col])
            color = block_colors.get(bid, "lightgray")
            ax.barh(y, xe - xs, left=xs, height=0.7,
                    color=color, edgecolor="black", linewidth=0.3, alpha=0.8)

        # Breakpoints
        if bp_df is not None and len(bp_df) > 0:
            chr_bps = bp_df[bp_df["chromosome"] == chromosome]
            for _, bp in chr_bps.iterrows():
                bp_mid = (bp["start"] + bp["end"]) / 2
                bp_type = bp["type"]
                colour = BP_COLORS.get(bp_type, "gray")
                ax.plot([bp_mid, bp_mid], [y - 0.42, y + 0.42],
                        color=colour, lw=1.8, alpha=0.85, zorder=5)

    ax.set_ylim(-0.5, y_pos + 0.3)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontsize=9)
    ax.set_xlabel("Position (bp)")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.2, axis="x")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot conserved blocks and classified breakpoints "
                    "for two genomes."
    )
    parser.add_argument("--conserved_blocks", required=True,
                        help="Directory with conserved_blocks_by_query{1,2}.tsv")
    parser.add_argument("--breakpoints_g1", required=True,
                        help="Classified breakpoints TSV for genome 1")
    parser.add_argument("--breakpoints_g2", required=True,
                        help="Classified breakpoints TSV for genome 2")
    parser.add_argument("--name_g1", default="Genome 1",
                        help="Display name for genome 1")
    parser.add_argument("--name_g2", default="Genome 2",
                        help="Display name for genome 2")
    parser.add_argument("--output", default="conserved_blocks_and_breakpoints.pdf",
                        help="Output PDF file")

    args = parser.parse_args()

    # ── Load data ────────────────────────────────────────────────────────
    print("Loading conserved blocks …", file=sys.stderr)
    blocks = load_conserved_blocks(args.conserved_blocks)

    print("Loading classified breakpoints …", file=sys.stderr)
    bp_g1 = load_classified_breakpoints(args.breakpoints_g1)
    bp_g2 = load_classified_breakpoints(args.breakpoints_g2)

    if not blocks:
        sys.exit("No conserved blocks found!")

    # ── Consistent block colours ─────────────────────────────────────────
    all_ids = set()
    for df in blocks.values():
        all_ids.update(df["block_id"].unique())
    n = max(len(all_ids), 3)
    cmap = plt.cm.tab20c(np.linspace(0, 1, n))
    block_colors = {bid: cmap[i % n] for i, bid in enumerate(sorted(all_ids))}

    print(f"Colour map for {len(block_colors)} conserved blocks", file=sys.stderr)

    # ── Figure ───────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(22, 10))

    # Genome 1 – use query1 columns from query1-sorted blocks
    q1_blocks = pick_blocks_df(blocks, "query1")
    plot_genome(axes[0], q1_blocks,
                chr_col="query1_chr", start_col="query1_start",
                end_col="query1_end",
                bp_df=bp_g1, block_colors=block_colors,
                title=f"Conserved Blocks & Breakpoints – {args.name_g1}")

    # Genome 2 – use query2 columns from query2-sorted blocks
    q2_blocks = pick_blocks_df(blocks, "query2")
    plot_genome(axes[1], q2_blocks,
                chr_col="query2_chr", start_col="query2_start",
                end_col="query2_end",
                bp_df=bp_g2, block_colors=block_colors,
                title=f"Conserved Blocks & Breakpoints – {args.name_g2}")

    # ── Legend ────────────────────────────────────────────────────────────
    legend_handles = []
    for bp_type, colour in BP_COLORS.items():
        legend_handles.append(
            mpatches.Patch(facecolor=colour, edgecolor="black",
                           linewidth=0.5, label=bp_type)
        )
    # Add a few representative block colours
    for bid in sorted(block_colors)[:8]:
        legend_handles.append(
            mpatches.Patch(facecolor=block_colors[bid], edgecolor="black",
                           linewidth=0.3, label=bid, alpha=0.8)
        )
    if len(block_colors) > 8:
        legend_handles.append(
            mpatches.Patch(facecolor="white", edgecolor="gray",
                           label=f"… +{len(block_colors) - 8} more blocks")
        )

    fig.legend(handles=legend_handles, loc="lower center", ncol=6,
               fontsize=8, framealpha=0.9, title="Breakpoint types  /  Conserved blocks")

    plt.tight_layout(rect=[0, 0.07, 1, 1])
    plt.savefig(args.output, format="pdf", dpi=300, bbox_inches="tight")
    print(f"\nWrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
