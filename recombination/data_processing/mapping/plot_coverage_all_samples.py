#!/usr/bin/env python3
import argparse
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_fai(fai_path):
    fai = pd.read_csv(
        fai_path,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["chrom", "length"],
    )
    fai = fai[fai["length"] > 0].copy()
    fai["offset"] = fai["length"].cumsum() - fai["length"]
    fai["midpoint"] = fai["offset"] + fai["length"] / 2.0
    return fai


def load_sample_coverage(path, chrom_offsets):
    sample_name = os.path.basename(path).replace(".regions.bed.gz", "")
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "depth"],
        compression="gzip",
    )
    df = df[df["chrom"].isin(chrom_offsets.keys())].copy()
    df["sample"] = sample_name
    df["mid"] = (df["start"] + df["end"]) / 2.0 + df["chrom"].map(chrom_offsets)
    return df[["chrom", "start", "end", "depth", "sample", "mid"]]


def build_combined_table(sample_tables):
    combined = pd.concat(sample_tables, ignore_index=True)
    summary = (
        combined.groupby(["chrom", "start", "end"], as_index=False)["depth"]
        .median()
        .rename(columns={"depth": "median_depth"})
    )
    return combined, summary


def build_individual_summary(combined):
    autosomal_mask = ~combined["chrom"].str.contains(
        r"(?:^|[^A-Za-z])(?:Z|W)(?:$|[^A-Za-z])",
        case=False,
        regex=True,
    )
    autosomal_data = combined[autosomal_mask]

    individual_summary = (
        combined.groupby("sample", as_index=False)
        .agg(
            windows=("depth", "size"),
            median_depth=("depth", "median"),
            mean_depth=("depth", "mean"),
        )
        .sort_values("median_depth", ascending=False)
    )

    if autosomal_data.empty:
        individual_summary["autosomal_median_depth"] = np.nan
    else:
        autosomal_summary = (
            autosomal_data.groupby("sample", as_index=False)["depth"]
            .median()
            .rename(columns={"depth": "autosomal_median_depth"})
        )
        individual_summary = individual_summary.merge(autosomal_summary, on="sample", how="left")

    return individual_summary


def read_coverage_list(list_path):
    with open(list_path, "r", encoding="utf-8") as handle:
        entries = [line.strip() for line in handle if line.strip() and not line.lstrip().startswith("#")]

    if not entries:
        raise FileNotFoundError(f"No coverage files listed in {list_path}")

    missing = [path for path in entries if not os.path.isfile(path)]
    if missing:
        missing_preview = ", ".join(missing[:5])
        if len(missing) > 5:
            missing_preview += ", ..."
        raise FileNotFoundError(f"Listed coverage files not found: {missing_preview}")

    return sorted(entries)


def build_sample_label_map(sample_names):
    # Remove a common trailing token suffix (split by . _ -) to keep labels concise.
    token_lists = [re.split(r"[._-]+", name) for name in sample_names]

    common_suffix_len = 0
    while True:
        idx_from_end = common_suffix_len + 1
        if any(len(tokens) < idx_from_end for tokens in token_lists):
            break

        trailing_tokens = [tokens[-idx_from_end] for tokens in token_lists]
        if len(set(trailing_tokens)) == 1:
            common_suffix_len += 1
        else:
            break

    raw_labels = []
    for name, tokens in zip(sample_names, token_lists):
        if common_suffix_len > 0 and len(tokens) > common_suffix_len:
            candidate = "_".join(tokens[:-common_suffix_len])
            raw_labels.append(candidate if candidate else name)
        else:
            raw_labels.append(name)

    # Keep labels unique after pruning.
    label_counts = {}
    unique_labels = []
    for label in raw_labels:
        label_counts[label] = label_counts.get(label, 0) + 1
        if label_counts[label] == 1:
            unique_labels.append(label)
        else:
            unique_labels.append(f"{label}_{label_counts[label]}")

    return dict(zip(sample_names, unique_labels))


def plot_coverage(
    combined,
    summary,
    fai,
    output_png,
    output_tsv,
    title,
    prune_sample_names,
):
    # Only plot chromosomes that are present in the coverage data
    present_chroms = set(combined["chrom"].unique())
    chrom_order = [c for c in fai["chrom"].tolist() if c in present_chroms]
    
    # Debug: show what chromosomes are included/excluded
    fai_chroms = set(fai["chrom"].tolist())
    missing_from_fai = present_chroms - fai_chroms
    missing_from_coverage = fai_chroms - present_chroms
    if missing_from_fai:
        print(f"Warning: chromosomes in coverage but NOT in FAI (excluded from plot): {sorted(missing_from_fai)}")
    if missing_from_coverage:
        print(f"Note: chromosomes in FAI but NOT in coverage data: {sorted(missing_from_coverage)}")
    print(f"Plotting {len(chrom_order)} chromosomes: {chrom_order}")

    # Rank samples by overall median coverage so rows are ordered by depth.
    ordering_medians = combined.groupby("sample")["depth"].median().sort_values(ascending=False)
    sample_order = ordering_medians.index.tolist()
    if prune_sample_names:
        sample_label_map = build_sample_label_map(sample_order)
    else:
        sample_label_map = {name: name for name in sample_order}

    sample_label_order = [sample_label_map[name] for name in sample_order]
    n_chroms = len(chrom_order)
    n_cols = min(4, n_chroms)
    n_rows = (n_chroms + n_cols - 1) // n_cols

    fig_width = n_cols * 7
    fig_height = max(len(sample_order) * 0.22 * n_rows, n_rows * 3)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(fig_width, fig_height),
        squeeze=False,
    )
    fig.suptitle(title, fontsize=13, y=1.01)

    # Normalise each sample by its own autosomal median depth so that sex
    # chromosome ploidy differences do not affect the baseline.
    autosomal_mask = ~combined["chrom"].str.contains(r"(?:^|[^A-Za-z])(?:Z|W)(?:$|[^A-Za-z])", case=False, regex=True)
    autosomal_data = combined[autosomal_mask]

    if autosomal_data.empty:
        print("Warning: no autosomal windows detected (Z/W filter). Falling back to all chromosomes for normalization.")
        norm_medians = combined.groupby("sample")["depth"].median()
    else:
        norm_medians = autosomal_data.groupby("sample")["depth"].median()
        norm_medians = norm_medians.reindex(sample_order)
        fallback_medians = combined.groupby("sample")["depth"].median()
        norm_medians = norm_medians.fillna(fallback_medians)

    combined = combined.copy()
    combined["norm_depth"] = combined["depth"] / combined["sample"].map(norm_medians)
    combined["sample_label"] = combined["sample"].map(sample_label_map)

    # Shared colour scale across all panels: cap at 2× normalised depth
    vmin, vmax = 0, 2

    for panel_idx, chrom in enumerate(chrom_order):
        row, col = divmod(panel_idx, n_cols)
        ax = axes[row][col]

        chrom_data = combined[combined["chrom"] == chrom].copy()
        chrom_data["pos_mb"] = (chrom_data["start"] + chrom_data["end"]) / 2e6

        # Build a 2-D matrix: rows = samples, columns = windows (sorted by position)
        pivot = (
            chrom_data
            .pivot_table(index="sample_label", columns="pos_mb", values="norm_depth", aggfunc="mean")
            .reindex(index=sample_label_order)
        )

        im = ax.imshow(
            pivot.values,
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
            cmap="RdYlBu_r",
            interpolation="nearest",
            extent=[pivot.columns.min(), pivot.columns.max(), len(sample_label_order) - 0.5, -0.5],
        )

        ax.set_title(chrom, fontsize=9)
        ax.set_xlabel("Position (Mb)", fontsize=7)
        ax.set_yticks(range(len(sample_label_order)))
        ax.set_yticklabels(sample_label_order, fontsize=5)
        ax.tick_params(axis="x", labelsize=6)

        # Colourbar per panel
        cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
        cbar.set_label("Norm. depth", fontsize=6)
        cbar.ax.tick_params(labelsize=5)

    # Hide unused panels
    for unused in range(n_chroms, n_rows * n_cols):
        row, col = divmod(unused, n_cols)
        axes[row][col].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close()

    summary.to_csv(output_tsv, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Overlay coverage from many mosdepth samples on one reference genome axis."
    )
    coverage_input = parser.add_mutually_exclusive_group(required=True)
    coverage_input.add_argument(
        "--coverage-dir",
        help="Directory containing *.regions.bed.gz files from mosdepth.",
    )
    coverage_input.add_argument(
        "--coverage-list",
        help="Text file with one .regions.bed.gz path per line (comments with # supported).",
    )
    parser.add_argument(
        "--fai",
        required=True,
        help="Reference FASTA index (.fai) used to define chromosome order and lengths.",
    )
    parser.add_argument(
        "--output-png",
        required=True,
        help="Path to save the combined coverage plot.",
    )
    parser.add_argument(
        "--output-tsv",
        required=True,
        help="Path to save median coverage per window across samples.",
    )
    parser.add_argument(
        "--output-individual-tsv",
        default=None,
        help="Optional path to save per-individual median coverage summary.",
    )
    parser.add_argument(
        "--title",
        default="Coverage across reference genome",
        help="Figure title.",
    )
    parser.add_argument(
        "--keep-full-sample-names",
        action="store_true",
        help="Disable automatic pruning of shared trailing sample-name suffixes.",
    )

    args = parser.parse_args()

    if args.coverage_dir:
        coverage_files = sorted(glob.glob(os.path.join(args.coverage_dir, "*.regions.bed.gz")))
        if not coverage_files:
            raise FileNotFoundError(f"No .regions.bed.gz files found in {args.coverage_dir}")
    else:
        coverage_files = read_coverage_list(args.coverage_list)

    fai = read_fai(args.fai)
    chrom_offsets = dict(zip(fai["chrom"], fai["offset"]))

    sample_tables = [load_sample_coverage(path, chrom_offsets) for path in coverage_files]
    combined, summary = build_combined_table(sample_tables)
    individual_summary = build_individual_summary(combined)

    if args.output_individual_tsv:
        output_individual_tsv = args.output_individual_tsv
    else:
        output_root, output_ext = os.path.splitext(args.output_tsv)
        if output_ext.lower() == ".tsv":
            output_individual_tsv = f"{output_root}.per_individual.tsv"
        else:
            output_individual_tsv = f"{args.output_tsv}.per_individual.tsv"

    plot_coverage(
        combined=combined,
        summary=summary,
        fai=fai,
        output_png=args.output_png,
        output_tsv=args.output_tsv,
        title=args.title,
        prune_sample_names=not args.keep_full_sample_names,
    )
    individual_summary.to_csv(output_individual_tsv, sep="\t", index=False)

    print(f"Loaded {len(coverage_files)} samples")
    print(f"Plot written to: {args.output_png}")
    print(f"Median table written to: {args.output_tsv}")
    print(f"Per-individual summary written to: {output_individual_tsv}")


if __name__ == "__main__":
    main()
