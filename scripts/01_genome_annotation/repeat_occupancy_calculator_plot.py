import pandas as pd
import numpy as np
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc
import argparse
import logging
import os
import re
from datetime import date
from matplotlib.patches import Patch
from matplotlib.colors import to_rgb

# Setup logging
logging.basicConfig(filename="repeat_occupancy.log", level=logging.INFO, 
                    format="%(asctime)s - %(levelname)s - %(message)s")



# Calculate repeat occupancy in genomic windows
def calculate_repeat_occupancy(windows_bedfile, repeats_bedfile):
    """
    Calculate the fraction of each genomic window covered by each repeat family and the total number of overlapping bases.

    Args:
        repeats_bedfile (str): Path to the BED file with repeat positions and classes in the 4th column.
        windows_bedfile (str): Path to the BED file with genomic windows.
        output_csv (str): Path to the output CSV file.
    """
    
    # Load the BED files as BedTool objects
    windows = pybedtools.BedTool(windows_bedfile).sort()
    repeats = pybedtools.BedTool(repeats_bedfile).sort()

    # Intersect repeats with windows, keeping the repeat class (column 4)
    intersection = windows.intersect(repeats, wa=True, wb=True)

    # Parse the intersection results
    data = []
    for line in intersection:
        window_chrom, window_start, window_end, repeat_chrom, repeat_start, repeat_end, repeat_class = (
            line[0], int(line[1]), int(line[2]), line[3], int(line[4]), int(line[5]), line[6]
        )
        overlap = min(window_end, repeat_end) - max(repeat_start, window_start)
        window_length = window_end - window_start
        fraction_covered = overlap / window_length
        # Extract the repeat family (everything before the '/')
        # If the class contains a comma, rename as 'overlapping'
        if ',' in repeat_class:
            repeat_family = 'Overlapping'
        else:
            repeat_family = repeat_class.split('/')[0]
        data.append((window_chrom, window_start, window_end, repeat_family, overlap, fraction_covered))

    # Convert to a DataFrame
    df = pd.DataFrame(data, columns=["chrom", "start", "end", "repeat_family", "overlap_bases", "fraction_covered"])

    # Group by window and repeat family, summing the overlaps and fractions
    grouped = df.groupby(["chrom", "start", "end", "repeat_family"]).agg(
        {"overlap_bases": "sum", "fraction_covered": "sum"}
    ).reset_index()

    # Save the grouped data in long format to CSV
    return grouped


def summarise_whole_genome_repeat_coverage(repeats_bedfile, genome_fai, output_summary_file):
    """
    Summarise genome-wide repeat coverage by repeat family and total repeats.

    Returns a DataFrame with columns:
      - repeat_family
      - repeat_count
      - covered_bases
      - genome_length
      - fraction_genome_covered

    The final row, 'Total_repeats', gives the total number of repeat records
    across all families and the combined fraction of the genome covered by any repeat.
    """
    genome_df = pd.read_csv(
        genome_fai,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["chrom", "length"],
    )
    total_genome_length = int(genome_df["length"].sum())

    bed = pybedtools.BedTool(repeats_bedfile).sort()
    records = []
    for line in bed:
        chrom, start, end, repeat_class = line[0], int(line[1]), int(line[2]), str(line[3])
        if "," in repeat_class:
            repeat_family = "Overlapping"
        else:
            repeat_family = repeat_class.split("/")[0]
        records.append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "repeat_family": repeat_family,
            "length": end - start,
        })

    summary_df = pd.DataFrame(records)
    if summary_df.empty:
        summary_df = pd.DataFrame([
            {"repeat_family": "Total_repeats", "repeat_count": 0, "covered_bases": 0, "genome_length": total_genome_length, "fraction_genome_covered": 0.0}
        ])
        summary_df.to_csv(output_summary_file, sep="\t", index=False)
        return summary_df

    family_summary = (
        summary_df.groupby("repeat_family", as_index=False)
        .agg(
            repeat_count=("repeat_family", "count"),
            covered_bases=("length", "sum"),
        )
    )
    family_summary["genome_length"] = total_genome_length
    family_summary["fraction_genome_covered"] = family_summary["covered_bases"] / total_genome_length
    family_summary = family_summary.sort_values(
        ["fraction_genome_covered", "repeat_family"],
        ascending=[False, True],
        kind="mergesort",
    ).reset_index(drop=True)

    merged_bed = bed.merge()
    total_covered_bases = sum(int(record[2]) - int(record[1]) for record in merged_bed)
    total_count = len(summary_df)
    total_row = pd.DataFrame([
        {
            "repeat_family": "Total_repeats",
            "repeat_count": total_count,
            "covered_bases": total_covered_bases,
            "genome_length": total_genome_length,
            "fraction_genome_covered": total_covered_bases / total_genome_length,
        }
    ])
    summary_df = pd.concat([family_summary, total_row], ignore_index=True)

    output_dir = os.path.dirname(output_summary_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    summary_df.to_csv(output_summary_file, sep="\t", index=False)

    return summary_df


def get_family_colours(data, family_colours_file):
    """
    Build a colour map for repeat families using a predefined TSV file.

    The file must be tab-separated with:
      column 1: repeat family
    column 2: hex colour code

    Families in the colour file but not present in the current dataset are ignored.
    Families present in the data but missing from the colour file are assigned fallback
    colours from glasbey to keep plotting robust.
    """
    unique_families = sorted(data["repeat_family"].dropna().unique())

    colours_df = pd.read_csv(
        family_colours_file,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["repeat_family", "hex_colour"],
        dtype={"repeat_family": str, "hex_colour": str}
    )
    colours_df = colours_df.dropna(subset=["repeat_family", "hex_colour"])

    predefined_colours = {}
    for _, row in colours_df.iterrows():
        family = row["repeat_family"].strip()
        hex_colour = row["hex_colour"].strip().strip('"').strip("'")
        if hex_colour and not hex_colour.startswith("#") and len(hex_colour) in (3, 6):
            hex_colour = f"#{hex_colour}"
        if family in unique_families:
            try:
                predefined_colours[family] = to_rgb(hex_colour)
            except ValueError:
                logging.warning(f"Invalid colour '{hex_colour}' for family '{family}'. Skipping.")

    missing_families = [f for f in unique_families if f not in predefined_colours]
    if missing_families:
        logging.warning(
            "Missing predefined colours for families: %s. Assigning fallback colours.",
            ", ".join(missing_families)
        )
        # Use a high-contrast palette for fallbacks to avoid pale/white-looking bars.
        fallback_palette = sns.color_palette(cc.glasbey, n_colors=len(missing_families))
        for family, colour in zip(missing_families, fallback_palette):
            predefined_colours[family] = colour

    return predefined_colours


def chromosome_sort_key(chrom_name):
    """
    Sort chromosomes by the chromosome label after the last underscore.

    Examples:
      taxa1_1_1  -> 1
      taxa1_1_W1 -> W1
      taxa1_1_Z  -> Z
      SUPER_2    -> 2

    Ordering:
      1) Pure numeric labels in numeric order (1, 2, ... 10)
      2) Sex-like labels X, Y, W, Z (and numbered forms like W1) in natural order
      3) Other alphanumeric labels in natural order
    """
    chrom_str = str(chrom_name)
    chr_label = chrom_str.rsplit("_", 1)[-1]
    chr_label_upper = chr_label.upper()

    # 1) Pure numeric labels
    if re.fullmatch(r"\d+", chr_label):
        return (0, int(chr_label), "", 0, chr_label_upper, chrom_str)

    # 2) Sex-like labels (X, Y, W, Z) with optional trailing digits (e.g., W1)
    sex_match = re.fullmatch(r"([WXYZ])(\d*)", chr_label_upper)
    if sex_match:
        base = sex_match.group(1)
        number = int(sex_match.group(2)) if sex_match.group(2) else 0
        sex_rank = {"W": 0, "Z": 1, "X": 2, "Y": 3}
        return (1, sex_rank[base], base, number, chr_label_upper, chrom_str)

    # 3) Other labels: natural token sort (A2 < A10)
    tokens = re.split(r"(\d+)", chr_label_upper)
    token_key = tuple(int(t) if t.isdigit() else t for t in tokens)
    return (2, token_key, chr_label_upper, chrom_str)


def make_family_legend_handles(family_order, family_colours):
    return [
        Patch(facecolor=family_colours[family], edgecolor="none", label=family)
        for family in family_order
        if family in family_colours
    ]


def get_family_order(data):
    return (
        data.groupby("repeat_family", sort=False)["fraction_covered"]
        .median()
        .sort_values(ascending=False)
        .index
        .tolist()
    )


# Plotting repeat occupancy
def box_plot(result_df, output_summary_plot, sample_name, family_colours_file):
    family_colours = get_family_colours(result_df, family_colours_file)
    sample_name = sample_name

    # Compute median (or mean) for sorting
    median_vals = result_df.groupby("repeat_family")["fraction_covered"].median()
    sorted_families = median_vals.sort_values(ascending=False).index.tolist()

    plt.figure(figsize=(12, 6))
    sns.violinplot(data=result_df, x="repeat_family", y="fraction_covered", hue="repeat_family",
                   palette=family_colours, dodge=False, cut=0, order=sorted_families)  # Use violin plot with consistent colours
    plt.title(f"Repeat Occupancy per Class in {sample_name}")
    plt.xticks(rotation=45)
    plt.ylabel("Fraction Covered")
    plt.tight_layout()
    plt.savefig(output_summary_plot)

# Plotting repeat occupancy
def box_plot_log(result_df, output_summary_plot, sample_name, family_colours_file):
    family_colours = get_family_colours(result_df, family_colours_file)
    sample_name = sample_name

    # Compute median (or mean) for sorting
    median_vals = result_df.groupby("repeat_family")["fraction_covered"].median()
    sorted_families = median_vals.sort_values(ascending=False).index.tolist()

    plt.figure(figsize=(12, 6))
    sns.violinplot(data=result_df, x="repeat_family", y="fraction_covered", hue="repeat_family", 
                   palette=family_colours, dodge=False, cut=0, order=sorted_families)  # Use violin plot with consistent colours
    plt.yscale('log')
    plt.title(f"Repeat Occupancy per Class in {sample_name}")
    plt.xticks(rotation=45)
    plt.ylabel("Fraction Covered (log scale)")    
    plt.tight_layout()
    plt.savefig(output_summary_plot)


# Define a plotting function for stacked bars
def stacked_barplot(data, family_colours, family_order, bar_width=0.8, **kwargs):
    # Ensure data is sorted and grouped by start positions
    unique_starts = data["start"].unique()
    unique_starts.sort()
    # Initialize bottom_values as a float dtype to avoid dtype incompatibility
    bottom_values = pd.Series([0.0] * len(unique_starts), index=unique_starts)

    for cls in family_order:
        if cls not in data["repeat_family"].values:
            continue
        cls_data = data[data["repeat_family"] == cls].sort_values("start")
        cls_data = cls_data.set_index("start").reindex(unique_starts, fill_value=0).reset_index()
        # Remove 'color' from kwargs to avoid duplication
        kwargs.pop("color", None)
        plt.bar(cls_data["start"], cls_data["fraction_covered"], 
                bottom=bottom_values.loc[cls_data["start"]].values, 
            label=cls,
            color=family_colours[cls],
            alpha=1.0,
            edgecolor="none",
            linewidth=0,
            width=bar_width,
            antialiased=False,
            **kwargs,
        )
        bottom_values.loc[cls_data["start"]] += cls_data["fraction_covered"].values

    # Fix the y-axis
    plt.ylim(0, 1)  # Set y-axis limits (e.g., 0 to 1 for fractions)
    plt.ylabel("Fraction Covered")  # Add y-axis label
    plt.yticks(np.linspace(0, 1, 6))  # Set y-axis ticks (e.g., 0, 0.2, 0.4, ..., 1)

# Plot using seaborn with facet wrap
def plot_repeat_occupancy(coverage_df, output_file, family_colours_file, family_order):
    coverage_df = coverage_df[~coverage_df["chrom"].str.contains("SCAFF|unloc", case=False, na=False)]
    # Cast to float before dividing to avoid dtype warnings
    coverage_df = coverage_df.copy()
    coverage_df.loc[:, "start"] = coverage_df["start"].astype(float) / 1_000_000

    # Estimate bar width from genomic window spacing (in Mb) so bars are visible in facet plots.
    unique_starts = np.sort(coverage_df["start"].dropna().unique())
    if len(unique_starts) > 1:
        spacing = np.diff(unique_starts)
        spacing = spacing[spacing > 0]
        if len(spacing) > 0:
                # Use a slight overlap so adjacent bars render without visible seams.
                bar_width = float(np.median(spacing)) * 1.01
        else:
            bar_width = 0.1
    else:
        bar_width = 0.1

    # Generate a consistent colour mapping for repeat families
    family_colours = get_family_colours(coverage_df, family_colours_file)

    # Plot using seaborn with separate plots per chromosome
    logging.info("Generating multi-chromosome plots")

    chromosomes = sorted(coverage_df["chrom"].unique(), key=chromosome_sort_key)
    n_chromosomes = len(chromosomes)
    cols = min(4, n_chromosomes)
    rows = int(np.ceil(n_chromosomes / cols))

    plt.figure(figsize=(cols * 4, rows * 3))
    sns.set(style="whitegrid")

    # Use FacetGrid for clear multi-chromosome view
    g = sns.FacetGrid(
        coverage_df,
        col="chrom",
        col_order=chromosomes,
        col_wrap=cols,
        sharey=False,
        height=3,
        aspect=1.5,
    )

    # Apply the plotting function to each facet
    g.map_dataframe(
        stacked_barplot,
        family_colours=family_colours,
        family_order=family_order,
        bar_width=bar_width,
    )
    g.set_titles("Chromosome: {col_name}")
    g.set_axis_labels("Genomic Position (Mb)", "Repeat Occupancy")
    legend_handles = make_family_legend_handles(family_order, family_colours)
    g.fig.legend(
        handles=legend_handles,
        labels=[handle.get_label() for handle in legend_handles],
        title="Repeat Family",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )
    g.tight_layout()
    g.savefig(output_file, dpi=300)

    logging.info("Multi-chromosome facet plot saved.")

# Plot repeat landscape using divergence bins as a stacked bar plot
def plot_repeat_landscape(repeats_file, output_landscape_plot, genome_file, family_colours_file, family_order, bin_size=1):
    """
    Plot the repeat landscape using divergence (Kimura distances) as a stacked bar plot.

    Args:
        repeats_file (str): Path to the BED file with repeat positions and divergence values (e.g., in the 5th column).
        output_landscape_plot (str): Path to save the repeat landscape plot.
        genome_file (str): Path to the genome file (FAI format) to calculate total genome length.
        bin_size (int): Size of the bins for divergence (default: 1%).
    """
    logging.info(f"Generating repeat landscape stacked bar plot from {repeats_file}")

    # Load the repeats file
    repeats_df = pd.read_csv(repeats_file, sep="\t", header=None, names=["chrom", "start", "end", "class", "divergence"])
    
    # Extract the repeat family (everything before the '/')
    repeats_df["family"] = repeats_df["class"].str.split('/').str[0]

    # Rename 'family' to 'repeat_family' for consistency with other plots
    repeats_df.rename(columns={"family": "repeat_family"}, inplace=True)

    # Load the genome file to calculate total genome length
    genome_df = pd.read_csv(genome_file, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"])
    total_genome_length = genome_df["length"].sum()

    # Bin the divergence values
    repeats_df["divergence_bin"] = (repeats_df["divergence"] // bin_size) * bin_size

    # Group by divergence bin and repeat family, summing the lengths
    repeats_df["length"] = repeats_df["end"] - repeats_df["start"]
    grouped = repeats_df.groupby(["divergence_bin", "repeat_family"])["length"].sum().reset_index()

    # Calculate the percentage of the genome length
    grouped["percentage"] = (grouped["length"] / total_genome_length) * 100

    # Pivot the data for stacked bar plotting
    pivoted = grouped.pivot(index="divergence_bin", columns="repeat_family", values="percentage").fillna(0)

    # Use the same family order as the occupancy plot, then append any landscape-only families.
    landscape_family_order = [family for family in family_order if family in pivoted.columns]
    landscape_family_order.extend([family for family in pivoted.columns if family not in landscape_family_order])

    # Generate a consistent colour mapping for repeat families
    family_colours = get_family_colours(grouped, family_colours_file)
    sample_name = os.path.basename(genome_file).split('.')[0]

    # Plot the stacked bar plot
    plt.figure(figsize=(12, 6))
    bottom_values = np.zeros(len(pivoted))
    for family in landscape_family_order:
        plt.bar(pivoted.index, pivoted[family], bottom=bottom_values, label=family, color=family_colours[family])
        bottom_values += pivoted[family]

    legend_handles = make_family_legend_handles(landscape_family_order, family_colours)

    plt.title(f"Repeat Landscape for {sample_name} (Divergence vs. Percentage of Genome)")
    plt.xlabel("Divergence (%)")
    plt.ylabel("Percentage of Genome")
    plt.legend(handles=legend_handles, title="Repeat Family", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_landscape_plot)
    logging.info(f"Repeat landscape stacked bar plot saved to {output_landscape_plot}")

# Command-line interface (CLI)
def main():
    parser = argparse.ArgumentParser(description="Calculate repeat occupancy in genomic windows.")
    parser.add_argument("--genome-fai", required=True, help="Path to the genome FAI file (e.g., genome.fai)")
    parser.add_argument("--prefix", required=True, help="Prefix for input and output files")
    parser.add_argument("--family-colours", required=True, help="Path to TSV with repeat family in column 1 and hex colour in column 2")
    parser.add_argument("--window", type=int, default=100000, help="Window size (default: 100kb)")
    parser.add_argument("--bin-size", type=int, default=1, help="Bin size for divergence (default: 1%)")

    args = parser.parse_args()
    logging.info("Validateing input files and parameters")
    # Validate input files
    genome_file = args.genome_fai
    if not os.path.exists(genome_file):
        logging.error(f"Genome file {genome_file} does not exist.")
        raise FileNotFoundError(f"Genome file {genome_file} does not exist.")   

    family_colours_file = args.family_colours
    if not os.path.exists(family_colours_file):
        logging.error(f"Colour mapping file {family_colours_file} does not exist.")
        raise FileNotFoundError(f"Colour mapping file {family_colours_file} does not exist.")
    
    repeats_file = f"{args.prefix}_repeats_merged.bed"
    if not os.path.exists(repeats_file):
        logging.error(f"Repeats file {repeats_file} does not exist.")
        raise FileNotFoundError(f"Repeats file {repeats_file} does not exist.") 
    
    repeats_file_unmerged = f"{args.prefix}_repeats.bed"
    if not os.path.exists(repeats_file_unmerged):
        logging.error(f"Repeats file {repeats_file_unmerged} does not exist.")
        raise FileNotFoundError(f"Repeats file {repeats_file_unmerged} does not exist.") 

    windows_df = f"{args.prefix}_{args.window}_wd.bed"
    if not os.path.exists(windows_df):
        logging.error(f"Windows file {windows_df} does not exist.")
        raise FileNotFoundError(f"Windows file {windows_df} does not exist.")   
    
    # Prepare output file names
    output_csv = f"{args.prefix}_{args.window}_repeat_occupancy.csv"
    today = date.today().strftime("%Y%m%d")
    output_summary_tsv = f"{args.prefix}_repeat_family_genome_summary_{today}.tsv"
    output_plot = f"{args.prefix}_{args.window}_repeat_occupancy_plot.png"
    output_summary_plot = f"{args.prefix}_{args.window}_repeat_summary_plot.png"
    output_summary_plot_log = f"{args.prefix}_{args.window}_repeat_summary_plot_log.png"
    output_landscape_plot = f"{args.prefix}_{args.window}_repeat_landscape_plot.png"
    logging.info("Starting repeat occupancy calculation")
    coverage_df = calculate_repeat_occupancy(windows_df, repeats_file)
    family_order = get_family_order(coverage_df)

    # Save CSV output
    logging.info(f"Saving repeat occupancy to CSV: {output_csv}")
    coverage_df.to_csv(output_csv, index=False)

    logging.info(f"Calculating whole-genome repeat family summary: {output_summary_tsv}")
    summarise_whole_genome_repeat_coverage(repeats_file, genome_file, output_summary_tsv)

    sample_name = os.path.basename(args.prefix).split('.')[0]
    # Plot
    plot_repeat_occupancy(coverage_df, output_plot, family_colours_file, family_order)
    box_plot(coverage_df, output_summary_plot, sample_name, family_colours_file)
    box_plot_log(coverage_df, output_summary_plot_log, sample_name, family_colours_file)
    plot_repeat_landscape(repeats_file_unmerged, output_landscape_plot, genome_file, family_colours_file, family_order, args.bin_size)
    logging.info("Repeat occupancy calculation completed")

if __name__ == "__main__":
    main()
