#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def zeropad_chrom(s, width=2):
    try:
        # Case: CAR_1, chr_12, etc.
        if '_' in s:
            prefix, num = s.split('_')
            return f"{prefix}_{int(num):0{width}d}"

        # Case: pure numeric like "1", "12"
        if s.isdigit():
            return f"{int(s):0{width}d}"

        # Case: letters like "Z", "X"
        return s

    except Exception:
        return s
    

# -----------------------------
# Breakpoint classification
# -----------------------------
def simplify_bp_class(bp_class_str):
    if not isinstance(bp_class_str, str):
        return None, set()

    if not bp_class_str or bp_class_str == 'NA':
        return None, set()

    parts = {x.strip() for x in bp_class_str.split(';') if x}

    classes = set()
    genomes = set()
    has_shared = False

    for val in parts:
        if val.startswith('fusion_'):
            classes.add('fusion')
            genomes.add(val.split('_')[1])

        elif val.startswith('fission_'):
            classes.add('fission')
            genomes.add(val.split('_')[1])

        elif val.startswith('inversion_'):
            classes.add('inversion')
            genomes.add(val.split('_')[1])

        elif val == 'shared_inversion':
            classes.add('shared_inversion')
            has_shared = True

        elif val == 'shared_fission':
            classes.add('shared_fission')
            has_shared = True

        elif val == 'shared_fusion':
            classes.add('shared_fusion')
            has_shared = True

        else:
            classes.add(val)

    # determine class
    if len(classes) == 1:
        final_class = next(iter(classes))
    else:
        final_class = 'complex'

    # determine genomes to plot on
    if has_shared:
        final_genomes = {'Q1', 'Q2'}
    elif genomes:
        final_genomes = genomes
    else:
        final_genomes = set()

    return final_class, final_genomes

# -----------------------------
# Plot one genome
# -----------------------------
def plot_genome(ax, df, chr_col, start_col, end_col, title, color_map):

    df = df[
    df[chr_col].notna() &
    (df[chr_col].astype(str).str.lower() != 'nan') &
    (df[chr_col] != '.') &
    (df[chr_col] != '')
    ].copy()

    chromosomes = sorted(df[chr_col].dropna().unique())
    genome_name = chr_col.split('_')[0].upper()  # e.g. 'q1_chr' → 'Q1'
    
    if genome_name == 'ANC':
        title = title # override title for ancestor
    else:
        title = df[title].dropna().iloc[0]

    y_ticks = []
    y_labels = []

    for i, chrom in enumerate(chromosomes):
        sub = df[df[chr_col] == chrom].sort_values(start_col)

        y = i
        y_ticks.append(y)
        y_labels.append(chrom)

        # backbone
        ax.barh(y, sub[end_col].max(), left=sub[start_col].min(), height=0.7, 
            color='white', edgecolor='black', linewidth=0.5, alpha=0.8)

        # genes + breakpoints
        for _, row in sub.iterrows():
            start = row[start_col]
            end = row[end_col]
            anc = row['anc_chr']

            # gene painting
            color = color_map.get(anc, 'lightgray')
            ax.barh(y, end - start, left=start, height=0.8,
                    color=color, edgecolor='none')

            # breakpoint overlay
            bp_type, bp_genome = simplify_bp_class(row.get('bp_class', 'NA'))

            if not bp_type:
                continue

            pos = (start + end) / 2
            if bp_type:

                print(bp_type, bp_genome, genome_name)

                draw = False
                style = {}

                # shared always on both
                if bp_type.startswith('shared'):
                    if genome_name in bp_genome:
                        draw = True
                        style = dict(color='grey', linewidth=3)

                # complex: respect genome assignment
                elif bp_type == 'complex':
                    if genome_name in bp_genome:
                        draw = True
                        style = dict(color='black', linewidth=3)

                elif bp_type == 'fusion':
                    if genome_name in bp_genome:
                        draw = True
                        style = dict(color='brown', linewidth=3)

                elif bp_type == 'fission':
                    if genome_name in bp_genome:
                        draw = True
                        style = dict(color='goldenrod', linewidth=3)

                elif bp_type == 'inversion' or bp_type == 'shared_inversion':
                    if genome_name in bp_genome:
                        draw = True
                        style = dict(color='purple', linewidth=3)

                if draw:
                    ax.plot([pos, pos], [y - 0.4, y + 0.4], **style)


    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel("Genomic position")
    ax.grid(axis='x', alpha=0.3)


# -----------------------------
# Main
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Chromosome painting by ancestor with breakpoint overlay"
    )
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', default='chromosome_painting.pdf')

    args = parser.parse_args()

    df = pd.read_csv(
        args.input,
        sep='\t',
        header=0,
        engine='python'
    )

    df['anc_chr'] = df['anc_chr'].astype(str).apply(zeropad_chrom)
    df['q1_chr'] = df['q1_chr'].astype(str).apply(zeropad_chrom)
    df['q2_chr'] = df['q2_chr'].astype(str).apply(zeropad_chrom)

    print(df.head())
    print(df.columns)
    print(df['anc_chr'].unique())
    print(df['q1_chr'].unique())
    print(df['q2_chr'].unique())

    # ensure numeric
    for col in ['anc_start','anc_end','q1_start','q1_end','q2_start','q2_end']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df = df.dropna(subset=['anc_chr'])

    # -----------------------------
    # Color map (ancestor-based)
    # -----------------------------
    anc_chrs = sorted(df['anc_chr'].unique())

    tab20 = plt.colormaps.get_cmap('tab20')
    extra_cmap = plt.colormaps.get_cmap('Set3')

    n = len(anc_chrs)

    colors = []

    for i in range(n):
        if i < 20:
            colors.append(tab20(i / 20))
        else:
            colors.append(extra_cmap((i - 20) / (n - 20)))
            
    colors = np.array(colors)

    color_map = dict(zip(anc_chrs, colors))

    # -----------------------------
    # Plot
    # -----------------------------
    fig, axes = plt.subplots(3, 1, figsize=(16, 12))

    plot_genome(
        axes[0], df,
        'anc_chr', 'anc_start', 'anc_end',
        'Ancestral genome', color_map
    )

    plot_genome(
        axes[1], df,
        'q1_chr', 'q1_start', 'q1_end',
        'q1_species', color_map
    )

    plot_genome(
        axes[2], df,
        'q2_chr', 'q2_start', 'q2_end',
        'q2_species', color_map
    )

    # -----------------------------
    # Legends
    # -----------------------------
    # ancestral colors
    anc_legend = [
        mpatches.Patch(color=color_map[c], label=c)
        for c in anc_chrs[:20]
    ]

    # breakpoint legend
    bp_legend = [
        mpatches.Patch(color='brown', label='Fusion'),
        mpatches.Patch(color='goldenrod', label='Fission'),
        mpatches.Patch(color='purple', label='Inversion'),
        mpatches.Patch(color='grey', label='Shared'),
        mpatches.Patch(color='black', label='Complex')
    ]

    axes[0].legend(
        handles=anc_legend,
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        title='Ancestral chromosomes',
        fontsize=8
    )

    axes[1].legend(
        handles=bp_legend,
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        title='Breakpoints',
        fontsize=9
    )

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"Saved figure → {args.output}")


if __name__ == '__main__':
    main()
