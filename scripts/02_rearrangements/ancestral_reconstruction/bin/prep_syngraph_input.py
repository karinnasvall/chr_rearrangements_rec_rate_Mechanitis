#!/usr/bin/env python3
"""
Prepare Syngraph input from OrthoFinder results.

Uses the per-species pairwise Orthologues TSV files
(results/Orthologues/{species}.tsv) to extract single-copy orthologues.
Rows with exactly 4 tab-separated columns and no commas in the gene
fields represent single-copy 1:1 orthologues.

Gene format in column 3 (focal species):
  busco_gene_name|CHROM_start-end|strand

Output per species: busco_id <TAB> chromosome <TAB> start <TAB> end
(start > end when gene is on minus strand)
"""

import os
import sys
import re


EXCLUDE_PATTERN = re.compile(r"SCAFF|UNLOC|unloc|scaffold", re.IGNORECASE)


def parse_gene_field(gene_field):
    """
    Parse a gene field like:
      607839at7088_278856_0_000e95|SUPER_17_2227585-2228774|+

    Returns (busco_id, chrom, start, end) or None.
    """
    parts = gene_field.split('|')
    if len(parts) < 3:
        return None

    gene_name = parts[0]
    location = parts[1]
    strand = parts[2].strip()

    # BUSCO ID = leading digits + "at" + digits
    busco_match = re.match(r'(\d+at\d+)', gene_name)
    if not busco_match:
        return None
    busco_id = busco_match.group(1)

    # Location: SUPER_17_2227585-2228774
    # Split on '_', last element has 'start-end'
    loc_parts = location.split('_')
    if len(loc_parts) < 2:
        return None

    coords_str = loc_parts[-1]
    if '-' not in coords_str:
        return None

    start_s, end_s = coords_str.split('-', 1)
    chrom = '_'.join(loc_parts[:-1])

    try:
        start = int(start_s)
        end = int(end_s)
    except ValueError:
        return None

    # Swap start/end for minus strand so start > end signals direction
    if strand == '-':
        start, end = end, start

    return busco_id, chrom, start, end


def process_species_tsv(tsv_path, output_path):
    """
    Read an OrthoFinder Orthologues/{species}.tsv and write syngraph input.

    Filters for:
      - Exactly 4 tab columns (Orthogroup, Species, FocalGene, QueryGene)
      - No commas in columns 3 or 4 (single-copy only)
      - Chromosome not matching scaffold/unloc patterns
    Deduplicates by BUSCO ID (keeps first occurrence).
    """
    seen = set()
    kept = 0

    with open(tsv_path) as fin, open(output_path, 'w') as fout:
        for line in fin:
            line = line.rstrip('\n')
            cols = line.split('\t')

            # Skip header and non-4-column rows
            if len(cols) != 4:
                continue
            if cols[0] == 'Orthogroup':
                continue

            focal_gene = cols[2]
            query_gene = cols[3]

            # Skip multi-copy (comma-separated genes)
            if ',' in focal_gene or ',' in query_gene:
                continue

            parsed = parse_gene_field(focal_gene)
            if not parsed:
                continue

            busco_id, chrom, start, end = parsed

            if EXCLUDE_PATTERN.search(chrom):
                continue

            # Deduplicate by BUSCO ID
            if busco_id in seen:
                continue
            seen.add(busco_id)

            fout.write(f"{busco_id}\t{chrom}\t{start}\t{end}\n")
            kept += 1

    return kept


def main():
    if len(sys.argv) < 4:
        print(
            "Usage: prep_syngraph_input.py <orthofinder_results_dir> "
            "<fasta_dir> <output_dir>",
            file=sys.stderr,
        )
        sys.exit(1)

    of_results = sys.argv[1]
    fasta_dir = sys.argv[2]   # kept for CLI compatibility (unused)
    output_dir = sys.argv[3]
    os.makedirs(output_dir, exist_ok=True)

    # Find per-species Orthologues TSVs
    ortho_dir = os.path.join(of_results, "Orthologues")
    if not os.path.isdir(ortho_dir):
        print(
            f"ERROR: Orthologues directory not found: {ortho_dir}",
            file=sys.stderr,
        )
        sys.exit(1)

    n_species = 0
    for fname in sorted(os.listdir(ortho_dir)):
        if not fname.endswith('.tsv'):
            continue
        species = fname.replace('.tsv', '')
        tsv_path = os.path.join(ortho_dir, fname)
        output_path = os.path.join(output_dir, f"{species}.tsv")

        kept = process_species_tsv(tsv_path, output_path)
        print(f"[syngraph prep] {species}: {kept} genes", file=sys.stderr)
        n_species += 1

    print(
        f"[syngraph prep] Created input for {n_species} species in {output_dir}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
