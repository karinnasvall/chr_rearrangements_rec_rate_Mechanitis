#!/usr/bin/env python3
"""
Detect breakpoints using ancestor-based comparisons

Takes PSL files for both genomes vs. the ancestor, identifies gaps
(breakpoints) that differ between the two genomes, and classifies them.
"""

import argparse
import pandas as pd
import sys
from collections import defaultdict


def parse_psl_to_blocks(psl_file, genome_name):
    """Parse PSL file and extract alignment blocks."""
    blocks = []
    
    with open(psl_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            fields = line.split('\t')
            if len(fields) < 21:
                continue
            
            q_chr = fields[9]
            q_start = int(fields[11])
            q_end = int(fields[12])
            t_chr = fields[13]
            t_start = int(fields[15])
            t_end = int(fields[16])
            strand = fields[8]
            
            block = {
                'query_chr': q_chr,
                'query_start': q_start,
                'query_end': q_end,
                'target_chr': t_chr,
                'target_start': t_start,
                'target_end': t_end,
                'strand': strand
            }
            blocks.append(block)
    
    return blocks


def detect_gaps_in_query(blocks):
    """
    Detect gaps (breakpoints) within query genome by finding:
    1. Chromosome boundaries
    2. Switch in target chromosome (fusion indicator)
    3. Switch in strand (inversion indicator)
    """
    if not blocks:
        return []
    
    # Sort by query chromosome and position
    sorted_blocks = sorted(blocks, key=lambda x: (x['query_chr'], x['query_start']))
    
    # Group by chromosome
    by_chr = defaultdict(list)
    for block in sorted_blocks:
        by_chr[block['query_chr']].append(block)
    
    gaps = []
    
    # Find gaps within each chromosome
    for q_chr in sorted(by_chr.keys()):
        chr_blocks = sorted(by_chr[q_chr], key=lambda x: x['query_start'])
        
        for i in range(len(chr_blocks) - 1):
            block1 = chr_blocks[i]
            block2 = chr_blocks[i + 1]
            
            gap_start = block1['query_end']
            gap_end = block2['query_start']
            
            if gap_end > gap_start:
                # Determine breakpoint type
                if block1['target_chr'] != block2['target_chr']:
                    bp_type = 'FUSION'
                elif block1['strand'] != block2['strand']:
                    bp_type = 'INVERSION'
                else:
                    bp_type = 'UNKNOWN'
                
                gap = {
                    'chromosome': q_chr,
                    'start': gap_start,
                    'end': gap_end,
                    'type': bp_type,
                    'gap_size': gap_end - gap_start,
                    'left_target_chr': block1['target_chr'],
                    'right_target_chr': block2['target_chr']
                }
                gaps.append(gap)
    
    return gaps


def detect_fission_in_target(blocks):
    """
    Detect fission breakpoints within target (ancestor) genome by finding:
    - Gaps in target where query chromosomes change
    This indicates the ancestor had a single chromosome that was split in the query
    """
    if not blocks:
        return []
    
    # Sort by target chromosome and position
    sorted_blocks = sorted(blocks, key=lambda x: (x['target_chr'], x['target_start']))
    
    # Group by target chromosome
    by_t_chr = defaultdict(list)
    for block in sorted_blocks:
        by_t_chr[block['target_chr']].append(block)
    
    fission_gaps = []
    
    # Find gaps within each target chromosome
    for t_chr in sorted(by_t_chr.keys()):
        chr_blocks = sorted(by_t_chr[t_chr], key=lambda x: x['target_start'])
        
        for i in range(len(chr_blocks) - 1):
            block1 = chr_blocks[i]
            block2 = chr_blocks[i + 1]
            
            gap_start = block1['target_end']
            gap_end = block2['target_start']
            
            if gap_end > gap_start:
                # Determine breakpoint type based on query chromosome continuity
                if block1['query_chr'] != block2['query_chr']:
                    bp_type = 'FISSION'
                elif block1['strand'] != block2['strand']:
                    bp_type = 'INVERSION'
                else:
                    bp_type = 'UNKNOWN'
                
                gap = {
                    'chromosome': t_chr,
                    'start': gap_start,
                    'end': gap_end,
                    'type': bp_type,
                    'gap_size': gap_end - gap_start,
                    'left_query_chr': block1['query_chr'],
                    'left_query_start': block1['query_start'],
                    'left_query_end': block1['query_end'],
                    'right_query_chr': block2['query_chr'],
                    'right_query_start': block2['query_start'],
                    'right_query_end': block2['query_end']
                }
                fission_gaps.append(gap)
    
    return fission_gaps


def main():
    parser = argparse.ArgumentParser(
        description='Detect breakpoints using ancestor-based comparisons'
    )
    parser.add_argument('--psl_g1', required=True, help='PSL file: genome1 vs ancestor')
    parser.add_argument('--psl_g2', required=True, help='PSL file: genome2 vs ancestor')
    parser.add_argument('--genome1', required=True, help='Genome 1 name')
    parser.add_argument('--genome2', required=True, help='Genome 2 name')
    parser.add_argument('--ancestor', required=True, help='Ancestor genome name')
    parser.add_argument('--min_block_size', type=int, default=50000,
                        help='Minimum block size in bp')
    
    args = parser.parse_args()
    
    # Parse PSL files
    print(f"Parsing {args.psl_g1}...", file=sys.stderr)
    blocks_g1 = parse_psl_to_blocks(args.psl_g1, args.genome1)
    print(f"  Parsed {len(blocks_g1)} blocks", file=sys.stderr)
    
    print(f"Parsing {args.psl_g2}...", file=sys.stderr)
    blocks_g2 = parse_psl_to_blocks(args.psl_g2, args.genome2)
    print(f"  Parsed {len(blocks_g2)} blocks", file=sys.stderr)
    
    # Detect breakpoints
    print(f"\nDetecting breakpoints in {args.genome1}...", file=sys.stderr)
    gaps_g1 = detect_gaps_in_query(blocks_g1)
    print(f"  Found {len(gaps_g1)} gaps in query genome", file=sys.stderr)
    
    fission_g1 = detect_fission_in_target(blocks_g1)
    print(f"  Found {len(fission_g1)} gaps in target genome (fission events)", file=sys.stderr)
    
    # Combine both types of breakpoints for genome1
    all_breakpoints_g1 = gaps_g1 + fission_g1
    
    print(f"\nDetecting breakpoints in {args.genome2}...", file=sys.stderr)
    gaps_g2 = detect_gaps_in_query(blocks_g2)
    print(f"  Found {len(gaps_g2)} gaps in query genome", file=sys.stderr)
    
    fission_g2 = detect_fission_in_target(blocks_g2)
    print(f"  Found {len(fission_g2)} gaps in target genome (fission events)", file=sys.stderr)
    
    # Combine both types of breakpoints for genome2
    all_breakpoints_g2 = gaps_g2 + fission_g2
    
    # Create dataframes and save
    if all_breakpoints_g1:
        df_g1 = pd.DataFrame(all_breakpoints_g1)
        df_g1.to_csv(f'{args.genome1}_breakpoints_classified.tsv', sep='\t', index=False)
        print(f"\nWrote {args.genome1}_breakpoints_classified.tsv ({len(df_g1)} total breakpoints)", file=sys.stderr)
    else:
        # Create empty file with headers
        pd.DataFrame(columns=['chromosome', 'start', 'end', 'type', 'gap_size', 
                             'left_target_chr', 'right_target_chr']).to_csv(
            f'{args.genome1}_breakpoints_classified.tsv', sep='\t', index=False)
        print(f"Wrote empty {args.genome1}_breakpoints_classified.tsv", file=sys.stderr)
    
    if all_breakpoints_g2:
        df_g2 = pd.DataFrame(all_breakpoints_g2)
        df_g2.to_csv(f'{args.genome2}_breakpoints_classified.tsv', sep='\t', index=False)
        print(f"Wrote {args.genome2}_breakpoints_classified.tsv ({len(df_g2)} total breakpoints)", file=sys.stderr)
    else:
        # Create empty file with headers
        pd.DataFrame(columns=['chromosome', 'start', 'end', 'type', 'gap_size',
                             'left_target_chr', 'right_target_chr']).to_csv(
            f'{args.genome2}_breakpoints_classified.tsv', sep='\t', index=False)
        print(f"Wrote empty {args.genome2}_breakpoints_classified.tsv", file=sys.stderr)


if __name__ == '__main__':
    main()
