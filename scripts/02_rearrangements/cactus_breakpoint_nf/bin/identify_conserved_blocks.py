#!/usr/bin/env python3
"""
Step 2: Identify conserved blocks from direct query1 vs query2 alignment

Step 2: Identify conserved blocks
- Reads direct query1 vs query2 alignment block file
- Creates two separate outputs:
  1. Sorted by query1_chr and query1_start, merged by query1_chr, query2_chr, orientation
  2. Sorted by query2_chr and query2_start, merged by query1_chr, query2_chr, orientation
- Both use the same conserved block IDs

Step 3: Detect unshared breakpoints
- Extracts start/end positions of each conserved block for both queries
- Creates separate dataframes for each query: query_chr, position, type (start/end), block_id
- Filters OUT chromosome boundaries shared between query1 and query2
- Keeps only query-specific breakpoints

Usage:
    python3 identify_conserved_blocks.py \
        --blocks query1_vs_query2_blocks.tsv \
        --outdir results
"""

import argparse
import pandas as pd
from pathlib import Path
import sys

def read_alignment_blocks(blocks_file):
    '''Read direct query1 vs query2 alignment blocks.'''
    print(f"\nReading query1 vs query2 alignment blocks from {blocks_file}")
    df = pd.read_csv(blocks_file, sep='\t')
    print(f"  Read {len(df)} blocks")
    return df


def sort_and_merge_blocks_by_query(df, query_num):
    '''
    Sort and merge blocks by specified query (1 or 2).
    Returns merged blocks sorted by that query's chromosome and start position.
    '''
    # Make a copy and rename columns for clarity
    df_sorted = df.copy()
    
    # For query2, need to swap the columns
    if query_num == 2:
        df_sorted = df_sorted.rename(columns={
            'query1_chr': 'query2_chr_temp',
            'query1_start': 'query2_start_temp',
            'query1_end': 'query2_end_temp',
            'query2_chr': 'query1_chr',
            'query2_start': 'query1_start',
            'query2_end': 'query1_end',
        })
        df_sorted = df_sorted.rename(columns={
            'query2_chr_temp': 'query2_chr',
            'query2_start_temp': 'query2_start',
            'query2_end_temp': 'query2_end',
        })
    
    # Sort by query1_chr with proper natural ordering:
    # autosomes (SUPER_1, SUPER_2, …) first, then sex chromosomes (SUPER_Z, SUPER_Z1, …)
    import re

    def _chr_sort_key(name):
        """Return (is_sex, numeric_part, full_name) for natural chromosome sorting."""
        suffix = re.sub(r'^SUPER_', '', str(name))
        is_sex = not suffix.isdigit()  # Z, Z1, Z2, W, etc. → True
        nums = re.findall(r'\d+', suffix)
        num = int(nums[0]) if nums else 0
        return (is_sex, num, suffix)

    df_sorted['_sort_key'] = df_sorted['query1_chr'].map(_chr_sort_key)
    df_sorted = df_sorted.sort_values(['_sort_key', 'query1_start']).reset_index(drop=True)
    df_sorted = df_sorted.drop('_sort_key', axis=1)
    
    # Merge consecutive rows with matching query1_chr, query2_chr, and orientation
    merged_records = []
    current_merge = None
    
    for idx, row in df_sorted.iterrows():
        q1_chr = row['query1_chr']
        q2_chr = row['query2_chr']
        orientation = row['strand']
        
        if current_merge is None:
            # Start new merge group
            current_merge = {
                'q1_chr': q1_chr,
                'q2_chr': q2_chr,
                'orientation': orientation,
                'blocks': [row]
            }
        else:
            # Check if this row continues the merge
            if (q1_chr == current_merge['q1_chr'] and 
                q2_chr == current_merge['q2_chr'] and 
                orientation == current_merge['orientation']):
                # Continue current merge
                current_merge['blocks'].append(row)
            else:
                # Save current merge and start new one
                merged_records.append(create_merged_block(current_merge))
                
                # Start new merge group
                current_merge = {
                    'q1_chr': q1_chr,
                    'q2_chr': q2_chr,
                    'orientation': orientation,
                    'blocks': [row]
                }
    
    # Don't forget the last merge group
    if current_merge:
        merged_records.append(create_merged_block(current_merge))
    
    df_merged = pd.DataFrame(merged_records)

    # If we sorted by query2, swap column names back so query1/query2 labels
    # match the original species (the data is just sorted by query2 coords)
    if query_num == 2:
        df_merged = df_merged.rename(columns={
            'query1_chr': 'query2_chr_temp',
            'query1_start': 'query2_start_temp',
            'query1_end': 'query2_end_temp',
            'query2_chr': 'query1_chr',
            'query2_start': 'query1_start',
            'query2_end': 'query1_end',
        })
        df_merged = df_merged.rename(columns={
            'query2_chr_temp': 'query2_chr',
            'query2_start_temp': 'query2_start',
            'query2_end_temp': 'query2_end',
        })
        # Reorder columns to standard order
        df_merged = df_merged[['query1_chr', 'query1_start', 'query1_end',
                                'query2_chr', 'query2_start', 'query2_end', 'strand']]

    return df_merged


def assign_block_ids(df):
    '''
    Assign block IDs based on unique chromosome pairs and strand orientation.
    The same chromosome pair with the same orientation gets the same block_id.
    This ensures consistency across query1 and query2 dataframes.
    '''
    # Create a canonical representation for each block
    # (chr1, chr2, strand) uniquely identifies a conserved block
    # For reversed blocks, we'll normalize them
    
    block_id_map = {}
    df_copy = df.copy()
    block_ids = []
    block_counter = 1
    
    for idx, row in df_copy.iterrows():
        chr1 = row['query1_chr']
        chr2 = row['query2_chr']
        strand = row['strand']
        
        # Create a key for this chromosome pair
        # Normalize so that comparisons across query1/query2 work
        # (e.g., SUPER_1 <-> SUPER_13 is the same as SUPER_13 <-> SUPER_1 but reversed)
        key = (chr1, chr2, strand)
        
        # Check if we've seen this exact combination before
        if key not in block_id_map:
            block_id_map[key] = f"conserved_block_{block_counter}"
            block_counter += 1
        
        block_ids.append(block_id_map[key])
    
    df_copy['block_id'] = block_ids
    return df_copy, block_id_map


def map_block_ids_for_query2(merged_by_q2, block_id_map):
    '''
    Map block IDs for query2 dataframe using the same mapping as query1.
    Handles the fact that query2 has chr1 and chr2 swapped.
    '''
    df_copy = merged_by_q2.copy()
    block_ids = []
    
    for idx, row in df_copy.iterrows():
        # In query2, chr1 and chr2 are swapped from the original alignment
        # So we need to look up (chr2, chr1, strand) in the map
        chr1 = row['query1_chr']
        chr2 = row['query2_chr']
        strand = row['strand']
        
        # Try to find this block in the map
        # It could be (chr1, chr2, strand) or (chr2, chr1, strand_reversed)
        key1 = (chr1, chr2, strand)
        key2 = (chr2, chr1, strand)  # Swapped version
        
        if key1 in block_id_map:
            block_ids.append(block_id_map[key1])
        elif key2 in block_id_map:
            block_ids.append(block_id_map[key2])
        else:
            # This is a new block not seen in query1
            # Assign a new ID
            new_id = f"conserved_block_{len(block_id_map) + len(block_ids) + 1}"
            block_ids.append(new_id)
    
    df_copy['block_id'] = block_ids
    return df_copy


def create_merged_block(merge_group):
    '''Create a merged block from consecutive blocks in a merge group.'''
    # merge_group is a dict with 'blocks' key containing list of rows
    block_list = merge_group['blocks']
    first_block = block_list[0]
    last_block = block_list[-1]
    
    record = {
        'query1_chr': first_block['query1_chr'],
        'query1_start': first_block['query1_start'],
        'query1_end': last_block['query1_end'],
        'query2_chr': first_block['query2_chr'],
        'query2_start': first_block['query2_start'],
        'query2_end': last_block['query2_end'],
        'strand': first_block['strand']
    }
    return record


def extract_breakpoints(merged_blocks_by_q1):
    '''
    Extract start/end positions of conserved blocks for both queries.
    Create separate dataframes for each query.
    Filter out only TRULY SHARED chromosome boundaries:
    - For same-strand blocks (+): start must match start, end must match end
    - For inverted blocks (-): start must match end, end must match start
    This preserves breakpoints that are specific to each genome even if both
    happen to be at chromosome boundaries (due to inversion or other rearrangement).
    '''
    print("\n" + "="*80)
    print("DETECTING UNSHARED BREAKPOINTS")
    print("="*80 + "\n")
    
    # Extract breakpoint positions for query1
    q1_records = []
    for idx, row in merged_blocks_by_q1.iterrows():
        block_id = row['block_id']
        q1_chr = row['query1_chr']
        q1_start = row['query1_start']
        q1_end = row['query1_end']
        
        q1_records.append({
            'query_chr': q1_chr,
            'position': q1_start,
            'type': 'start',
            'block_id': block_id
        })
        q1_records.append({
            'query_chr': q1_chr,
            'position': q1_end,
            'type': 'end',
            'block_id': block_id
        })
    
    df_q1 = pd.DataFrame(q1_records)
    
    # Extract breakpoint positions for query2
    q2_records = []
    for idx, row in merged_blocks_by_q1.iterrows():
        block_id = row['block_id']
        q2_chr = row['query2_chr']
        q2_start = row['query2_start']
        q2_end = row['query2_end']
        
        q2_records.append({
            'query_chr': q2_chr,
            'position': q2_start,
            'type': 'start',
            'block_id': block_id
        })
        q2_records.append({
            'query_chr': q2_chr,
            'position': q2_end,
            'type': 'end',
            'block_id': block_id
        })
    
    df_q2 = pd.DataFrame(q2_records)
    
    # Identify which breakpoints are at chromosome boundaries in each query
    # A position is a chromosome boundary if it's the first or last position in that chromosome's blocks
    q1_chr_boundaries = {}
    for chr_name in df_q1['query_chr'].unique():
        chr_data = df_q1[df_q1['query_chr'] == chr_name]
        q1_chr_boundaries[chr_name] = {
            'first': chr_data['position'].min(),  # First position (start of first block)
            'last': chr_data['position'].max()     # Last position (end of last block)
        }
    
    q2_chr_boundaries = {}
    for chr_name in df_q2['query_chr'].unique():
        chr_data = df_q2[df_q2['query_chr'] == chr_name]
        q2_chr_boundaries[chr_name] = {
            'first': chr_data['position'].min(),  # First position (start of first block)
            'last': chr_data['position'].max()     # Last position (end of last block)
        }
    
    # Find SHARED chromosome boundaries
    # A boundary is shared if it's a chromosome boundary in BOTH query1 AND query2 for the SAME block
    shared_block_ids = set()
    for idx, row in merged_blocks_by_q1.iterrows():
        block_id = row['block_id']
        q1_chr = row['query1_chr']
        q1_start = row['query1_start']
        q1_end = row['query1_end']
        q2_chr = row['query2_chr']
        q2_start = row['query2_start']
        q2_end = row['query2_end']
        strand = row['strand']
        

        # Check if block is at the very start/end of chromosome in BOTH queries
        q1_at_chr_start = (q1_start == q1_chr_boundaries[q1_chr]['first'])
        q1_at_chr_end = (q1_end == q1_chr_boundaries[q1_chr]['last'])
        
        q2_at_chr_start = (q2_start == q2_chr_boundaries[q2_chr]['first'])
        q2_at_chr_end = (q2_end == q2_chr_boundaries[q2_chr]['last'])
        
        # Account for chromosome inversion (strand orientation)
        # Strand format: '++' (same), '+-' (inverted)
        # Extract individual strand components
        strand_chars = list(strand) if isinstance(strand, str) else ['+', '+']
        q1_strand = strand_chars[0] if len(strand_chars) > 0 else '+'
        q2_strand = strand_chars[1] if len(strand_chars) > 1 else '+'
        
        is_shared = False
        
        if q1_strand == q2_strand:
            # Same orientation: starts match, ends match
            if (q1_at_chr_start and q2_at_chr_start) or (q1_at_chr_end and q2_at_chr_end):
                is_shared = True
        else:
            # Inverted orientation: start matches end, end matches start
            if (q1_at_chr_start and q2_at_chr_end) or (q1_at_chr_end and q2_at_chr_start):
                is_shared = True
        
        if is_shared:
            shared_block_ids.add(block_id)
    
    # Filter out breakpoints from shared boundary blocks (keep other breakpoints)
    df_q1_filtered = df_q1[~df_q1['block_id'].isin(shared_block_ids)]
    df_q2_filtered = df_q2[~df_q2['block_id'].isin(shared_block_ids)]
    
    print(f"Query1:")
    print(f"  All breakpoint positions: {len(df_q1)} (start/end of {len(df_q1)//2} blocks)")
    print(f"  Shared chromosome boundaries (conserved): {len(df_q1) - len(df_q1_filtered)}")
    print(f"  True breakpoints: {len(df_q1_filtered)}\n")
    
    print(f"Query2:")
    print(f"  All breakpoint positions: {len(df_q2)} (start/end of {len(df_q2)//2} blocks)")
    print(f"  Shared chromosome boundaries (conserved): {len(df_q2) - len(df_q2_filtered)}")
    print(f"  True breakpoints: {len(df_q2_filtered)}\n")
    
    return df_q1_filtered, df_q2_filtered


def main():
    parser = argparse.ArgumentParser(
        description='Identify conserved blocks from direct query1 vs query2 alignment'
    )
    parser.add_argument(
        '--blocks',
        required=True,
        help='Alignment blocks TSV file (query1 vs query2)'
    )
    parser.add_argument(
        '--outdir',
        default='.',
        help='Output directory [Default: current directory]'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    Path(args.outdir).mkdir(exist_ok=True, parents=True)
    
    print("="*80)
    print("IDENTIFY CONSERVED BLOCKS AND UNSHARED BREAKPOINTS")
    print("="*80)
    
    # Read input blocks
    print(f"\nReading alignment blocks from {args.blocks}")
    blocks_df = read_alignment_blocks(args.blocks)
    
    # Step 2a: Sort and merge by query1
    print("\n" + "="*80)
    print("STEP 2a: MERGE BLOCKS SORTED BY QUERY1")
    print("="*80)
    merged_by_q1 = sort_and_merge_blocks_by_query(blocks_df, query_num=1)
    merged_by_q1, block_id_map = assign_block_ids(merged_by_q1)
    print(f"Merged: {len(blocks_df)} blocks → {len(merged_by_q1)} conserved blocks")
    print(f"Unique conserved blocks (chromosome pairs): {len(block_id_map)}")
    
    # Step 2b: Sort and merge by query2
    print("\n" + "="*80)
    print("STEP 2b: MERGE BLOCKS SORTED BY QUERY2")
    print("="*80)
    merged_by_q2 = sort_and_merge_blocks_by_query(blocks_df, query_num=2)
    merged_by_q2 = map_block_ids_for_query2(merged_by_q2, block_id_map)
    print(f"Merged: {len(blocks_df)} blocks → {len(merged_by_q2)} conserved blocks")
    print(f"NOTE: Block IDs are now consistent with query1 based on chromosome pairs")
    
    # Save merged blocks
    output_q1 = Path(args.outdir) / 'conserved_blocks_by_query1.tsv'
    output_q2 = Path(args.outdir) / 'conserved_blocks_by_query2.tsv'
    
    merged_by_q1.to_csv(output_q1, sep='\t', index=False)
    merged_by_q2.to_csv(output_q2, sep='\t', index=False)
    
    print(f"\nSaved conserved blocks (sorted by query1) to {output_q1}")
    print(f"Saved conserved blocks (sorted by query2) to {output_q2}")
    
    # Step 3: Extract breakpoints
    print("\n" + "="*80)
    print("STEP 3: DETECT UNSHARED BREAKPOINTS")
    print("="*80)
    df_q1_bp, df_q2_bp = extract_breakpoints(merged_by_q1)
    
    # Save breakpoint dataframes
    output_q1_bp = Path(args.outdir) / 'query1_breakpoints.tsv'
    output_q2_bp = Path(args.outdir) / 'query2_breakpoints.tsv'
    
    df_q1_bp.to_csv(output_q1_bp, sep='\t', index=False)
    df_q2_bp.to_csv(output_q2_bp, sep='\t', index=False)
    
    print(f"\nSaved query1 breakpoints to {output_q1_bp}")
    print(f"Saved query2 breakpoints to {output_q2_bp}")
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)


if __name__ == '__main__':
    main()
