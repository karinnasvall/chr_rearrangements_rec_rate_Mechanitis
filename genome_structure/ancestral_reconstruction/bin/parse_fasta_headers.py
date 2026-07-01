#!/usr/bin/env python3
"""
Parse FASTA headers and convert to TSV format for AGORA.

Input header format:
>109721at7088_69820_0:000838|SUPER_9:3985955-3986581|+ <unknown description>

Output TSV format (no header):
chromosome  start  end  orientation  gene_id
"""

import sys


def parse_header(header):
    """
    Parse FASTA header and extract gene information.

    Returns:
        tuple: (chromosome, start, end, orientation, gene_id)
        or None if parsing fails
    """
    # Remove '>' and split by whitespace to get just the header part
    header = header.lstrip('>')
    header_part = header.split()[0] if ' ' in header else header

    # Split by pipe '|' to extract chromosome, coordinates, and orientation
    parts = header_part.split('|')
    if len(parts) < 3:
        return None

    # Parse chromosome and coordinates from second part
    chr_coords = parts[1]
    if ':' not in chr_coords:
        return None
    chromosome, coords = chr_coords.split(':', 1)

    if '-' not in coords:
        return None
    start, end = coords.split('-', 1)

    # Parse orientation from third part
    orientation = parts[2].strip()
    if orientation == '+':
        orientation = '1'
    elif orientation == '-':
        orientation = '-1'
    else:
        return None

    # Gene ID: everything except the orientation part, with special chars replaced
    gene_id = '|'.join(parts[:-1])
    gene_id = gene_id.replace('|', '_').replace(':', '_')

    try:
        start = int(start)
        end = int(end)
    except ValueError:
        return None

    return chromosome, start, end, orientation, gene_id


def main():
    """Process FASTA file and output TSV."""
    if len(sys.argv) < 2:
        print("Usage: parse_fasta_headers.py <fasta_file> [output_file]", file=sys.stderr)
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None

    output = open(output_file, 'w') if output_file else sys.stdout

    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    result = parse_header(line.strip())
                    if result:
                        chromosome, start, end, orientation, gene_id = result
                        output.write(f"{chromosome}\t{start}\t{end}\t{orientation}\t{gene_id}\n")
    finally:
        if output_file:
            output.close()


if __name__ == '__main__':
    main()
