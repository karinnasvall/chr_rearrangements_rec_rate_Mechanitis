import pandas as pd
import argparse
import logging
import os

# Setup logging
logging.basicConfig(filename="repeatmasker_to_bed.log", level=logging.INFO, 
                    format="%(asctime)s - %(levelname)s - %(message)s")


# Function to convert RepeatMasker .out to BED
def out_to_bed(out_file, output_bed):
    logging.info(f"Reading RepeatMasker .out file: {out_file}")
    bed_records = []

    with open(out_file, 'r') as file:
        for line in file:
            if line.strip() == "" or line.startswith("SW") or line.startswith("score"):
                continue

            parts = line.strip().split()
            if len(parts) < 15:
                continue

            chrom = parts[4]
            start = int(parts[5]) - 1  # BED is 0-based
            end = int(parts[6])
            class_family = parts[10]  # Extracting the class/family
            percent_div = parts[1]   # Extracting the divergence

            bed_records.append([chrom, start, end, class_family, percent_div])

    # Save to BED
    logging.info(f"Writing BED file: {output_bed}")
    bed_df = pd.DataFrame(bed_records)
    bed_df.to_csv(output_bed, sep="\t", header=False, index=False)

# Command-line interface (CLI)
def main():
    parser = argparse.ArgumentParser(description="Convert RepeatMasker .out to BED format.")
    parser.add_argument("--prefix", required=True, help="Prefix for input and output files")
    parser.add_argument("--input_file", required=True, help="Input file, ie the output (.out) from RepeatMasker")
    parser.add_argument("--output_dir", required=True, help="Directory to save output BED file")

    args = parser.parse_args()

    if args.input_file:
        out_file = args.input_file
        if not args.prefix:
            args.prefix = os.path.splitext(os.path.basename(out_file))[0]
    else:
        if not args.prefix:
            print("Error: Please specify either --prefix or --input_file.")
            return
        out_file = os.path.join(f"{args.prefix}.out")

    output_bed = os.path.join(args.output_dir, f"{args.prefix}_repeats.bed")

    if not os.path.exists(out_file):
        logging.error(f"Input file does not exist: {out_file}")
        print(f"Error: Input file does not exist: {out_file}")
        return

    os.makedirs(args.output_dir, exist_ok=True)

    logging.info("Starting RepeatMasker .out to BED conversion")
    out_to_bed(out_file, output_bed)
    logging.info("Conversion completed.")

if __name__ == "__main__":
    main()
