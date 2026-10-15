# Genome annotation workflow

This directory contains a small collection of batch scripts and helper programs used to prepare, annotate, and summarise genome assemblies before downstream synteny and rearrangement analyses.

- `helixer.bsub`  
  LSF batch job for Helixer annotation. Runs annotation from genome FASTA to GFF3, extracts CDS sequences with gffread, and runs BUSCO on the annotated protein set.

- `merge_hap_index.bsub`  
  Merges the primary and alternate haplotype assembly chunks into a single FASTA for downstream use, then indexes the merged assembly with samtools.

- `pantera_ml.bsub`  
  Builds a repeat library with PANTEra, classifies repeats, runs RepeatMasker on the assembly, converts the RepeatMasker output to BED, creates genomic windows, estimates repeat occupancy, and summarises repeat content per chromosome.

- `pggb.bsub`  
  Runs pggb on a merged genome assembly to generate a pangenome graph for a set of species/assemblies.

- `get_windowbased_content.sh`  
  Generic script to compute feature density across fixed genomic windows for any GFF feature type (for example CDS, repeat_region, gene). It creates windows, counts overlaps per window, calculates base-pair occupancy, and generates per-chromosome summaries.

- `get_windowbased_content_run.bsub`  
  LSF wrapper that runs the feature-density pipeline for a list of genomes, using CDS as the default feature and a 100 kb window size.

- `repeatmasker_to_bed.py`  
  Converts RepeatMasker `.out` output to a BED-like table with chromosome, start, end, repeat class/family, and divergence information.

- `repeat_occupancy_calculator_plot.py`  
  Estimates repeat occupancy across fixed windows and generates summary plots.

- `summarise_repeat_per_chr.sh`  
  Collapses repeat intervals by chromosome and repeat class and reports total covered bases and counts using `awk` and `bedtools groupby`.

- `bedtools_nuc.bsub`  
  Calculates window-based GC-content using bedtools nuc.

## Workflow

Repeat annotation
1. Merge haplotype files into a single assembly using `merge_hap_index.bsub`.
2. Build a pangenome with `pggb.bsub` to use as input for pantera.
3. Build an annotation/repeat library with `pantera_ml.bsub`.
4. Run RepeatMasker and generate a `.out` file.
5. Convert the `.out` file to BED using `repeatmasker_to_bed.py`.
6. Summarise repeat content per chromosome with `summarise_repeat_per_chr.sh`.

Gene annotation
7. Use `helixer.bsub` to annotate gene models in fasta (hap1 incl W)
8. For broader feature-density summaries, run `get_windowbased_content.sh` (wrapped by `get_windowbased_content_run.bsub`) to compute CDS or other feature counts/occupancy in fixed-size windows.
9. Use the repeat occupancy summaries, feature-density tables, and BED outputs as inputs for downstream rearrangement/assembly analyses.

## Expected inputs

These scripts assume upstream data are available in project directories such as:

- genome FASTAs in the project data folders,
- index files (`.fai`),
- per-species assembly lists

The exact paths are hard-coded in the job scripts, so check the variables near the top of each `.bsub`.

## Notes

- These are primarily cluster job scripts for an LSF environment (`bsub`).
- Several scripts reference project-specific directories and environment modules, so they are not meant to be run without adapting paths and environment names.
