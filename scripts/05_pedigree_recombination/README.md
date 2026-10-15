# 05_pedigree_recombination

Linkage-map construction and crossover-based recombination rate estimation from the F2 pedigrees, using Lep-MAP3 (Rastas, 2017, *Bioinformatics* 33:3726–3732). Scripts are separated by species:

| Species | Pedigree  | Folder |
|---|---|---|
| *M. lysimnia* | two families, 48 and 49 F2 offspring, analysed jointly | `scripts_lysimnia/` |
| *M. polymnia* | one family, 28 F2 offspring | `scripts_polymnia/` |

Only males recombine in Lepidoptera: female recombination is fixed to zero and the W chromosome is excluded.

## Workflow

Numbers in file names give the run order. Step 04 Join Singles was not performed as it did not yield more markers.

| Step | Script(s) | Lep-MAP3 module / task |
|---|---|---|
| 01 | `01_*` | `Pileup2Likelihoods` (Genotype likelihoods from deduplicated bam-files) |
| 02 | `02_parentcall2.bsub`,  `02_Filtering2.bsub` | `ParentCall2` / `Filtering2` (Infer parental haplotypes) |
| 03 | `03_Separate_chromosomes2.bsub` | `SeparateChromosomes2` (Linkage groups) |
| - | `scripts/create_lod_summary.sh` | Summarises LOD-limit tests (1–15) |
| - | `scripts/plot_sep_chr.R` | Plots linkage-group separation |
| 05 | `05_order_*` | `OrderMarkers2` (Infer marker order, recombination fractions) |
| 06 | `06_fitstepfunction.bsub` | Step-function fit of the ordered maps to physical positions |
| – | R scripts | Plots, crossover counts, recombination rates |
| - | `scripts/map_physical_pos_to_ordered_maps.sh` | Maps markers to phyical position |
| - | `scripts/lepmap_plot_order.R`, `lepmap_plot_map_fitted.R` | Plots of marker order and fitted maps |
| - | `scripts/plot_marey_rec_rate_cobs_smoothing.R` | Marey maps; local rate per 100 kb with cobs |
| - | `scripts/summarize_rec_events.R` | Crossovers per chromosome per individual |
| - | `scripts/compare_pop_ped_rec_rate.R` | Spearman correlation of pedigree and pyrho rates |


Additional notes:
- Genotype likelihoods: samtools `mpileup` (MAPQ ≥ 20, baseQ ≥ 20) and `Pileup2Likelihoods` on bam-fiels mapped to repeat-masked references.
- One *M. lysimnia* F2 (CAM073484) was excluded after the `IBD` check.
- *M. polymnia* chromosome 2 shows a 25 cM gap; markers informative in both parents give 16 linkage groups, with chromosome 2 split into 2A and 2B. The recombination analysis uses the split state.

## Inputs and outputs

- Input: BAM files (see `../04_data_processing`), repeat-masked references, pedigree file.
- Output: linkage maps in `data/results/linkage_maps/`; per chromosome and per-window (100kb) rates (cM/Mb)).

## Software

Lep-MAP3, samtools, R (cobs, ggplot2, tidyverse), Nextflow (pipelines for BAMs).

