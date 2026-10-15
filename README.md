# Chromosomal rearrangements remodel the recombination landscape in holocentric butterflies

Code, linkage maps and recombination maps accompanying the manuscript:

> Näsvall K, Sánchez-Carvajal MJ, Arias-Cruz AP, Chandi F, Namuli S, Groudko A, Rueda-M N, Salazar-Carrión P, van der Heijden ESM, Brooks K, Zilov D, Wood JMD, Salazar C, Elias M, Willmott K, Bacquet C, Meier JI. *Chromosomal rearrangements remodel the recombination landscape in holocentric butterflies.* In review.

Species: *Mechanitis lysimnia* and *M. polymnia* (Ithomiini, Lepidoptera).
Archive: [Zenodo DOI: **will be added at release**] · Raw reads for population and pedigree data in `data/metadata/` (ENA project [PRJEB63068](https://www.ebi.ac.uk/ena/browser/view/PRJEB63068)). Assembly accessions in `data/metadata/`.

---

## Repository layout

```
.
├── README.md
├── LICENSE
├── CITATION.cff
├── data/
│   ├── metadata/           # Metadata from assemblies, population and pedigree data
│   └── results/            # Final linkage maps and recombination maps, breakpoints
├── scripts/                # All analysis code, numbered in run order (see below)

```

---

## Scripts directory structure

Folders are numbered in analysis order.

| Folder | Contents | README |
|---|---|---|
| `scripts/01_genome_annotation/` | Gene and repeat annotation and window-based feature densities |
| `scripts/02_rearrangements/` | Breakpoint detection and ancestral reconstruction |
| `scripts/03_chromatin/` | Hi-C processing and A/B compartments |
| `scripts/04_data_processing/` | Trimming, mapping and variant calling of population and pedigree data |
| `scripts/05_pedigree_recombination/` | Lep-MAP3 linkage maps and crossover-based rates from the F2 pedigrees |
| `scripts/06_population_recombination/` | Population recombination maps (pyrho) |
| `scripts/07_recombination_models/` | Recombination rate as a function of genome structure and genomic features |
| `scripts/08_pop_gen/` | Population structure (ADMIXTURE) and diversity/divergence statistics |
| `scripts/09_combined_analysis/` | Main analyses and figures |

```text
scripts
├── 01_genome_annotation
│   ├── get_windowbased_content_run.bsub
│   ├── get_windowbased_content.sh
│   ├── helixer.bsub
│   ├── merge_hap_index.bsub
│   ├── pantera_ml.bsub
│   ├── pggb.bsub
│   ├── README.md
│   ├── repeatmasker_to_bed.py
│   ├── repeat_occupancy_calculator_plot.py
│   └── summarise_repeat_per_chr.sh
├── 02_rearrangements
│   ├── analysis_steps.txt
│   ├── ancestral_reconstruction
│   ├── cactus_breakpoint_nf
│   ├── cactus.bsub
│   ├── minimap
│   ├── odp
│   └── README.md
├── 03_chromatin
├── 04_data_processing
│   ├── 01_trimming
│   ├── 02_mapping
│   ├── 03_variant_calling
│   └── 04_filter_pop_data
├── 05_pedigree_recombination
│   ├── 01_concatenate_pileup2lik.bsub
│   ├── 01_pileup2lik_from_nextflow_bams.bsub
│   ├── ibd.bsub
│   ├── README.md
│   ├── scripts_lysimnia
│   └── scripts_polymnia
├── 06_population_recombination
│   ├── 01a_bcftools_filter.bsub
│   ├── 01b_smc_pipeline.bsub
│   ├── 02_pyrho_make_table.bsub
│   ├── 03_pyrho_hyperparam.bsub
│   └── 04_pyrho_optimize.bsub
├── 07_recombination_models
│   └── model_rec_rate.Rmd
├── 08_pop_gen
│   ├── admixture
│   └── popgen
└── 09_combined_analysis
    ├── bin
    ├── data_preparation_summary.Rmd
    ├── Figure_2_genome_tracks.R
    └── figure_3_combined_plot_261007.Rmd
```

### Mapping manuscript sections to code

| Manuscript section (Methods) | Folder | Main tools | Key outputs |
|---|---|---|---|
| Gene annotation, repeats, 100 kb window features | `01_genome_annotation/` | Helixer, pantera, RepeatMasker | gene and repeat occupancy per window |
| Characterisation of chromosomal rearrangements | `02_rearrangements/` | syngraph, AGORA, minimap2, ODP | Breakpoint table (type, polarity, coordinates), Figure 1 |
| Chromatin compartments | `03_chromatin/` | nf-core/hic, HiC-Pro, Cooler, cooltools | A/B compartment calls (100 kb) |
| Read processing, alignment and variant calling | `04_data_processing/` | fastp, bwa-mem2, samtools, bcftools | Filtered bam-files, VCFs, QC and coverage summaries |
| Pedigree recombination | `05_pedigree_recombination/` | Lep-MAP3 | Linkage maps, crossover counts |
| Population recombination | `06_population_recombination/` | pyrho | Population recombination maps |
| Recombination across the genome | `07_recombination_models/` | R: lme4  | Model objects, diagnostics |
| Differentiation and gene flow | `08_pop_gen/` | ADMIXTURE, genomics_general | π, *D*<sub>XY</sub>, *F*<sub>ST</sub> per 100 kb |
| Figures and supplementary tables | `09_combined_analysis/` | R | Figures 2–3 and supplementary material |



Scripts were written for an LSF cluster (`.bsub`); adjust resource requests for other schedulers.


## Software versions

Principal versions reported in the manuscript: BUSCO 6.0.0, OrthoFinder 3.1.0, BEDTools 2.31.1, Helixer 0.3.5, OMArk 0.3.1, nf-core/hic 2.1.0, fastp 0.23.4, bwa-mem2 2.2.1, bcftools 1.20.



## Citation

Please cite the article above and the archived release of this repository ([Zenodo DOI]). Machine-readable metadata are in `CITATION.cff`.

## Licence

Code: [MIT / GPL-3.0: choose]. Data in `data/results/`: [CC BY 4.0: choose].

## Contact

Karin Näsvall (karin.nasvall@sanger.ac.uk)
