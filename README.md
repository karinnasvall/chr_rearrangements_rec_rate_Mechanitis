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

| Folder | Contents | Main tools |
|---|---|---|
| `scripts/01_genome_annotation/` | Gene and repeat annotation and 100 kb window features | Helixer, pantera, RepeatMasker |
| `scripts/02_rearrangements/` | Characterisation of chromosomal rearrangements, ancestral reconstruction and breakpoint detection | syngraph, AGORA, minimap2, ODP |
| `scripts/03_chromatin/` | Hi-C processing and A/B compartments | nf-core/hic, HiC-Pro, Cooler, cooltools |
| `scripts/04_data_processing/` | Read processing. Trimming, mapping and variant calling of population and pedigree data | fastp, bwa-mem2, samtools, bcftools |
| `scripts/05_pedigree_recombination/` | Lep-MAP3 linkage maps and crossover-based rates from the F2 pedigrees | Lep-MAP3 |
| `scripts/06_population_recombination/` | Population recombination maps | pyrho |
| `scripts/07_recombination_models/` | Recombination rate as a function of genome structure and genomic features | R: lme4 |
| `scripts/08_pop_gen/` | Population structure and diversity/divergence statistics | ADMIXTURE, genomics_general |
| `scripts/09_combined_analysis/` | Main analyses and figures | R |


Scripts were written for an LSF cluster (`.bsub`); adjust resource requests for other schedulers.


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


## Citation

Please cite the article above and the archived release of this repository ([Zenodo DOI]). Machine-readable metadata are in `CITATION.cff`.

## Licence

Code: [MIT / GPL-3.0: choose]. Data in `data/results/`: [CC BY 4.0: choose].

## Contact

Karin Näsvall (karin.nasvall@sanger.ac.uk)
