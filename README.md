# chr_rearrangements_rec_rate_Mechanitis

This repository contain scripts used in the analyses in the manuscript "Chromosomal rearrangements remodel the recombination landscape in holocentric butterflies" 

Authors: Karin Näsvall, María José Sánchez-Carvajal, Alex P. Arias-Cruz, Franz Chandi, Shakirah Namuli, Alice Groudko, Nicol Rueda-M, Patricio Salazar-Carrión, Eva S.M. van der Heijden, Karen Brooks, Danil Zilov, Jonathan M. D. Wood´, Camilo Salazar, Marianne Elias, Keith Willmott, Caroline Bacquet, Joana I. Meier


Project structure for rearrangement and recombination analyses.

## Directory overview

- `data/`: metadata, linkage maps and recombination data
- `scripts/`: analysis scripts organised by analysis workflow

## Scripts directory structure

The `scripts/` directory currently contains the early analysis stages that are already established in this project tree. The remaining analysis modules are kept in the sibling directory `chr_rearrangements_rec_rate_Mechanitis_extra/` until they are migrated or recreated here.

```text
scripts/
├── 01_genome_annotation/
│   ├── README.md
│   ├── helixer.bsub
│   ├── merge_hap_index.bsub
│   ├── pantera_ml.bsub
│   ├── pggb.bsub
│   ├── get_windowbased_content.sh
│   ├── get_windowbased_content_run.bsub
│   ├── repeatmasker_to_bed.py
│   ├── repeat_occupancy_calculator_plot.py
│   └── summarise_repeat_per_chr.sh
│
├── 02_rearrangements/
│   ├── README.md
│   ├── analysis_steps.txt
│   ├── ancestral_reconstruction/
│   ├── cactus_breakpoint_nf/
│   ├── minimap/
│   ├── odp/
│   ├── cactus.bsub
│   └── ...
│
└── [future stages]

```

## Notes

- `01_genome_annotation/` contains the genome assembly, repeat annotation and feature-density workflows.
- `02_rearrangements/` contains the breakpoint and ancestral reconstruction workflows.



