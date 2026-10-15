# Rearrangements Analysis (02_rearrangements)

This directory contains workflows and helper scripts for ancestral reconstruction, breakpoint detection, whole-genome alignment refinement, and synteny plotting.

## Directory Map

- [analysis_steps.txt](analysis_steps.txt): high-level step log for the rearrangements workflow.
- [cactus.bsub](cactus.bsub): batch script to run cactus alignment.
- [ancestral_reconstruction](ancestral_reconstruction): BUSCO/Orthofinder/Syngraph/AGORA workflow and helper scripts.
- [cactus_breakpoint_nf](cactus_breakpoint_nf): Nextflow pipeline for extracting conserved blocks and ancestral breakpoints from cactus output.
- [minimap](minimap): pairwise minimap-based synteny inputs and plotting script.
- [odp](odp): Multiple synteny alignment with ODP config/list files and batch scripts.

## Key Files

### Ancestral reconstruction workflow

- [main.nf](ancestral_reconstruction/main.nf)
- [nextflow.config](ancestral_reconstruction/nextflow.config)
- [run_nf.bsub](ancestral_reconstruction/run_nf.bsub)
- [syngraph_agora_results.Rmd](ancestral_reconstruction/syngraph_agora_results.Rmd)

### Ancestral reconstruction helper scripts (`bin`)

- [run_breakpoint_analysis.sh](ancestral_reconstruction/bin/run_breakpoint_analysis.sh)
- [detect_breakpoints_agora.py](ancestral_reconstruction/bin/detect_breakpoints_agora.py)
- [refine_breakpoints.py](ancestral_reconstruction/bin/refine_breakpoints.py)
- [plot_breakpoints_agora.py](ancestral_reconstruction/bin/plot_breakpoints_agora.py)
- [plot_tree.py](ancestral_reconstruction/bin/plot_tree.py)
- [run_syntenyplotter_busco_agora.R](ancestral_reconstruction/bin/run_syntenyplotter_busco_agora.R)
- [correct_breakpoints.R](ancestral_reconstruction/bin/correct_breakpoints.R)

### Cactus breakpoint Nextflow

- [main.nf](cactus_breakpoint_nf/main.nf)
- [nextflow.config](cactus_breakpoint_nf/nextflow.config)
- [run_nf.bsub](cactus_breakpoint_nf/run_nf.bsub)
- [identify_conserved_blocks.py](cactus_breakpoint_nf/bin/identify_conserved_blocks.py)
- [detect_ancestral_breakpoints.py](cactus_breakpoint_nf/bin/detect_ancestral_breakpoints.py)
- [plot_results.py](cactus_breakpoint_nf/bin/plot_results.py)

## Path conventions

- In this module, helper scripts are expected under local `bin/` directories:
  - [ancestral_reconstruction/bin](ancestral_reconstruction/bin)
  - [cactus_breakpoint_nf/bin](cactus_breakpoint_nf/bin)


