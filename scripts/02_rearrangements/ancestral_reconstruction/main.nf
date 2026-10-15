#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * =====================================================
 *  BREAKPOINT DETECTION PIPELINE
 * =====================================================
 *  Steps:
 *   1. Prepare genomes (move W, filter scaffolds, index)
 *   2. BUSCO on both haplotypes
 *   3. OrthoFinder on hap1 single-copy BUSCOs
 *   4. Syngraph ancestral genome reconstruction
 *   5. AGORA ancestral genome ordering
 * =====================================================
 */

// =====================================================
// Parameter validation
// =====================================================

if (!params.genome_list) {
    error "ERROR: --genome_list is required (text file with one genome path per line)"
}
if (!params.ref_genome) {
    error "ERROR: --ref_genome is required (reference tolid for Syngraph, e.g. ilForEqui1)"
}

log.info """
=====================================================
 BREAKPOINT DETECTION PIPELINE
=====================================================
 genome_list       : ${params.genome_list}
 outdir            : ${params.outdir}
 min_scaffold_size : ${params.min_scaffold_size}
 move_w            : ${params.move_w}
 busco_lineage     : ${params.busco_lineage}
 skip_busco        : ${params.skip_busco}
 skip_orthofinder  : ${params.skip_orthofinder}
 ref_genome        : ${params.ref_genome}
 user_tree         : ${params.user_tree ?: 'auto (OrthoFinder inferred)'}
 syngraph_m_range  : ${params.syngraph_m_min}-${params.syngraph_m_max}
=====================================================
"""

// =====================================================
// Processes
// =====================================================

/*
 * Prepare genomes: move W chromosomes, filter scaffolds, and index
 * - If move_w: extracts W scaffolds from hap1, renames to hap2, appends to hap2
 * - Filters out scaffold/SCAFF/unloc patterns and scaffolds < min_scaffold_size
 * - Indexes with samtools faidx
 * Outputs filtered+indexed FASTA per haplotype
 */
process PREPARE_GENOMES {
    tag "${tolid}"
    label 'small'
    publishDir "${params.outdir}/01_genomes", mode: 'copy'

    input:
    tuple val(tolid), path(hap1), path(hap2)

    output:
    tuple val(tolid),
          path("${tolid}_hap1.fa"), path("${tolid}_hap1.fa.fai"),
          path("${tolid}_hap2.fa"), path("${tolid}_hap2.fa.fai")

    script:
    def move_w = params.move_w ? "true" : "false"
    """
    # ---- Move W chromosomes if requested ----
    if [ "${move_w}" = "true" ]; then
        seqkit grep -r -p "_W" ${hap1} | seqkit replace -p "_1_" -r "_2_" > W_scaffolds.fa || true
        seqkit grep -r -p "_W" -v ${hap1} > hap1_noW.fa
        if [ -s W_scaffolds.fa ]; then
            cat ${hap2} W_scaffolds.fa > hap2_plusW.fa
        else
            cp ${hap2} hap2_plusW.fa
        fi
    else
        cp ${hap1} hap1_noW.fa
        cp ${hap2} hap2_plusW.fa
    fi

    # ---- Filter scaffolds by name and minimum size ----
    seqkit grep -v -r -p "scaffold|SCAFFOLD|SCAFF|unloc|UNLOC" hap1_noW.fa | \
        seqkit seq -m ${params.min_scaffold_size} -o ${tolid}_hap1.fa
    seqkit grep -v -r -p "scaffold|SCAFFOLD|SCAFF|unloc|UNLOC" hap2_plusW.fa | \
        seqkit seq -m ${params.min_scaffold_size} -o ${tolid}_hap2.fa

    # ---- Index ----
    samtools faidx ${tolid}_hap1.fa
    samtools faidx ${tolid}_hap2.fa

    # ---- Cleanup intermediates ----
    rm -f hap1_noW.fa hap2_plusW.fa W_scaffolds.fa
    """
}

/*
 * Run BUSCO on a genome
 * Collects: full_table, single-copy FAA, summary JSON
 */
process BUSCO {
    tag "${tolid}_${hap_name}"
    label 'busco'
    publishDir "${params.outdir}/02_busco/${tolid}_${hap_name}", mode: 'copy'

    input:
    tuple val(tolid), val(hap_name), path(genome)

    output:
    tuple val(tolid), val(hap_name), path("full_table_${tolid}_${hap_name}.tsv"),      emit: full_table
    tuple val(tolid), val(hap_name), path("${tolid}_${hap_name}_single_copy.faa"),     emit: single_copy_faa
    tuple val(tolid), val(hap_name), path("summary_${tolid}_${hap_name}.json"),        emit: summary

    script:
    """
    busco -i ${genome} \
        -l ${params.busco_lineage} \
        -m geno \
        -o busco_out/${tolid}_${hap_name} \
        -c ${task.cpus}

    # Collect key outputs with consistent naming
    cp busco_out/${tolid}_${hap_name}/run_${params.busco_lineage}/full_table.tsv full_table_${tolid}_${hap_name}.tsv

    cat busco_out/${tolid}_${hap_name}/run_${params.busco_lineage}/busco_sequences/single_copy_busco_sequences/*.faa \
        > ${tolid}_${hap_name}_single_copy.faa

    cp busco_out/${tolid}_${hap_name}/short_summary.specific.${params.busco_lineage}.${tolid}_${hap_name}.json \
        summary_${tolid}_${hap_name}.json

    tar -czf busco_out/${tolid}_${hap_name}.tar.gz --remove-files busco_out/${tolid}_${hap_name}
    """
}


/*
 * Stage FAA files into an OrthoFinder input directory
 * Renames {tolid}_single_copy.faa -> {tolid}.fasta
 */
process STAGE_ORTHOFINDER_INPUT {
    label 'small'

    input:
    path(faa_files)

    output:
    path("of_input/")

    script:
    """
    mkdir -p of_input
    for f in ${faa_files}; do
        species=\$(basename \$f _hap1_single_copy.faa)
        cp \$f of_input/\${species}.fasta
    done
    """
}

/*
 * Run OrthoFinder
 * If --user_tree is provided, uses it with the -s flag
 * Outputs: full results directory + species trees
 */
process ORTHOFINDER {
    label 'orthofinder'
    publishDir "${params.outdir}/03_orthofinder", mode: 'copy'

    input:
    path(fasta_dir)

    output:
    path("results/"),                        emit: all_results
    path("results/latest/"),                 emit: results
    path("species_tree_rooted.txt"),         emit: species_tree
    path("species_tree_node_labels.txt"),    emit: node_label_tree

    script:
    def tree_opt = params.user_tree ? "-s ${params.user_tree}" : ""
    """
    orthofinder -f ${fasta_dir} -t ${task.cpus} ${tree_opt}

    # Keep all date-stamped results (preserves earlier runs)
    mkdir -p results
    mv ${fasta_dir}/OrthoFinder/Results_* results/

    # Find the latest run (most recently modified Results_* directory)
    LATEST=\$(ls -td results/Results_* | head -1)

    # Copy key outputs to a fixed 'results/latest/' for downstream processes
    mkdir -p results/latest
    cp -r \$LATEST/Orthologues                          results/latest/
    cp -r \$LATEST/Orthogroups                          results/latest/
    cp -r \$LATEST/Species_Tree                         results/latest/
    cp -r \$LATEST/Phylogenetic_Hierarchical_Orthogroups results/latest/
    cp -r \$LATEST/Comparative_Genomics_Statistics      results/latest/

    # Copy species trees for convenient downstream access
    cp results/latest/Species_Tree/SpeciesTree_rooted.txt species_tree_rooted.txt
    cp results/latest/Species_Tree/SpeciesTree_rooted_node_labels.txt species_tree_node_labels.txt

    # cleanup intermediate results (optional, can be commented out to keep all runs)
    rm -rf results/Results_*/WorkingDirectory
    """
}

/*
 * Prepare Syngraph input from OrthoFinder single-copy orthologues
 * Filters FASTA headers to genes in single-copy orthogroups,
 * extracts positional information, and outputs per-species TSVs:
 *   busco_id <TAB> chromosome <TAB> start <TAB> end
 */
process PREP_SYNGRAPH_INPUT {
    label 'small'
    publishDir "${params.outdir}/04_syngraph/input", mode: 'copy'

    input:
    path(orthofinder_results)
    path(fasta_dir)

    output:
    path("sg_input/")

    script:
    """
    prep_syngraph_input.py ${orthofinder_results} ${fasta_dir} sg_input
    """
}


/*
 * Run Syngraph: build, tabulate, and infer at multiple m values
 * Uses the OrthoFinder rooted species tree (without node labels)
 */
process SYNGRAPH {
    label 'syngraph'
    publishDir "${params.outdir}/04_syngraph", mode: 'copy'

    input:
    path(input_dir)
    path(species_tree)

    output:
    path("syngraph_out/"),    emit: results

    script:
    """
    mkdir -p syngraph_out

    # Build synteny graph and tabulate
    ${params.syngraph_path} build -d ${input_dir} -m -o syngraph_out/syngraph_build
    ${params.syngraph_path} tabulate -g syngraph_out/syngraph_build.pickle -o syngraph_out/syngraph_build

    # Infer ancestral genomes at multiple m values
    for m in \$(seq -w ${params.syngraph_m_min} ${params.syngraph_m_max}); do
        mkdir -p syngraph_out/syngraph_m\${m}
        ${params.syngraph_path} infer \
            -g syngraph_out/syngraph_build.pickle \
            -t ${species_tree} \
            -r 2 --use_dist \
            -s ${params.ref_genome} \
            -m \${m} \
            -o syngraph_out/syngraph_m\${m}/infer_m\${m}
        ${params.syngraph_path} tabulate \
            -g syngraph_out/syngraph_m\${m}/infer_m\${m}.with_ancestors.pickle \
            -o syngraph_out/syngraph_m\${m}/infer_m\${m}
    done
    """
}

/*
 * Prepare AGORA input files from OrthoFinder output
 * Creates:
 *   - species_tree.nwk    (with node labels, required by AGORA)
 *   - orthologyGroups/     (one file per phylogenetic node)
 *   - genes/               (one file per species, from FASTA headers)
 */
process PREP_AGORA_INPUT {
    label 'small'
    publishDir "${params.outdir}/05_agora/input", mode: 'copy'

    input:
    path(orthofinder_results)
    path(fasta_dir)

    output:
    path("species_tree.nwk"),      emit: species_tree
    path("orthologyGroups/"),      emit: ortho_groups
    path("genes/"),                emit: genes

    script:
    """
    # Species tree with node labels (required by AGORA)
    cp ${orthofinder_results}/Species_Tree/SpeciesTree_rooted_node_labels.txt species_tree.nwk

    # Create orthologyGroups from Phylogenetic Hierarchical Orthogroups
    # For each node: remove header, filter SCAFF/unloc, keep gene columns (4+),
    # strip leading whitespace, sanitise gene IDs (remove |+, |-, replace |: with _)
    mkdir -p orthologyGroups
    for tsv in ${orthofinder_results}/Phylogenetic_Hierarchical_Orthogroups/*.tsv; do
        node=\$(basename \$tsv .tsv)
        tail -n +2 \$tsv | \
            grep -v -E "SCAFF|unloc" | \
            cut -f4- | \
            sed 's/^[ \\t]*//' | \
            sed 's/|+//g;s/|-//g;s/[|:]/_/g' \
            > orthologyGroups/orthologyGroups.\${node}.list
    done

    # N0 orthogroups from Orthogroups.tsv (newer OrthoFinder 3.x format)
    if [ -f ${orthofinder_results}/Orthogroups/Orthogroups.tsv ]; then
        cut -f2- ${orthofinder_results}/Orthogroups/Orthogroups.tsv | \
            tail -n +2 | \
            grep -v -E "SCAFF|unloc" | \
            sed 's/^[ \\t]*//' | \
            sed 's/|+//g;s/|-//g;s/[|:]/_/g' \
            > orthologyGroups/orthologyGroups.N0.list
    fi

    # Create gene lists from FASTA headers
    mkdir -p genes
    for fasta in ${fasta_dir}/*.fasta; do
        species=\$(basename \$fasta .fasta)
        parse_fasta_headers.py \$fasta genes/genes.\${species}.list
    done
    """
}

/*
 * Run AGORA ancestral genome reconstruction
 */
process AGORA {
    label 'agora'
    publishDir "${params.outdir}/05_agora/output", mode: 'copy'

    input:
    path(species_tree)
    path(ortho_groups)
    path(genes)

    output:
    path("agora_output/"),  emit: results

    script:
    """
    mkdir -p agora_output
    ${params.agora_path} ${species_tree} \
        ${ortho_groups}/orthologyGroups.%s.list \
        ${genes}/genes.%s.list \
        -workingDir=agora_output \
        -nbThreads=${task.cpus}
    """
}

// =====================================================
// Workflow
// =====================================================

workflow {

    // ----- Parse genome list -----
    // Input: text file with one genome path per line (includes both hap1 and hap2)
    // Output channel: [tolid, hap1_path, hap2_path]
    ch_raw = Channel.fromPath(params.genome_list)
        .splitText()
        .map { it.trim() }
        .filter { it }
        .map { path ->
            def fname = file(path).name
            def tolid = fname.split('\\.')[0]
            def hap   = fname.contains('hap1') ? 'hap1' : 'hap2'
            tuple(tolid, hap, file(path))
        }
        .groupTuple()
        .map { tolid, haps, paths ->
            def idx1 = haps.indexOf('hap1')
            def idx2 = haps.indexOf('hap2')
            def hap1 = idx1 >= 0 ? paths[idx1] : null
            def hap2 = idx2 >= 0 ? paths[idx2] : null
            tuple(tolid, hap1, hap2)
        }
        .filter { tolid, hap1, hap2 -> hap1 != null && hap2 != null }

    // ----- Step 1/2: Prepare genomes + BUSCO (or use precomputed BUSCO) -----
    if (params.skip_busco) {
        if (!params.precomputed_busco_dir) {
            error "ERROR: --precomputed_busco_dir is required when --skip_busco is set"
        }
        ch_faa_files = Channel.fromPath("${params.precomputed_busco_dir}/**/*_hap1_single_copy.faa").collect()
        if (!ch_faa_files) {
            error "ERROR: no *_hap1_single_copy.faa files found under --precomputed_busco_dir (${params.precomputed_busco_dir})"
        }
    } else {
        // Step 1: Prepare genomes (move W, filter, index)
        PREPARE_GENOMES(ch_raw)

        // Split into per-haplotype channels for BUSCO
        ch_hap1 = PREPARE_GENOMES.out.map { tolid, h1, h1i, h2, h2i -> tuple(tolid, 'hap1', h1) }
        ch_hap2 = PREPARE_GENOMES.out.map { tolid, h1, h1i, h2, h2i -> tuple(tolid, 'hap2', h2) }

        // Step 2: BUSCO on both haplotypes
        BUSCO(ch_hap1.mix(ch_hap2))
        ch_faa_files = BUSCO.out.single_copy_faa
            .filter { tolid, hap_name, faa -> hap_name == 'hap1' }
            .map { tolid, hap_name, faa -> faa }
            .collect()
    }

    // ----- Step 3: Prepare OrthoFinder input (hap1 only) -----
    STAGE_ORTHOFINDER_INPUT(ch_faa_files)

    // ----- Step 4: OrthoFinder (or use precomputed) -----
    if (params.skip_orthofinder) {
        if (!params.precomputed_orthofinder_dir) {
            error "ERROR: --precomputed_orthofinder_dir is required when --skip_orthofinder is set"
        }
        ch_orthofinder_results  = Channel.fromPath(params.precomputed_orthofinder_dir, type: 'dir')
        ch_species_tree         = Channel.fromPath("${params.precomputed_orthofinder_dir}/Species_Tree/SpeciesTree_rooted.txt")
        ch_node_label_tree      = Channel.fromPath("${params.precomputed_orthofinder_dir}/Species_Tree/SpeciesTree_rooted_node_labels.txt")
    } else {
        ORTHOFINDER(STAGE_ORTHOFINDER_INPUT.out)
        ch_orthofinder_results  = ORTHOFINDER.out.results
        ch_species_tree         = ORTHOFINDER.out.species_tree
        ch_node_label_tree      = ORTHOFINDER.out.node_label_tree
    }

    // ----- Step 5a: Prepare Syngraph input from OrthoFinder SCOs -----
    PREP_SYNGRAPH_INPUT(
        ch_orthofinder_results,
        STAGE_ORTHOFINDER_INPUT.out
    )

    // ----- Step 5b: Syngraph -----
    // Uses the rooted species tree (without node labels)
    SYNGRAPH(
        PREP_SYNGRAPH_INPUT.out,
        ch_species_tree
    )

    // ----- Step 6: Prepare AGORA input -----
    PREP_AGORA_INPUT(
        ch_orthofinder_results,
        STAGE_ORTHOFINDER_INPUT.out
    )

    // ----- Step 7: AGORA -----
    AGORA(
        PREP_AGORA_INPUT.out.species_tree,
        PREP_AGORA_INPUT.out.ortho_groups,
        PREP_AGORA_INPUT.out.genes
    )
}
