#!/usr/bin/env nextflow

/**
 * Breakpoint Detection Pipeline
 * 
 * Identifies conserved blocks between two query genomes using direct pairwise alignment,
 * detects lineage-specific breakpoints, and generates visualizations.
 *
 * Workflow:
 * 1. Extract PSL alignments: query1 vs query2, query1 vs ancestor, query2 vs ancestor
 * 2. Identify conserved blocks from direct query1 vs query2 alignment
 * 3. Detect breakpoints in both genomes
 * 4. Plot conserved blocks and breakpoints
 */

nextflow.enable.dsl = 2

params.hal_file = null
params.genome1 = null
params.genome2 = null
params.ancestor = null
params.outdir = "results"
params.min_block_size = 50000
params.min_scaffold_length = 1000000
params.min_alignment_length = 100000

// Validate inputs
if (!params.hal_file || !params.genome1 || !params.genome2 || !params.ancestor) {
    error """
    Missing required parameters.
    
    Usage:
    nextflow run main.nf \\
        --hal_file alignment.hal \\
        --genome1 species_A \\
        --genome2 species_B \\
        --ancestor ancestor_name \\
        --outdir results \\
        --min_block_size 50000
    
    Required parameters:
        --hal_file  : HAL alignment file
        --genome1   : First genome name
        --genome2   : Second genome name
        --ancestor  : Ancestor genome name
    
    Optional parameters:
        --min_block_size : Minimum block size in bp [Default: 50000]
        --outdir         : Output directory [Default: results]
    """
}

// ============================================================================
// Process 1: Extract PSL alignments using halSynteny
// ============================================================================
process extract_psl_alignments {
    publishDir "${params.outdir}/01_psl", mode: 'copy'
    
    input:
    val hal_file
    tuple val(genome_query), val(genome_target)
    
    output:
    path "${genome_query}_vs_${genome_target}.psl"
    
    shell:
    '''
    echo "Extracting alignment: !{genome_query} vs !{genome_target}..."
    halSynteny --queryGenome !{genome_query} --targetGenome !{genome_target} !{hal_file} !{genome_query}_vs_!{genome_target}.psl
    echo "PSL extraction complete: !{genome_query} vs !{genome_target}"
    '''
}

// ============================================================================
// Process 2: Convert query1 vs query2 PSL to alignment blocks
// ============================================================================
process psl_to_alignment_blocks {
    publishDir "${params.outdir}/02_blocks", mode: 'copy'
    
    input:
    path psl_file
    
    output:
    path "alignment_blocks.tsv"
    
    script:
    """
    python3 << 'EOF'
import pandas as pd
import sys

def parse_psl_to_blocks(psl_file):
    '''Convert PSL to alignment blocks TSV format'''
    records = []
    
    with open('${psl_file}', 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            fields = line.split('\\t')
            if len(fields) < 21:
                continue
            
            # PSL format columns
            # 0:matches, 1:mismatches, 8:strand, 9:qName, 11:qStart, 12:qEnd,
            # 13:tName, 15:tStart, 16:tEnd
            q_chr = fields[9]
            q_start = int(fields[11])
            q_end = int(fields[12])
            t_chr = fields[13]
            t_start = int(fields[15])
            t_end = int(fields[16])
            strand = fields[8]
            
            record = {
                'query1_chr': q_chr,
                'query1_start': q_start,
                'query1_end': q_end,
                'query2_chr': t_chr,
                'query2_start': t_start,
                'query2_end': t_end,
                'strand': strand
            }
            records.append(record)
    
    df = pd.DataFrame(records)
    return df

# Parse PSL
df = parse_psl_to_blocks('${psl_file}')
print(f"Parsed {len(df)} alignment blocks from PSL", file=sys.stderr)

# Filter out small scaffolds (query1 and query2)
min_scaffold_length = ${params.min_scaffold_length}
query1_scaffold_lengths = df.groupby('query1_chr')[['query1_start', 'query1_end']].apply(
    lambda x: max(x['query1_end'].max(), x['query1_start'].max())
)
query2_scaffold_lengths = df.groupby('query2_chr')[['query2_start', 'query2_end']].apply(
    lambda x: max(x['query2_end'].max(), x['query2_start'].max())
)

valid_query1_chrs = query1_scaffold_lengths[query1_scaffold_lengths >= min_scaffold_length].index
valid_query2_chrs = query2_scaffold_lengths[query2_scaffold_lengths >= min_scaffold_length].index

df_filtered = df[(df['query1_chr'].isin(valid_query1_chrs)) & (df['query2_chr'].isin(valid_query2_chrs))]

# remove short alignements
min_alignment_length = ${params.min_alignment_length}
df_filtered_2 = df_filtered[(df_filtered['query1_end'] - df_filtered['query1_start'] >= min_alignment_length) &
                          (df_filtered['query2_end'] - df_filtered['query2_start'] >= min_alignment_length)]

print(f"After filtering scaffolds < {min_scaffold_length}bp: {len(df_filtered)} blocks retained", file=sys.stderr)
print(f"Removed {len(df) - len(df_filtered)} blocks from small scaffolds", file=sys.stderr)
print(f"After filtering alignments < {min_alignment_length}bp: {len(df_filtered_2)} blocks retained", file=sys.stderr)
print(f"Removed {len(df_filtered) - len(df_filtered_2)} blocks from short alignments", file=sys.stderr)

# Write output
df_filtered_2.to_csv('alignment_blocks.tsv', sep='\\t', index=False)
print(f"Wrote alignment_blocks.tsv with {len(df_filtered_2)} blocks", file=sys.stderr)
EOF
    """
}

// ============================================================================
// Process 3: Identify conserved blocks from direct alignment
// ============================================================================
process identify_conserved_blocks {
    publishDir "${params.outdir}/03_conserved_blocks", mode: 'copy'
    
    input:
    path blocks_tsv
    path identify_script
    
    output:
    path "conserved_blocks_by_*.tsv", emit: conserved_blocks
    path "query1_breakpoints.tsv", emit: unshared_bp_g1
    path "query2_breakpoints.tsv", emit: unshared_bp_g2
    
    script:
    """
    python3 ${identify_script} --blocks ${blocks_tsv}
    """
}

// ============================================================================
// Process 4: Classify breakpoints using ancestor-based comparisons
// ============================================================================
process classify_breakpoints_with_ancestors {
    publishDir "${params.outdir}/04_classified_breakpoints", mode: 'copy'
    
    input:
    path psl_g1_anc
    path psl_g2_anc
    val genome1
    val genome2
    val ancestor
    path classify_script
    
    output:
    path "${genome1}_breakpoints_classified.tsv", emit: g1_classified
    path "${genome2}_breakpoints_classified.tsv", emit: g2_classified
    
    script:
    """
    python3 ${classify_script} \
        --psl_g1 ${psl_g1_anc} \
        --psl_g2 ${psl_g2_anc} \
        --genome1 ${genome1} \
        --genome2 ${genome2} \
        --ancestor ${ancestor} \
        --min_block_size ${params.min_block_size}
    """
}

// ============================================================================
// Process 5: plot alignemnts and breakpoints
// ============================================================================

process plot_results {
    publishDir "${params.outdir}/06_plots", mode: 'copy'
    
    input:
    path conserved_blocks
    path classified_bp_g1
    path classified_bp_g2
    val genome1
    val genome2
    path plot_script
    
    output:
    path "*.pdf"
    
    script:
    """
    mkdir -p blocks_data
    cp ${conserved_blocks} blocks_data/ || true
    
    python3 ${plot_script} \
        --conserved_blocks blocks_data \
        --breakpoints_g1 ${classified_bp_g1} \
        --breakpoints_g2 ${classified_bp_g2} \
        --name_g1 ${genome1} \
        --name_g2 ${genome2} \
        --output conserved_blocks_and_breakpoints.pdf
    """
}

// ============================================================================
// Main Workflow
// ============================================================================
workflow {
    // Step 1: Create channel with alignment pairs and extract PSL files in parallel
    alignment_pairs = Channel.of(
        [params.genome1, params.genome2],  // direct alignment
        [params.genome1, params.ancestor], // genome1 vs ancestor
        [params.genome2, params.ancestor]  // genome2 vs ancestor
    )
    
    psl_files = extract_psl_alignments(
        params.hal_file,
        alignment_pairs
    )
    
    // Separate the outputs by pair for downstream processes
    psl_direct = psl_files.filter { it.name.contains("${params.genome1}_vs_${params.genome2}") }
    psl_g1_anc = psl_files.filter { it.name.contains("${params.genome1}_vs_${params.ancestor}") }
    psl_g2_anc = psl_files.filter { it.name.contains("${params.genome2}_vs_${params.ancestor}") }
    
    // Step 2: Convert direct alignment PSL to blocks
    alignment_blocks = psl_to_alignment_blocks(psl_direct)
    
    // Step 3: Identify conserved blocks
    identify_script = file("${workflow.projectDir}/bin/identify_conserved_blocks.py")
    conserved = identify_conserved_blocks(
        alignment_blocks,
        identify_script
    )
    
    // Step 4: Classify breakpoints using ancestor-based comparisons
    classify_script = file("${workflow.projectDir}/bin/detect_ancestral_breakpoints.py")
    classified = classify_breakpoints_with_ancestors(
        psl_g1_anc,
        psl_g2_anc,
        params.genome1,
        params.genome2,
        params.ancestor,
        classify_script
    )
    
    
    // Step 5: Plot results with classified breakpoints
    plot_script = file("${workflow.projectDir}/bin/plot_results.py")
    plot_results(
        conserved.conserved_blocks.collect(),
        classified.g1_classified,
        classified.g2_classified,
        params.genome1,
        params.genome2,
        plot_script
    )
}

// ============================================================================
// Workflow Summary
// ============================================================================
workflow.onComplete {
    log.info """
    ==========================================
    Breakpoint Detection Pipeline Complete
    ==========================================
    Status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }
    Duration: $workflow.duration
    
    Results in: ${params.outdir}
      - 01_psl/                   : PSL alignments from halSynteny
      - 02_blocks/                : Alignment blocks in TSV format
      - 03_conserved_blocks/      : Conserved blocks and unshared breakpoints
      - 04_classified_breakpoints/: Ancestry-based breakpoint classification
      - 05_plots/                 : Visualizations of conserved blocks and breakpoints
    
    """.stripIndent()
}

workflow.onError {
    log.error """
    ==========================================
    Breakpoint Detection Pipeline Failed
    ==========================================
    Error: $workflow.errorMessage
    """.stripIndent()
}
