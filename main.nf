/*
 * Phosphorylation Data Analysis Pipeline
 * 
 * A Nextflow pipeline for processing phosphorylation site data across multiple species.
 * This pipeline orchestrates data cleaning, BLAST analysis, and sequence alignment.
 */

// Workflow parameters
params.species = 'mouse'
params.input_data = null
params.data_dir = null  // Path to data directory (default: ${projectDir}/data)
params.output_dir = 'results'
params.help = false

// Print help message
if (params.help) {
    log.info """
    Phosphorylation Data Analysis Pipeline
    
    Usage:
        nextflow run main.nf --species <mouse|rat> [options]
    
    Required parameters:
        --species          Species to process (mouse or rat)
        --input_data       Path to input data file (e.g., data/processed/mouse/full_data.csv)
    
    Optional parameters:
        --data_dir         Path to data directory containing blast_dbs, mouse_blast, rat_blast (default: ${projectDir}/data)
        --output_dir       Output directory (default: results)
        --help             Show this help message
    
    Examples:
        # Process mouse data
        nextflow run main.nf --species mouse --input_data data/mouse_full_data.csv
        
        # Process rat data
        nextflow run main.nf --species rat --input_data data/rat_full_data.csv
    """
    exit 0
}

// Validate required parameters
if (!params.species) {
    error "Species parameter is required. Use --species mouse or --species rat"
}
if (params.species != 'mouse' && params.species != 'rat') {
    error "Species must be either 'mouse' or 'rat'"
}

// Set default data directory if not provided (use NF_pipeline/data)
if (!params.data_dir) {
    params.data_dir = "${projectDir}/data"
}

// Set default input data path if not provided
if (!params.input_data) {
    params.input_data = "${params.data_dir}/${params.species}_full_data.csv"
}

// Define output directory structure
results_dir = file(params.output_dir)
cleaned_dir = file("${params.output_dir}/cleaned")
blast_dir = file("${params.output_dir}/blast")
aligned_dir = file("${params.output_dir}/aligned")
visualizations_dir = file("${params.output_dir}/visualizations")

log.info "========================================="
log.info "Phosphorylation Data Analysis Pipeline"
log.info "========================================="
log.info "Species: ${params.species}"
log.info "Input data: ${params.input_data}"
log.info "Data directory: ${params.data_dir}"
log.info "Output directory: ${params.output_dir}"
log.info "========================================="

// Process definitions

process CREATE_BLAST_DB {
    tag "${species}"
    label 'blast_db_creation'
    
    publishDir "${params.output_dir}/blast_dbs", mode: 'copy', pattern: '*.*'
    
    input:
    val species
    path scripts_dir
    path data_dir
    
    output:
    path "${species}_full_blast.*", emit: full_blast_db
    path "${species}_paper_blast.*", emit: paper_blast_db
    
    script:
    """
    # Create blast_dbs directory if it doesn't exist
    mkdir -p \${data_dir}/blast_dbs
    
    # Map species to FASTA files
    if [ "${species}" == "mouse" ]; then
        FULL_FASTA="\${data_dir}/mouse_blast/uniprotkb_taxonomy_id_10090_2025_10_06.fasta"
        PAPER_FASTA="\${data_dir}/mouse_blast/UniprotKB_mouse_paper.fasta"
        FULL_DB="\${data_dir}/blast_dbs/mouse_full_blast"
        PAPER_DB="\${data_dir}/blast_dbs/mouse_paper_blast"
    else
        FULL_FASTA="\${data_dir}/rat_blast/uniprotkb_taxonomy_id_10116_2025_10_06.fasta"
        PAPER_FASTA="\${data_dir}/rat_blast/UniprotKB_rat_paper.fasta"
        FULL_DB="\${data_dir}/blast_dbs/rat_full_blast"
        PAPER_DB="\${data_dir}/blast_dbs/rat_paper_blast"
    fi
    
    # Create full BLAST database if FASTA exists and DB doesn't exist
    if [ -f "\$FULL_FASTA" ] && [ ! -f "\${FULL_DB}.phr" ]; then
        echo "Creating full BLAST database from: \$FULL_FASTA"
        python \${scripts_dir}/create_blast_dbs.py \$FULL_FASTA --db-path \${FULL_DB}
    else
        if [ -f "\${FULL_DB}.phr" ]; then
            echo "Full BLAST database already exists: \${FULL_DB}"
        else
            echo "Warning: Full FASTA file not found: \$FULL_FASTA"
        fi
    fi
    
    # Create paper BLAST database if FASTA exists and DB doesn't exist
    if [ -f "\$PAPER_FASTA" ] && [ ! -f "\${PAPER_DB}.phr" ]; then
        echo "Creating paper BLAST database from: \$PAPER_FASTA"
        python \${scripts_dir}/create_blast_dbs.py \$PAPER_FASTA --db-path \${PAPER_DB}
    else
        if [ -f "\${PAPER_DB}.phr" ]; then
            echo "Paper BLAST database already exists: \${PAPER_DB}"
        else
            echo "Warning: Paper FASTA file not found: \$PAPER_FASTA"
        fi
    fi
    
    # Copy database files to output (for publishing)
    cp \${FULL_DB}.* . 2>/dev/null || true
    cp \${PAPER_DB}.* . 2>/dev/null || true
    """
}

process CLEAN_DATA {
    tag "${species}"
    label 'data_cleaning'
    
    publishDir "${params.output_dir}/cleaned", mode: 'copy', pattern: '*.csv'
    
    input:
    path input_file
    val species
    path scripts_dir
    
    output:
    path "cleaned_${species}_data.csv", emit: cleaned_data
    
    script:
    """
    python \${scripts_dir}/clean_data.py \\
        -i ${input_file} \\
        -s ${species} \\
        -o cleaned_${species}_data.csv
    """
}

process PAPER_BLAST {
    tag "${species}"
    label 'blast_analysis'
    
    publishDir "${params.output_dir}/blast", mode: 'copy', pattern: '*.csv'
    
    input:
    path cleaned_data
    val species
    path scripts_dir
    path data_dir
    
    output:
    path "${species}_paper_blast.csv", emit: paper_blast
    
    script:
    """
    # Create symlink to data directory so scripts can find blast_dbs
    ln -sf \${data_dir}/blast_dbs data/blast_dbs || true
    
    python \${scripts_dir}/paper_blast.py \\
        -i ${cleaned_data} \\
        -s ${species} \\
        -o ${species}_paper_blast.csv
    """
}

process TOTAL_BLAST {
    tag "${species}"
    label 'blast_analysis'
    
    publishDir "${params.output_dir}/blast", mode: 'copy', pattern: '*.csv'
    
    input:
    path paper_blast
    val species
    path scripts_dir
    path data_dir
    
    output:
    path "${species}_total_blast.csv", emit: total_blast
    
    script:
    """
    # Create symlink to data directory so scripts can find blast_dbs
    ln -sf \${data_dir}/blast_dbs data/blast_dbs || true
    
    python \${scripts_dir}/total_blast.py \\
        -i ${paper_blast} \\
        -s ${species} \\
        -o ${species}_total_blast.csv
    """
}

process ALIGN_SEQUENCES {
    tag "${species}"
    label 'sequence_alignment'
    
    publishDir "${params.output_dir}/aligned", mode: 'copy', pattern: '*.csv'
    
    input:
    path total_blast
    val species
    path scripts_dir
    
    output:
    path "${species}_total_blast_aligned.csv", emit: aligned_data
    
    script:
    """
    python \${scripts_dir}/align_to_full_seq.py \\
        -i ${total_blast} \\
        -s ${species} \\
        -o ${species}_total_blast_aligned.csv
    """
}

process GENERATE_VISUALIZATIONS {
    tag "${species}"
    label 'visualization'
    
    publishDir "${params.output_dir}/visualizations", mode: 'copy', pattern: '*.{png,pdf}'
    
    input:
    path aligned_data
    val species
    path scripts_dir
    
    output:
    path "*.{png,pdf}", emit: visualizations
    
    script:
    """
    python \${scripts_dir}/generate_visualizations.py \\
        -i ${aligned_data} \\
        -s ${species} \\
        -o .
    """
}

// Main workflow
workflow {
    // Check if input file exists
    input_file = file(params.input_data)
    if (!input_file.exists()) {
        error "Input file does not exist: ${params.input_data}"
    }
    
    // Check if data directory exists (needed for BLAST databases)
    data_dir = file(params.data_dir)
    if (!data_dir.exists()) {
        error "Data directory does not exist: ${params.data_dir}. Please ensure the data directory contains mouse_blast/, rat_blast/, and blast_dbs/ subdirectories."
    }
    
    // Set up scripts directory (local scripts in NF_pipeline/scripts)
    scripts_dir = file("${projectDir}/scripts")
    if (!scripts_dir.exists()) {
        error "Scripts directory does not exist: ${scripts_dir}. Please ensure scripts are in the scripts/ directory."
    }
    
    // Run the pipeline
    // Step 1: Create BLAST databases from FASTA files
    CREATE_BLAST_DB(params.species, scripts_dir, data_dir)
    
    // Step 2: Clean input data
    CLEAN_DATA(input_file, params.species, scripts_dir)
    
    // Step 3: Paper BLAST (depends on BLAST databases)
    PAPER_BLAST(CLEAN_DATA.out.cleaned_data, params.species, scripts_dir, data_dir)
    
    // Step 4: Total BLAST
    TOTAL_BLAST(PAPER_BLAST.out.paper_blast, params.species, scripts_dir, data_dir)
    
    // Step 5: Align sequences
    ALIGN_SEQUENCES(TOTAL_BLAST.out.total_blast, params.species, scripts_dir)
    
    // Step 6: Generate visualizations
    GENERATE_VISUALIZATIONS(ALIGN_SEQUENCES.out.aligned_data, params.species, scripts_dir)
    
    // Print summary
    ALIGN_SEQUENCES.out.aligned_data
        .view { "Pipeline completed successfully! Final output: ${it}" }
    
    GENERATE_VISUALIZATIONS.out.visualizations
        .view { "Visualization generated: ${it}" }
}

