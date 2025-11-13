/*
 * New Nextflow Pipeline 
 * Create Blast databases and then processes data through scripts
 * Processes multiple species in parallel
 */

// Parameters
params.species_list = ['mouse', 'rat']  // List of species to process
params.outdir = "${projectDir}/results"

// Create Results Directory
new File(params.outdir).mkdirs()

// Clean Data (CSV)
process clean_data {
    publishDir "${params.outdir}/${species}", mode: 'copy'

    input:
    tuple val(species), path(input_csv)

    output:
    tuple val(species), path("cleaned_data.csv"), emit: cleaned_data

    script:
    """
    # Set environment variable for project directory
    export NF_PIPELINE_DIR="${projectDir}"
    
    # Create symlinks to project directories so scripts can use relative paths
    # This allows scripts to use relative paths like 'data/' and 'configs/' 
    # regardless of where the project is located
    if [ ! -e data ]; then
        ln -sf "${projectDir}/data" data
    fi
    if [ ! -e configs ]; then
        ln -sf "${projectDir}/configs" configs
    fi
    
    "${projectDir}/venv/bin/python3" "${projectDir}/scripts/clean_data.py" --input ${input_csv} --output cleaned_data.csv --species ${species}
    """
}

// Paper BLAST 
process paper_blast {
    publishDir "${params.outdir}/${species}", mode: 'copy'

    input:
    tuple val(species), path(cleaned_data)

    output:
    tuple val(species), path("paper_blast_results.csv"), emit: paper_blast_results

    script:
    """
    # Set environment variable for project directory
    export NF_PIPELINE_DIR="${projectDir}"
    
    # Create symlinks to project directories so scripts can use relative paths
    # This allows scripts to use relative paths like 'data/' and 'configs/' 
    # regardless of where the project is located
    if [ ! -e data ]; then
        ln -sf "${projectDir}/data" data
    fi
    if [ ! -e configs ]; then
        ln -sf "${projectDir}/configs" configs
    fi
    
    "${projectDir}/venv/bin/python3" "${projectDir}/scripts/paper_blast.py" --input ${cleaned_data} --output paper_blast_results.csv --species ${species}
    """
}

// Total BLAST
process total_blast {
    publishDir "${params.outdir}/${species}", mode: 'copy'

    input:
    tuple val(species), path(paper_blast_results)

    output:
    tuple val(species), path("total_blast_results.csv"), emit: total_blast_results

    script:
    """
    # Set environment variable for project directory
    export NF_PIPELINE_DIR="${projectDir}"
    
    # Create symlinks to project directories so scripts can use relative paths
    # This allows scripts to use relative paths like 'data/' and 'configs/' 
    # regardless of where the project is located
    if [ ! -e data ]; then
        ln -sf "${projectDir}/data" data
    fi
    if [ ! -e configs ]; then
        ln -sf "${projectDir}/configs" configs
    fi
    
    "${projectDir}/venv/bin/python3" "${projectDir}/scripts/total_blast.py" --input ${paper_blast_results} --output total_blast_results.csv --species ${species}
    """
}

// Align Sequences
process align_sequences {
    publishDir "${params.outdir}/${species}", mode: 'copy'

    input:
    tuple val(species), path(total_blast_results)

    output:
    tuple val(species), path("aligned_sequences.csv"), emit: aligned_sequences

    script:
    """
    # Set environment variable for project directory
    export NF_PIPELINE_DIR="${projectDir}"
    
    # Create symlinks to project directories so scripts can use relative paths
    # This allows scripts to use relative paths like 'data/' and 'configs/' 
    # regardless of where the project is located
    if [ ! -e data ]; then
        ln -sf "${projectDir}/data" data
    fi
    if [ ! -e configs ]; then
        ln -sf "${projectDir}/configs" configs
    fi
    
    "${projectDir}/venv/bin/python3" "${projectDir}/scripts/align_to_full_seq.py" --input ${total_blast_results} --output aligned_sequences.csv --species ${species}
    """
}

// Generate Visualizations
process visualization {
    publishDir "${params.outdir}/${species}", mode: 'copy'

    input:
    tuple val(species), path(aligned_sequences)
    
    output:
    path "*.{png}", emit: visualization

    script:
    """
    # Set environment variable for project directory
    export NF_PIPELINE_DIR="${projectDir}"
    
    # Create symlinks to project directories so scripts can use relative paths
    # This allows scripts to use relative paths like 'data/' and 'configs/' 
    # regardless of where the project is located
    if [ ! -e data ]; then
        ln -sf "${projectDir}/data" data
    fi
    if [ ! -e configs ]; then
        ln -sf "${projectDir}/configs" configs
    fi
    
    "${projectDir}/venv/bin/python3" "${projectDir}/scripts/visualization.py" --input ${aligned_sequences} --output_dir . --species ${species}
    """
}

// Workflow
workflow {
    // Create channel from species list for parallel processing
    species_ch = Channel.from(params.species_list)
    
    // Create input file channels for each species [species, file]
    input_files = species_ch.map { species ->
        tuple(species, file("${projectDir}/data/${species}_test_data.csv"))
    }
    
    // Step 1: Clean data (parallel for each species)
    // Input: tuple(species, file) -> Output: tuple(species, cleaned_data.csv)
    clean_data(input_files)
    
    // Output already includes species, so no need to re-pair
    cleaned_with_species = clean_data.out.cleaned_data
    
    // Step 2: Paper BLAST (parallel for each species)
    paper_blast(cleaned_with_species)
    
    // Output already includes species
    paper_blast_with_species = paper_blast.out.paper_blast_results
    
    // Step 3: Total BLAST (parallel for each species)
    total_blast(paper_blast_with_species)
    
    // Output already includes species
    total_blast_with_species = total_blast.out.total_blast_results
    
    // Step 4: Align sequences (parallel for each species)
    align_sequences(total_blast_with_species)
    
    // Output already includes species
    aligned_with_species = align_sequences.out.aligned_sequences
    
    // Step 5: Generate visualizations
   visualization(aligned_with_species)
}