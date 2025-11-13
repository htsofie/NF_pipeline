/*
 * New Nextflow Pipeline 
 * Create Blast databases and then processes data through scripts
 * Processes multiple species in parallel
 * Portable across Linux and macOS systems
 */

// Parameters (can be overridden in nextflow.config or command line)
params.species_list = ['mouse', 'rat']  // List of species to process
params.outdir = "${projectDir}/results"

pythonExec = "${projectDir}/venv/bin/python3"

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
    "${pythonExec}" "${projectDir}/scripts/clean_data.py" \
        --input ${input_csv} \
        --output cleaned_data.csv \
        --species ${species} \
        --project_dir "${projectDir}"
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
    "${pythonExec}" "${projectDir}/scripts/paper_blast.py" \
        --input ${cleaned_data} \
        --output paper_blast_results.csv \
        --species ${species} \
        --project_dir "${projectDir}"
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
    "${pythonExec}" "${projectDir}/scripts/total_blast.py" \
        --input ${paper_blast_results} \
        --output total_blast_results.csv \
        --species ${species} \
        --project_dir "${projectDir}"
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
    "${pythonExec}" "${projectDir}/scripts/align_to_full_seq.py" \
        --input ${total_blast_results} \
        --output aligned_sequences.csv \
        --species ${species} \
        --project_dir "${projectDir}"
    """
}

// Generate Visualizations
process visualization {
    publishDir "${params.outdir}/${species}", mode: 'copy'

    input:
    tuple val(species), path(aligned_sequences)
    
    output:
    path "*.png", emit: visualization

    script:
    """
    "${pythonExec}" "${projectDir}/scripts/visualization.py" \
        --input ${aligned_sequences} \
        --output_dir . \
        --species ${species} \
        --project_dir "${projectDir}"
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
    
    clean_data(input_files)
    
    cleaned_with_species = clean_data.out.cleaned_data
    
    paper_blast(cleaned_with_species)
    
    paper_blast_with_species = paper_blast.out.paper_blast_results

    total_blast(paper_blast_with_species)
    
    total_blast_with_species = total_blast.out.total_blast_results
    
    align_sequences(total_blast_with_species)
    
    aligned_with_species = align_sequences.out.aligned_sequences
    
   visualization(aligned_with_species)
}