/*
 * New Nextflow Pipeline 
 * Create Blast databases and then processes data through scripts
 * Processes multiple species in parallel
 * Portable across Linux and macOS systems
 */

// Parameters (can be overridden in nextflow.config or command line)
params.species_list = ['mouse', 'rat']  // List of species to process
params.outdir = "${projectDir}/results"
params.python_env = 'auto'  // 'auto', 'venv', or 'conda'
params.python_exec = null  // Optional: explicit path to Python executable

// Function to detect Python executable
// Works on both Linux and macOS, supports venv and conda
// This will be evaluated in each process context
def getPythonExec() {
    // If explicitly set, use it
    if (params.python_exec) {
        return params.python_exec
    }
    
    // Try venv first (most common)
    def venvPython3 = new File("${projectDir}/venv/bin/python3")
    def venvPython = new File("${projectDir}/venv/bin/python")
    if (venvPython3.exists()) {
        return venvPython3.toString()
    } else if (venvPython.exists()) {
        return venvPython.toString()
    }
    
    // Try conda environment if specified or auto
    if (params.python_env == 'conda' || params.python_env == 'auto') {
        // Return conda command - will be checked at runtime
        return "conda run -n nf_phospho_pipeline python"
    }
    
    // Fall back to system python3 (works on both Linux and macOS)
    return "python3"
}

// Get Python executable - will be resolved in process scripts
pythonExec = getPythonExec()

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