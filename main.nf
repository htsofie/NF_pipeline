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
    python3 ${projectDir}/scripts/clean_data.py --input ${input_csv} --output cleaned_data.csv --species ${species}
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
    # Create symlinks to data directories so scripts can find FASTA files and blast_dbs
    mkdir -p data
    ln -sf ${projectDir}/data/blast_dbs data/blast_dbs || true
    ln -sf ${projectDir}/data/mouse_paper_blast data/mouse_paper_blast || true
    ln -sf ${projectDir}/data/rat_paper_blast data/rat_paper_blast || true
    
    python3 ${projectDir}/scripts/paper_blast.py --input ${cleaned_data} --output paper_blast_results.csv --species ${species}
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
    # Create symlinks to data directories so scripts can find FASTA files and blast_dbs
    mkdir -p data
    ln -sf ${projectDir}/data/blast_dbs data/blast_dbs || true
    ln -sf ${projectDir}/data/mouse_total_blast data/mouse_total_blast || true
    ln -sf ${projectDir}/data/rat_total_blast data/rat_total_blast || true
    
    python3 ${projectDir}/scripts/total_blast.py --input ${paper_blast_results} --output total_blast_results.csv --species ${species}
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
    python3 ${projectDir}/scripts/align_to_full_seq.py --input ${total_blast_results} --output aligned_sequences.csv --species ${species}
    """
}

// Generate Visualizations (single species)
process generate_visualizations {
    publishDir "${params.outdir}/${species}", mode: 'copy'

    input:
    tuple val(species), path(aligned_sequences)
    
    output:
    path "*.{png,pdf}", emit: visualizations

    script:
    """
    python3 ${projectDir}/scripts/generate_visualizations.py --input ${aligned_sequences} --output_dir . --species ${species}
    """
}

// Generate Comparison Visualizations (both species)
process generate_comparison_visualizations {
    publishDir "${params.outdir}/comparison", mode: 'copy'

    input:
    path mouse_aligned, stageAs: "mouse_aligned.csv"
    path rat_aligned, stageAs: "rat_aligned.csv"
    
    output:
    path "*.{png,pdf}", emit: comparison_visualizations

    script:
    """
    # Generate comparison plot with both species
    python3 ${projectDir}/scripts/generate_visualizations.py \\
        --input mouse_aligned.csv \\
        --species mouse \\
        --comparison_input rat_aligned.csv \\
        --comparison_species rat \\
        --output_dir .
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
    
    // Step 5: Generate individual visualizations (parallel for each species)
    generate_visualizations(aligned_with_species)
    
    // Step 6: Generate comparison visualizations (after both species complete)
    // aligned_with_species has tuples of (species, file)
    // Filter to get mouse and rat separately
    mouse_aligned = aligned_with_species
        .filter { species, file -> species == 'mouse' }
        .map { species, file -> file }
    
    rat_aligned = aligned_with_species
        .filter { species, file -> species == 'rat' }
        .map { species, file -> file }
    
    generate_comparison_visualizations(mouse_aligned, rat_aligned)
}