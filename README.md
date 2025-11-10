# Phosphorylation Data Analysis Pipeline - Nextflow

A standalone Nextflow pipeline for processing phosphorylation site data across multiple species (mouse, rat). This pipeline orchestrates BLAST database creation, data cleaning, BLAST analysis, sequence alignment, and visualization. All scripts and data are self-contained within the `NF_pipeline` directory.

## Overview

This Nextflow pipeline automates the complete phosphorylation data analysis workflow:

1. **BLAST Database Creation** - Create BLAST databases from FASTA files in `mouse_blast/` and `rat_blast/` directories
2. **Data Cleaning** - Filter and clean phosphorylation site data
3. **Paper BLAST** - Map phosphosites using curated paper databases
4. **Total BLAST** - Comprehensive mapping using full UniProtKB databases
5. **Sequence Alignment** - Align phosphosite sequences to full proteins
6. **Visualization** - Generate plots showing percentage of sites, number of proteins, and S/T/Y distributions

## Features

- **Workflow Orchestration**: Automatically manages dependencies between pipeline steps
- **Error Handling**: Built-in retry logic and error reporting
- **Reproducibility**: Tracks all pipeline parameters and versions
- **Parallel Processing**: Can process multiple species simultaneously
- **Resume Capability**: Resume failed runs from the last successful step
- **Comprehensive Reporting**: Generates HTML reports, timelines, and DAGs

## Prerequisites

1. **Nextflow**: Install Nextflow (version 22.04.0 or later)
   ```bash
   curl -s https://get.nextflow.io | bash
   ```

2. **Python Environment**: Ensure Python 3.8+ is installed with required packages
   - The pipeline uses Python scripts from `NF_pipeline/scripts/` (included in this directory)
   - Required packages: pandas, numpy, BioPython, PyYAML, pyxlsb, matplotlib

3. **BLAST+**: Ensure BLAST+ command-line tools are installed and in PATH

4. **Data**: The pipeline is fully standalone and uses data from `NF_pipeline/data/`:
   - Input data files: `data/mouse_full_data.csv` or `data/rat_full_data.csv`
   - FASTA files for full databases: `data/mouse_blast/uniprotkb_taxonomy_id_10090_2025_10_06.fasta` and `data/rat_blast/uniprotkb_taxonomy_id_10116_2025_10_06.fasta`
   - Paper FASTA files: `data/mouse_blast/UniprotKB_mouse_paper.fasta` and `data/rat_blast/UniprotKB_rat_paper.fasta`
   - BLAST databases will be created automatically in `data/blast_dbs/` during pipeline execution
   - The data directory can be specified with `--data_dir` parameter (default: `NF_pipeline/data`)

## Quick Start

### Basic Usage

```bash
# Process mouse data
nextflow run main.nf \
    --species mouse \
    --input_data ../phospho_root/data/processed/mouse/full_data.csv

# Process rat data
nextflow run main.nf \
    --species rat \
    --input_data ../phospho_root/data/processed/rat/full_data.csv
```

### With Custom Paths

```bash
# Specify custom phospho_root path
nextflow run main.nf \
    --species mouse \
    --input_data ../phospho_root/data/processed/mouse/full_data.csv \
    --phospho_root /path/to/phospho_root \
    --output_dir results_mouse
```

## Parameters

### Required Parameters

- `--species`: Species to process (`mouse` or `rat`)
- `--input_data`: Path to input data file (CSV format)

### Optional Parameters

- `--data_dir`: Path to data directory containing `mouse_blast/`, `rat_blast/`, and `blast_dbs/` subdirectories (default: `NF_pipeline/data`)
- `--output_dir`: Output directory for results (default: `results`)
- `--help`: Show help message

## Pipeline Workflow

```
Input Data (CSV)
    ↓
[CREATE_BLAST_DB] (creates databases from FASTA files)
    ↓
[CLEAN_DATA]
    ↓
Cleaned Data
    ↓
[PAPER_BLAST]
    ↓
Paper BLAST Results
    ↓
[TOTAL_BLAST]
    ↓
Total BLAST Results
    ↓
[ALIGN_SEQUENCES]
    ↓
Aligned Sequences
    ↓
[GENERATE_VISUALIZATIONS]
    ↓
Visualization Plots (PNG/PDF)
```

## Output Structure

The pipeline generates the following output structure:

```
results/
├── blast_dbs/
│   ├── {species}_full_blast.*
│   └── {species}_paper_blast.*
├── cleaned/
│   └── cleaned_{species}_data.csv
├── blast/
│   ├── {species}_paper_blast.csv
│   └── {species}_total_blast.csv
├── aligned/
│   └── {species}_total_blast_aligned.csv
├── visualizations/
│   ├── {species}_summary_statistics.png
│   ├── {species}_summary_statistics.pdf
│   ├── {species}_sty_distribution.png
│   ├── {species}_sty_distribution.pdf
│   ├── {species}_tissue_sty_distribution.png
│   └── {species}_tissue_sty_distribution.pdf
├── pipeline_report.html
├── pipeline_timeline.html
├── pipeline_trace.txt
└── pipeline_dag.svg
```

## Process Details

### 1. CREATE_BLAST_DB
- **Purpose**: Create BLAST databases from FASTA files in `mouse_blast/` and `rat_blast/` directories
- **Input**: FASTA files from `data/mouse_blast/` or `data/rat_blast/` (both full and paper databases)
- **Output**: BLAST database files (`.phr`, `.pin`, `.psq`) in `data/blast_dbs/`
- **Resources**: 2 CPUs, 4 GB memory, 2 hours
- **Note**: Databases are created automatically if they don't exist

### 2. CLEAN_DATA
- **Purpose**: Clean and filter phosphorylation data
- **Input**: Raw processed data CSV
- **Output**: Cleaned data CSV
- **Resources**: 1 CPU, 2 GB memory, 1 hour

### 3. PAPER_BLAST
- **Purpose**: BLAST against curated paper databases
- **Input**: Cleaned data CSV
- **Output**: Paper BLAST results CSV
- **Resources**: 4 CPUs, 8 GB memory, 4 hours

### 4. TOTAL_BLAST
- **Purpose**: BLAST against full UniProtKB databases
- **Input**: Paper BLAST results CSV
- **Output**: Total BLAST results CSV
- **Resources**: 4 CPUs, 8 GB memory, 4 hours

### 5. ALIGN_SEQUENCES
- **Purpose**: Align sequences to full proteins
- **Input**: Total BLAST results CSV
- **Output**: Aligned sequences CSV
- **Resources**: 2 CPUs, 4 GB memory, 2 hours

### 6. GENERATE_VISUALIZATIONS
- **Purpose**: Generate visualization plots for analysis results
- **Input**: Aligned sequences CSV
- **Output**: PNG and PDF plots showing:
  - Summary statistics (total sites, aligned sites, unique proteins, alignment percentage)
  - S, T, Y distribution (pie chart and bar chart)
  - S, T, Y distribution by tissue type
- **Resources**: 2 CPUs, 4 GB memory, 1 hour

## Advanced Usage

### Resume Failed Runs

```bash
# Resume from last successful step
nextflow run main.nf -resume \
    --species mouse \
    --input_data ../phospho_root/data/processed/mouse/full_data.csv
```

### Run with Custom Configuration

```bash
# Use custom nextflow.config
nextflow run main.nf -c custom.config \
    --species mouse \
    --input_data ../phospho_root/data/processed/mouse/full_data.csv
```

### Process Multiple Species

```bash
# Process mouse
nextflow run main.nf --species mouse --input_data ../phospho_root/data/processed/mouse/full_data.csv -with-dag results_mouse/dag.svg

# Process rat (in parallel or sequentially)
nextflow run main.nf --species rat --input_data ../phospho_root/data/processed/rat/full_data.csv -with-dag results_rat/dag.svg
```

## Monitoring and Reports

The pipeline automatically generates:

- **HTML Report**: `pipeline_report.html` - Comprehensive execution report
- **Timeline**: `pipeline_timeline.html` - Visual timeline of process execution
- **Trace File**: `pipeline_trace.txt` - Detailed execution trace
- **DAG**: `pipeline_dag.svg` - Workflow dependency graph

View reports:
```bash
# Open HTML report in browser
firefox results/pipeline_report.html
```

## Troubleshooting

### Common Issues

1. **Python Scripts Not Found**
   - All scripts are included in `NF_pipeline/scripts/` directory
   - Ensure the `scripts/` directory exists and contains all required Python scripts

2. **BLAST Databases Missing**
   - The pipeline automatically creates BLAST databases from FASTA files in `mouse_blast/` and `rat_blast/` directories
   - Ensure full FASTA files exist in `data/mouse_blast/` and `data/rat_blast/` (uniprotkb_taxonomy_id_*.fasta)
   - Ensure paper FASTA files exist in `data/mouse_blast/` and `data/rat_blast/` (UniprotKB_mouse_paper.fasta and UniprotKB_rat_paper.fasta)
   - Databases are created automatically during the CREATE_BLAST_DB process

3. **Memory Issues**
   - Adjust memory limits in `nextflow.config`
   - Use `-with-report` to identify memory-intensive processes

4. **Permission Errors**
   - Ensure write permissions for output directory
   - Check Python script execution permissions

### Debug Mode

```bash
# Run with debug output
nextflow run main.nf -with-report -with-trace \
    --species mouse \
    --input_data ../phospho_root/data/processed/mouse/full_data.csv
```

## Integration with Existing Scripts

This pipeline uses Python scripts included in the `NF_pipeline/scripts/` directory:

- `scripts/clean_data.py` - Data cleaning
- `scripts/paper_blast.py` - Paper BLAST analysis
- `scripts/total_blast.py` - Total BLAST analysis
- `scripts/align_to_full_seq.py` - Sequence alignment
- `scripts/generate_visualizations.py` - Visualization generation

All scripts are self-contained in the `NF_pipeline` directory for standalone submission.

The pipeline maintains compatibility with the original shell scripts (`mouse_blast.sh`, `rat_blast.sh`) while adding:

- Automatic dependency management
- Error handling and retry logic
- Resource management
- Comprehensive reporting
- Resume capability

## Configuration

Edit `nextflow.config` to customize:

- Resource limits (CPU, memory, time)
- Error handling strategies
- Executor settings (local, SLURM, SGE, etc.)
- Report generation options

## Example: Complete Workflow

```bash
# 1. Navigate to pipeline directory
cd NF_pipeline

# 2. Run pipeline for mouse
nextflow run main.nf \
    --species mouse \
    --input_data ../phospho_root/data/processed/mouse/full_data.csv \
    --output_dir results_mouse

# 3. Check results
ls -lh results_mouse/

# 4. View report
firefox results_mouse/pipeline_report.html
```

## Future Enhancements

- Support for additional species (human, etc.)
- Parallel processing of multiple input files
- Integration with cloud executors (AWS, GCP)
- Container support (Docker, Singularity)
- Multi-sample processing capabilities

## Citation

If you use this pipeline in your research, please cite:

- Nextflow: Di Tommaso et al. (2017) Nature Biotechnology
- Original phospho_root project: [Your citation]

## License

[Specify license]

## Contact

For questions or issues, please contact: [Your contact information]

