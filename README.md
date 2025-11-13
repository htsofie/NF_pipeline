# Phosphorylation Data Analysis Pipeline - Nextflow

A standalone Nextflow pipeline for processing phosphorylation site data across multiple species (mouse, rat). This pipeline orchestrates BLAST database creation, data cleaning, BLAST analysis, sequence alignment, and visualization.

## Overview

This Nextflow pipeline automates the complete phosphorylation data analysis workflow:

1. **BLAST Database Creation** - Create BLAST databases from FASTA files in `mouse_blast/` and `rat_blast/` directories
2. **Data Cleaning** - Filter and clean phosphorylation site data
3. **Paper BLAST** - Map phosphosites using curated paper databases
4. **Total BLAST** - Comprehensive mapping using full UniProtKB databases
5. **Sequence Alignment** - Align phosphosite sequences to full proteins
6. **Visualization** - Generate plots showing percentage of sites, number of proteins, and S/T/Y distributions

## Features

- **Cross-Platform Compatibility**: Works on both Linux and macOS systems
- **Flexible Python Environment**: Supports venv, conda, or system Python
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

2. **Python Environment**: Python 3.8+ is required
   - The pipeline uses Python scripts from `NF_pipeline/scripts/` (included in this directory)
   - See [Environment Setup](#environment-setup) section below for detailed setup instructions

3. **BLAST+**: Ensure BLAST+ command-line tools are installed and in PATH
   - macOS: `brew install blast`
   - Ubuntu/Debian: `sudo apt-get install ncbi-blast+`
   - Or download from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

4. **Data**: The pipeline is fully standalone and uses data from `NF_pipeline/data/`:
   - Input data files: `data/mouse_full_data.csv` or `data/rat_full_data.csv`
   - FASTA files for full databases: `data/mouse_blast/uniprotkb_taxonomy_id_10090_2025_10_06.fasta` and `data/rat_blast/uniprotkb_taxonomy_id_10116_2025_10_06.fasta`
   - Paper FASTA files: `data/mouse_blast/UniprotKB_mouse_paper.fasta` and `data/rat_blast/UniprotKB_rat_paper.fasta`
   - BLAST databases will be created automatically in `data/blast_dbs/` during pipeline execution
   - The data directory can be specified with `--data_dir` parameter (default: `NF_pipeline/data`)

## Portability

This pipeline is designed to be portable across different systems:

- **Linux and macOS**: Fully tested and supported on both platforms
- **Python Environment Detection**: Automatically detects and uses:
  - Virtual environment (`venv/`) if present
  - Conda environment (`nf_phospho_pipeline`) if available
  - System Python as fallback
- **Path Handling**: Uses relative paths and `projectDir` for portability
- **No Hardcoded Paths**: All paths are relative to the project directory

### Python Environment Options

The pipeline supports multiple Python environment setups:

1. **Virtual Environment (venv)** - Recommended for most users
   - Created automatically by `setup_environment.sh`
   - Located at `venv/` in the project directory
   - Works identically on Linux and macOS

2. **Conda Environment** - Alternative option
   - Use `environment.yml` to create: `conda env create -f environment.yml`
   - Environment name: `nf_phospho_pipeline`
   - Specify with: `--python_env conda` or let auto-detection handle it

3. **System Python** - Fallback option
   - Uses system `python3` if no venv/conda found
   - Requires packages to be installed system-wide

## Environment Setup

### Quick Setup (Recommended)

Use the provided setup script to create a virtual environment with all required packages:

```bash
# Make the script executable (if not already)
chmod +x setup_environment.sh

# Run the setup script
./setup_environment.sh
```

This will:
- Detect your operating system (Linux or macOS)
- Optionally create a Conda environment if Conda is available
- Create a Python virtual environment (`venv/`) if not using Conda
- Handle broken/incomplete virtual environments automatically
- Install all required packages with specific versions
- Test the installation
- Check for BLAST+ installation

**Note**: On Linux systems, if you encounter venv creation errors, you may need to install the `python3-venv` package:
```bash
# Ubuntu/Debian
sudo apt install python3.12-venv  # Replace 3.12 with your Python version
```

### Manual Setup

If you prefer to set up manually:

```bash
# Create virtual environment
python3 -m venv venv

# Activate virtual environment
source venv/bin/activate  # On Linux/macOS
# or
venv\Scripts\activate  # On Windows

# Upgrade pip
pip install --upgrade pip setuptools wheel

# Install requirements
pip install -r requirements.txt
```

### Required Python Packages

The pipeline requires the following Python packages (see `requirements.txt` for specific versions):

- **pandas** (>=2.0.0): Data manipulation and analysis
- **numpy** (>=1.24.0): Numerical computing
- **biopython** (>=1.80): Bioinformatics tools for sequence analysis
- **PyYAML** (>=6.0.0): Configuration file parsing
- **matplotlib** (>=3.7.0): Visualization and plotting

### Using the Environment

After setup, activate the virtual environment before running the pipeline:

```bash
# Activate virtual environment
source venv/bin/activate

# Run the pipeline
nextflow run main.nf --species mouse --input_data data/mouse_full_data.csv

# Deactivate when done
deactivate
```

### Verifying Installation

Test that all packages are installed correctly:

```bash
python3 -c "
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import yaml
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Blast import NCBIXML
print('✓ All packages imported successfully!')
"
```

## Quick Start

### 1. Set Up Environment

First, set up the Python environment:

```bash
# Option 1: Use the setup script (recommended)
./setup_environment.sh

# Option 2: Manual setup
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Option 3: Use Conda
conda env create -f environment.yml
conda activate nf_phospho_pipeline
```

### 2. Run the Pipeline

### Basic Usage

```bash
# Activate virtual environment (if using venv)
source venv/bin/activate

# Process mouse data
nextflow run main.nf \
    --species mouse \
    --input_data data/mouse_full_data.csv

# Process rat data
nextflow run main.nf \
    --species rat \
    --input_data data/rat_full_data.csv
```

### With Custom Paths

```bash
# Specify custom data directory
nextflow run main.nf \
    --species mouse \
    --input_data data/mouse_full_data.csv \
    --data_dir /path/to/data \
    --output_dir results_mouse
```

## Parameters

### Required Parameters

- `--species_list`: List of species to process (default: `['mouse', 'rat']`)
- Input files are automatically detected from `data/{species}_test_data.csv`

### Optional Parameters

- `--outdir`: Output directory for results (default: `results`)
- `--python_env`: Python environment type - `auto`, `venv`, or `conda` (default: `auto`)
- `--python_exec`: Explicit path to Python executable (overrides auto-detection)
- `--help`: Show help message

### Python Environment Parameters

The pipeline automatically detects your Python environment, but you can override this:

```bash
# Force use of conda environment
nextflow run main.nf --python_env conda --species_list mouse

# Specify explicit Python path
nextflow run main.nf --python_exec /path/to/python --species_list mouse

# Use auto-detection (default)
nextflow run main.nf --python_env auto --species_list mouse
```

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

2. **Python Environment Not Detected**
   - The pipeline auto-detects Python environments, but if it fails:
     - Check that `venv/bin/python3` exists (if using venv)
     - Check that conda environment exists: `conda env list`
     - Specify explicitly: `--python_exec /path/to/python`
   - On Linux, if venv creation fails, install: `sudo apt install python3-venv`

3. **BLAST Databases Missing**
   - The pipeline automatically creates BLAST databases from FASTA files in `mouse_blast/` and `rat_blast/` directories
   - Ensure full FASTA files exist in `data/mouse_blast/` and `data/rat_blast/` (uniprotkb_taxonomy_id_*.fasta)
   - Ensure paper FASTA files exist in `data/mouse_blast/` and `data/rat_blast/` (UniprotKB_mouse_paper.fasta and UniprotKB_rat_paper.fasta)
   - Databases are created automatically during the CREATE_BLAST_DB process

4. **Memory Issues**
   - Adjust memory limits in `nextflow.config`
   - Use `-with-report` to identify memory-intensive processes

5. **Permission Errors**
   - Ensure write permissions for output directory
   - Check Python script execution permissions

6. **Cross-Platform Issues**
   - If moving between Linux and macOS, ensure:
     - Python environment is recreated on the new system
     - File permissions are correct (scripts should be executable)
     - Line endings are Unix-style (LF, not CRLF)

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

# 2. Set up environment (first time only)
./setup_environment.sh
source venv/bin/activate

# 3. Run pipeline for mouse
nextflow run main.nf \
    --species mouse \
    --input_data data/mouse_full_data.csv \
    --output_dir results_mouse

# 4. Check results
ls -lh results_mouse/

# 5. View report
firefox results_mouse/pipeline_report.html

# 6. Deactivate environment when done
deactivate
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

