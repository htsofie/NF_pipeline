# Phosphorylation Site Mapping and Analysis Across Multi Tissues and Species - Standard Operating Procedure

## Workflow Rationale

In order to make use of phosphorylation data from older experiments, you need to confront the issue of mapping old protein IDs and their phosphorylation sites to the current updated sets of protein IDs. Furthermore, lots of old data isn't even in Uniprot accession ID form and in a form of ID such as IPI that has long been extinct. In order to map phosphorylation sites to current IDs you have to go through multiple steps: 

1. Data cleaning - Need to remove any sites that are unlocalized, or rows that are not fully complete.
2. Paper BLAST - Search the phosphorylation window against a FASTA file of proteins published in that paper (which can be extracted from uniprot) to identify possible protein matches and record the type of ID they are (either swissprot (reviewed) or tremble (unreviewed)). 
3. Total BLAST - Search any unmatched sequences against a FASTA file of all proteins published in for that species (can download taxonomy specific FASTA files from uniprot).
4. Sequence Alignment and ID matching - Align the phosphorylation sequence window to the full sequence of possible proteins and evaluate match to select best one. This would exact alignments, followed by any swissprot identified match first, then the rest evauluated on distance of new aligned position from old reported position. 
5. Visualization - Generate graphs to visualize distribution of phosphorylation sites between types and tissues.
   
Manual execution of these steps is error-prone, time-consuming, and difficult to reproduce. This Nextflow pipeline automates the complete workflow, ensuring reproducibility, efficiency, accuracy and portability.

## Prerequisites

### Required Software

1. **Nextflow** (version 22.04.0 or later)
   ```bash
   curl -s https://get.nextflow.io | bash
   ```

2. **Java** (version 17 or later)
   - macOS: `brew install openjdk@17`
   - Linux: `sudo apt-get install openjdk-17-jdk`

3. **BLAST+** command-line tools
   - macOS: `brew install blast`
   - Linux: `sudo apt-get install ncbi-blast+`

4. **Python 3.8+** (for virtual environment)

### Required Data Files

The pipeline expects the following data structure:

```
data/
├── mouse_test_data.csv
├── rat_test_data.csv
├── mouse_paper_blast/
│   └── UniprotKB_mouse_paper.fasta
├── mouse_total_blast/
│   └── uniprotkb_taxonomy_id_10090_2025_10_06.fasta
├── rat_paper_blast/
│   └── UniprotKB_rat_paper.fasta
└── rat_total_blast/
    └── uniprotkb_taxonomy_id_10116_2025_10_06.fasta
```
These files are already included in the package so that you can test the pipeline. THe rat_test_data.csv and mouse_test_data.csv are smaller versions of what you could use normally so that running the test pipeline does not take too long. 

File descriptions so you can put your own/new files in:
- **rat_test_data.csv** or **mouse_test_data.csv**: A CSV file of experimental phosphorylation data with at least these column titles: Protein, site_motif, position, amino_acid (for S,T,Y visualization plots), Tissue columns (for tissue-specific plots).
-  **UniprotKB_mouse_paper.fasta** and **UniprotKB_rat_paper.fasta**: in the directory of {species}_paper_blast must contain one fasta file downloaded from uniprot from the experimental paper the phosphorylations were published from.
-  **uniprotkb_taxonomy_id_10090_2025_10_06.fasta** and **uniprotkb_taxonomy_id_10116_2025_10_06.fasta**: in the directory of {species}_total_balst must contain one fasta file downloaded from uniprot containing all proteins found for the species (can search by taxonomy and download).

## Standard Operating Procedure

### Step 1: Environment Setup

Run the setup script to create the Python virtual environment and install all dependencies:

```bash
chmod +x setup_environment.sh
./setup_environment.sh
# To activate environment:
source venv/bin/activate
```

### Step 2: Run the Pipeline

Execute the pipeline with the default settings (processes both mouse and rat):

```bash
./nextflow run main.nf
```

**No additional parameters are required** - the pipeline automatically:
- Detects input files from `data/{species}_test_data.csv`
- Processes both species in parallel
- Creates BLAST databases as needed
- Outputs results to `results/{species}/`

### Step 3: Verify Results

Check that all output files are generated:

```bash
# Check mouse results
ls results/mouse/
# Expected files:
# - cleaned_data.csv
# - paper_blast_results.csv
# - total_blast_results.csv
# - aligned_sequences.csv
# - mouse_sty_distribution.png
# - mouse_tissue_overlap_counts.png
# - mouse_tissue_phosphosite_percentages.png
# - mouse_tissue_total_phospho_counts.png

# Check rat results
ls results/rat/
# Expected files (same structure as mouse)
```

## Expected Results

### Output Files

For each species (`mouse` and `rat`), the pipeline generates:

1. **cleaned_data.csv**: Filtered and standardized phosphorylation data
   - Contains columns: Protein, cleaned_site_motif, motif_position, swissprot_id, trembl_id, ensembl

2. **paper_blast_results.csv**: BLAST results against curated paper databases
   - Contains columns: ID_matches, ID_types, length, full_sequence, match_method

3. **total_blast_results.csv**: BLAST results against full UniProtKB databases
   - Same structure as paper_blast_results.csv

4. **aligned_sequences.csv**: Phosphosites aligned to full protein sequences
   - Contains alignment positions and context sequences

5. **Visualization files** (PNG format):
   - `{species}_sty_distribution.png`: Bar chart of Ser/Thr/Tyr percentages
   - `{species}_tissue_overlap_counts.png`: Tissue overlap analysis
   - `{species}_tissue_phosphosite_percentages.png`: Phosphosite percentages by tissue
   - `{species}_tissue_total_phospho_counts.png`: Total phosphosite counts by tissue

## Pipeline Workflow

```
Input Data (CSV)
    ↓
[clean_data] → Cleaned data with standardized IDs
    ↓
[paper_blast] → BLAST against curated databases
    ↓
[total_blast] → BLAST against full UniProtKB
    ↓
[align_sequences] → Align phosphosites to full proteins
    ↓
[visualization] → Generate analysis plots
    ↓
Output Files (CSV + PNG)
```

## Additional Information

- **Pipeline version**: Nextflow 25.10.0
- **Python packages**: See `requirements.txt`
- **Configuration**: Default settings in `main.nf` work without modification
- **Resume capability**: Use `-resume` flag to continue from last successful step
- **Reports**: HTML reports generated in `results/` directory