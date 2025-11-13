# Generate visualizations that summarize phosphorylation data

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import sys


def get_tissue_columns(species):
    if species.lower() == 'mouse':
        return ['brain', 'brownfat', 'heart', 'kidney', 'liver', 'lung', 
                'pancreas', 'spleen', 'testis']
    elif species.lower() == 'rat':
        return ['all_brain', 'testicle',
                'pancreas', 'stomach', 'liver', 'fat', 'intestine', 'kidney',
                'spleen', 'thymus', 'lung', 'muscle', 'heart', 'blood']
    else:
        raise ValueError(f"Unknown species: {species}")

def calculate_overall_sty_percentages(df):
    if 'amino_acid' not in df.columns:
        return None
    
    # Filter to only S, T, Y phosphosites
    df_sty = df[df['amino_acid'].isin(['S', 'T', 'Y'])].copy()
    if len(df_sty) == 0:
        return None
    
    total = len(df_sty)
    return {
        "Ser": (df_sty["amino_acid"] == "S").sum() / total * 100,
        "Thr": (df_sty["amino_acid"] == "T").sum() / total * 100,
        "Tyr": (df_sty["amino_acid"] == "Y").sum() / total * 100
    }


def plot_sty_bar_graph(df, species, output_dir, figsize=(12, 8)):
    percentages = calculate_overall_sty_percentages(df)
    if percentages is None:
        print("Warning: No S/T/Y data found. Cannot create plot.")
        return
    
    residues = ['Ser', 'Thr', 'Tyr']
    pct_values = [percentages["Ser"], percentages["Thr"], percentages["Tyr"]]
    
    species_colors = {
        'mouse': '#e02b35',
        'rat': '#f0c571'
    }

    # Get color for additional species
    species_lower = species.lower()
    if species_lower not in species_colors: 
        additional_colors = ['#2E8B57', '#FF6347', '#4169E1']
        color_idx = hash(species) % len(additional_colors)
        species_colors[species_lower] = additional_colors[color_idx]
    
    color = species_colors[species_lower]

    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar(residues, pct_values, color=color, alpha=0.8, width=0.6)
    ax.set_ylabel('Percentage of total phosphosites (%)')
    ax.set_xlabel('Residue type')
    ax.set_title(f'Percentage of Ser, Thr, and Tyr phosphosites in {species.title()}')
    ax.set_ylim([0, max(pct_values) * 1.15])
    ax.grid(True, axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar, pct in zip(bars, pct_values):
        height = bar.get_height()
        ax.annotate(f'{pct:.1f}%',
                   xy=(bar.get_x() + bar.get_width()/2, height),
                   xytext=(0, 3),
                   textcoords="offset points",
                   ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, f'{species}_sty_distribution.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved S/T/Y distribution plot: {output_path}")
    
    plt.close()

def calculate_tissue_phosphosite_percentages(df, species):
   
    tissue_columns = get_tissue_columns(species)
    results = {}
    
    for tissue in tissue_columns:
        if tissue in df.columns:
            # Filter rows where this tissue has a value > 0 (phosphosite found in this tissue)
            tissue_data = df[df[tissue] > 0]
            
            if len(tissue_data) > 0:
                # Calculate percentages for S, T, Y in this tissue
                total_in_tissue = len(tissue_data)
                s_count = (tissue_data['amino_acid'] == 'S').sum()
                t_count = (tissue_data['amino_acid'] == 'T').sum()
                y_count = (tissue_data['amino_acid'] == 'Y').sum()
                
                results[tissue] = {
                    'Ser': (s_count / total_in_tissue) * 100,
                    'Thr': (t_count / total_in_tissue) * 100,
                    'Tyr': (y_count / total_in_tissue) * 100,
                    'total_count': total_in_tissue
                }
    
    return results


def create_tissue_phosphosite_plot(results, species, output_dir, figsize=(12, 8)):
   
    if not results:
        print("Warning: No tissue data found. Skipping tissue phosphosite plot.")
        return
    
    # Prepare data for plotting
    tissues = list(results.keys())
    ser_percentages = [results[tissue]['Ser'] for tissue in tissues]
    thr_percentages = [results[tissue]['Thr'] for tissue in tissues]
    tyr_percentages = [results[tissue]['Tyr'] for tissue in tissues]
    
    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)
    
    x = np.arange(len(tissues))
    
    # Create lines with dots
    ax.plot(x, ser_percentages, 'o-', label='Ser', color='#2E8B57', linewidth=2, markersize=6, alpha=0.8)
    ax.plot(x, thr_percentages, 's-', label='Thr', color='#FF6347', linewidth=2, markersize=6, alpha=0.8)
    ax.plot(x, tyr_percentages, '^-', label='Tyr', color='#4169E1', linewidth=2, markersize=6, alpha=0.8)
    
    # Customize the plot
    ax.set_xlabel('Tissue Type')
    ax.set_ylabel('Percentage of Phosphosites (%)')
    ax.set_title(f'Percentage of S, T, Y Phosphosites by Tissue Type - {species.title()}')
    ax.set_xticks(x)
    ax.set_xticklabels(tissues, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Set y-axis to start from 0 for better visualization
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(output_dir, f'{species}_tissue_phosphosite_percentages.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved tissue phosphosite percentages plot: {output_path}")
    
    plt.close()

def compute_tissue_counts(df, tissue_cols):
    counts = {}
    for col in tissue_cols:
        if col in df.columns:
            count = (df[col] > 0).sum()
            counts[col] = count / 1000.0  # Convert to thousands
    return counts


def plot_tissue_counts(counts, species, output_dir, figsize=(12, 6)):
    if not counts:
        print("Warning: No tissue counts found. Skipping tissue counts plot.")
        return
    
    tissues = list(counts.keys())
    values = [counts[t] for t in tissues]

    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar(tissues, values, color='#4C78A8')
    ax.set_ylabel('Total phosphorylation sites (thousands)')
    ax.set_xlabel('Tissue type')
    ax.set_title(f'Total phosphorylation sites per tissue - {species.title()}')
    ax.set_xticks(range(len(tissues)))
    ax.set_xticklabels(tissues, rotation=45, ha='right')
    ax.grid(True, axis='y', alpha=0.3)

    # Calculate average total number of phosphorylation sites across all tissues
    average_value = np.mean(values)
    
    # Add horizontal line for average
    ax.axhline(y=average_value, color='red', linestyle='--', linewidth=2, 
               label=f'Average ({average_value:.1f}k)')

    # Add labels on all bars
    for b in bars:
        h = b.get_height()
        ax.annotate(f'{h:.1f}k', (b.get_x() + b.get_width()/2, h),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=8)

    # Add legend
    ax.legend(loc='upper right')

    plt.tight_layout()
    
    # Save plot
    output_path = os.path.join(output_dir, f'{species}_tissue_total_phospho_counts.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved tissue total phospho counts plot: {output_path}")
    
    plt.close()
    
    # Print average information
    print(f"\n{species.title()} average total phosphorylation sites per tissue: {average_value:.1f}k")


def calculate_tissue_overlap_data(df, tissue_cols):
    """Calculate how many phosphosites are found in 1, 2, 3, etc. tissues"""
    # Count how many tissues each phosphosite is found in
    df['tissue_count'] = df[tissue_cols].apply(lambda row: (row > 0).sum(), axis=1)
    
    # Count phosphosites by number of tissues they're found in
    tissue_overlap_counts = df['tissue_count'].value_counts().sort_index()
    
    # Fill in missing tissue counts with 0
    max_tissues = len(tissue_cols)
    full_counts = pd.Series(0, index=range(1, max_tissues + 1))
    full_counts.update(tissue_overlap_counts)
    
    return full_counts


def plot_tissue_overlap(counts, species, output_dir, figsize=(10, 6)):
    if counts.sum() == 0:
        print("Warning: No tissue overlap data found. Skipping tissue overlap plot.")
        return
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Convert to thousands
    counts_k = counts / 1000.0
    
    # Create bar chart
    bars = ax.bar(counts.index.astype(str), counts_k, color='skyblue', alpha=0.8)
    
    # Customize plot
    ax.set_xlabel('Number of Tissues')
    ax.set_ylabel('Number of Sites (Thousands)')
    ax.set_title(f'Phosphorylation Sites Found in N Tissues - {species.title()}')
    ax.grid(True, axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax.annotate(f'{height:.1f}k',
                       xy=(bar.get_x() + bar.get_width()/2, height),
                       xytext=(0, 3),
                       textcoords="offset points",
                       ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    
    # Save plot
    output_path = os.path.join(output_dir, f'{species}_tissue_overlap_counts.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved tissue overlap plot: {output_path}")
    
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Generate visualizations for phosphorylation data')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file (aligned data)')
    parser.add_argument('-s', '--species', required=True, choices=['mouse', 'rat'], help='Species name')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for plots')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read input CSV
    print(f"Reading input file: {args.input}")
    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        sys.exit(1)
    
    print(f"Loaded {len(df)} rows")
    
    # Generate visualizations
    print("Generating visualizations...")
    tissue_columns = get_tissue_columns(args.species)
    available_tissues = [col for col in tissue_columns if col in df.columns]
    if len(available_tissues) == 0:
        print("Warning: No tissue columns found in data.")
    
    # 1. S/T/Y distribution plot
    if 'amino_acid' in df.columns:
        print("\n1. Creating S/T/Y distribution plot...")
        plot_sty_bar_graph(df, args.species, args.output_dir)
    
    # 2. Tissue phosphosite percentages plot
    if 'amino_acid' in df.columns and available_tissues:
        print("\n2. Creating tissue phosphosite percentages plot...")
        tissue_results = calculate_tissue_phosphosite_percentages(df, args.species)
        if tissue_results:
            create_tissue_phosphosite_plot(tissue_results, args.species, args.output_dir)
    
    # 3. Total phosphosites per tissue plot
    if available_tissues:
        print("\n3. Creating total phosphosites per tissue plot...")
        tissue_counts = compute_tissue_counts(df, available_tissues)
        if tissue_counts:
            plot_tissue_counts(tissue_counts, args.species, args.output_dir)
    
    # 4. Tissue overlap plot
    if available_tissues:
        print("\n4. Creating tissue overlap plot...")
        overlap_counts = calculate_tissue_overlap_data(df, available_tissues)
        plot_tissue_overlap(overlap_counts, args.species, args.output_dir)
    
    print("\nVisualization generation complete!")
if __name__ == '__main__':
    main() 