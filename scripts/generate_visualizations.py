#!/usr/bin/env python3
"""
Generate visualizations for phosphorylation data analysis pipeline.

This script creates visualizations matching the percent_phosphosites.ipynb notebook:
1. S, T, Y distribution by tissue (line plot with dots)
2. Total phosphorylation sites per tissue (bar chart with outlier detection)
3. Tissue overlap analysis (sites found in N tissues)
4. Tissue specificity analysis (tissue-specific vs multi-tissue)
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import sys


def get_tissue_columns(species):
    """Get tissue columns based on species."""
    if species.lower() == 'mouse':
        return ['brain', 'brownfat', 'heart', 'kidney', 'liver', 'lung', 
                'pancreas', 'spleen', 'testis']
    elif species.lower() == 'rat':
        return ['all_brain', 'cortex', 'brainstem', 'cerebellum', 'testicle',
                'pancreas', 'stomach', 'liver', 'fat', 'intestine', 'kidney',
                'spleen', 'thymus', 'lung', 'muscle', 'heart', 'blood']
    else:
        raise ValueError(f"Unknown species: {species}")


def calculate_tissue_phosphosite_percentages(df, species):
    """
    Calculate percentage of S, T, Y phosphosites for each tissue type.
    For each tissue, calculate: (S/T/Y count in tissue) / (total phosphosites in tissue) * 100
    """
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
    """
    Create a line plot with dots showing S, T, Y percentages for each tissue
    (matching Cell 6 of the notebook)
    """
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
    
    output_path_pdf = os.path.join(output_dir, f'{species}_tissue_phosphosite_percentages.pdf')
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved tissue phosphosite percentages plot (PDF): {output_path_pdf}")
    
    plt.close()


def compute_tissue_counts(df, tissue_cols):
    """For each tissue, compute total count of rows with tissue > 0 (in thousands)."""
    counts = {}
    for col in tissue_cols:
        if col in df.columns:
            count = (df[col] > 0).sum()
            counts[col] = count / 1000.0  # Convert to thousands
    return counts


def plot_tissue_counts(counts, species, output_dir, figsize=(12, 6)):
    """
    Create bar chart showing total phosphorylation sites per tissue
    """
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
    
    output_path_pdf = os.path.join(output_dir, f'{species}_tissue_total_phospho_counts.pdf')
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved tissue total phospho counts plot (PDF): {output_path_pdf}")
    
    plt.close()
    
    # Print average information
    print(f"\n{species.title()} average total phosphorylation sites per tissue: {average_value:.1f}k")


def calculate_tissue_overlap_data(df, tissue_cols):
    """
    Calculate how many phosphosites are found in 1, 2, 3, etc. tissues
    (matching Cell 10)
    """
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
    """
    Create a bar chart showing phosphosite counts by number of tissues (matching Cell 10)
    """
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
    
    output_path_pdf = os.path.join(output_dir, f'{species}_tissue_overlap_counts.pdf')
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved tissue overlap plot (PDF): {output_path_pdf}")
    
    plt.close()


def calculate_tissue_specific_percentage(df, tissue_cols, species):
    """
    Calculate the percentage of phosphosites that are tissue-specific (found in only one tissue)
    (matching Cell 11)
    """
    # Count how many tissues each phosphosite is found in
    df['tissue_count'] = df[tissue_cols].apply(lambda row: (row > 0).sum(), axis=1)
    
    # Calculate percentages
    total_phosphosites = len(df)
    tissue_specific_count = (df['tissue_count'] == 1).sum()
    tissue_specific_percentage = (tissue_specific_count / total_phosphosites) * 100
    
    # Also calculate multi-tissue phosphosites
    multi_tissue_count = (df['tissue_count'] > 1).sum()
    multi_tissue_percentage = (multi_tissue_count / total_phosphosites) * 100
    
    # Get breakdown by tissue count
    tissue_count_breakdown = df['tissue_count'].value_counts().sort_index()
    
    return {
        'total_phosphosites': total_phosphosites,
        'tissue_specific_count': tissue_specific_count,
        'tissue_specific_percentage': tissue_specific_percentage,
        'multi_tissue_count': multi_tissue_count,
        'multi_tissue_percentage': multi_tissue_percentage,
        'tissue_count_breakdown': tissue_count_breakdown
    }


def plot_tissue_specificity(results, species, output_dir):
    """
    Create a pie chart showing tissue-specific vs multi-tissue phosphosites (matching Cell 11)
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Pie chart for tissue specificity
    labels = ['Tissue-specific\n(1 tissue)', 'Multi-tissue\n(2+ tissues)']
    sizes = [results['tissue_specific_percentage'], results['multi_tissue_percentage']]
    colors = ['#FF6B6B', '#4ECDC4']
    
    ax1.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    ax1.set_title(f'Tissue Specificity - {species.title()}')
    
    # Bar chart showing distribution by tissue count
    tissue_counts = results['tissue_count_breakdown']
    counts_k = tissue_counts.values / 1000.0  # Convert to thousands
    
    bars = ax2.bar(tissue_counts.index.astype(str), counts_k, color='skyblue', alpha=0.8)
    ax2.set_xlabel('Number of Tissues')
    ax2.set_ylabel('Number of Sites (Thousands)')
    ax2.set_title(f'Distribution by Tissue Count - {species.title()}')
    ax2.grid(True, axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax2.annotate(f'{height:.1f}k',
                        xy=(bar.get_x() + bar.get_width()/2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    
    # Save plot
    output_path = os.path.join(output_dir, f'{species}_tissue_specificity.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved tissue specificity plot: {output_path}")
    
    output_path_pdf = os.path.join(output_dir, f'{species}_tissue_specificity.pdf')
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved tissue specificity plot (PDF): {output_path_pdf}")
    
    plt.close()


def calculate_phosphosite_percentages(df):
    """Calculate S, T, Y percentages from a dataframe."""
    if 'amino_acid' not in df.columns:
        return None
    
    df_sty = df[df['amino_acid'].isin(['S', 'T', 'Y'])].copy()
    if len(df_sty) == 0:
        return None
    
    total = len(df_sty)
    return {
        "Ser": (df_sty["amino_acid"] == "S").sum() / total * 100,
        "Thr": (df_sty["amino_acid"] == "T").sum() / total * 100,
        "Tyr": (df_sty["amino_acid"] == "Y").sum() / total * 100
    }


def create_species_comparison_sty_plot(df, species, output_dir, comparison_df=None, comparison_species=None):
    """
    Create a bar chart showing S, T, Y percentages.
    If comparison_df is provided, creates a grouped bar chart comparing both species.
    Otherwise, creates a single-species plot.
    """
    if 'amino_acid' not in df.columns:
        print("Warning: 'amino_acid' column not found. Skipping S/T/Y distribution plot.")
        return
    
    # Calculate percentages for current species
    percentages = calculate_phosphosite_percentages(df)
    if percentages is None:
        print("Warning: No S/T/Y data found. Skipping S/T/Y distribution plot.")
        return
    
    # If comparison data is provided, create grouped bar chart
    if comparison_df is not None and comparison_species is not None:
        comparison_percentages = calculate_phosphosite_percentages(comparison_df)
        if comparison_percentages is None:
            print("Warning: No S/T/Y data in comparison dataset. Creating single-species plot.")
            comparison_df = None
    
    if comparison_df is not None and comparison_species is not None:
        # Create grouped bar chart (matching Cell 3 of notebook)
        # Determine colors based on species (mouse='#e02b35', rat='#f0c571')
        species_color = '#e02b35' if species.lower() == 'mouse' else '#f0c571'
        comparison_color = '#e02b35' if comparison_species.lower() == 'mouse' else '#f0c571'
        
        data = pd.DataFrame({
            "residue": ["Ser", "Thr", "Tyr"],
            species: [percentages["Ser"], percentages["Thr"], percentages["Tyr"]],
            comparison_species: [comparison_percentages["Ser"], comparison_percentages["Thr"], comparison_percentages["Tyr"]]
        })
        
        fig, ax = plt.subplots(figsize=(6, 4))
        bar_width = 0.35
        x = range(len(data["residue"]))
        
        ax.bar([i - bar_width/2 for i in x], data[species], width=bar_width, 
               label=species.title(), color=species_color)
        ax.bar([i + bar_width/2 for i in x], data[comparison_species], width=bar_width, 
               label=comparison_species.title(), color=comparison_color)
        
        ax.set_xticks(x)
        ax.set_xticklabels(data["residue"])
        ax.set_ylabel("Percentage of total phosphosites (%)")
        ax.set_xlabel("Residue type")
        ax.set_title("Percentage of Ser, Thr, and Tyr phosphosites of total phosphosites in Mouse vs Rat")
        ax.legend()
        ax.grid(True, axis='y', alpha=0.3)
        
        plt.tight_layout()
        
        # Save the plot
        output_path = os.path.join(output_dir, 'mouse_rat_sty_comparison.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved S/T/Y comparison plot: {output_path}")
        
        output_path_pdf = os.path.join(output_dir, 'mouse_rat_sty_comparison.pdf')
        plt.savefig(output_path_pdf, bbox_inches='tight')
        print(f"Saved S/T/Y comparison plot (PDF): {output_path_pdf}")
        
        plt.close()
    else:
        # Create single-species plot
        residues = ['Ser', 'Thr', 'Tyr']
        pct_values = [percentages["Ser"], percentages["Thr"], percentages["Tyr"]]
        
        # Use colors matching the notebook: mouse='#e02b35', rat='#f0c571'
        color = '#e02b35' if species.lower() == 'mouse' else '#f0c571'
        
        fig, ax = plt.subplots(figsize=(6, 4))
        bars = ax.bar(residues, pct_values, color=color, alpha=0.8, width=0.6)
        ax.set_ylabel('Percentage of total phosphosites (%)')
        ax.set_xlabel('Residue type')
        ax.set_title(f'Percentage of Ser, Thr, and Tyr phosphosites of total phosphosites in {species.title()}')
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
        
        # Save the plot
        output_path = os.path.join(output_dir, f'{species}_sty_distribution.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved S/T/Y distribution plot: {output_path}")
        
        output_path_pdf = os.path.join(output_dir, f'{species}_sty_distribution.pdf')
        plt.savefig(output_path_pdf, bbox_inches='tight')
        print(f"Saved S/T/Y distribution plot (PDF): {output_path_pdf}")
        
        plt.close()


def main():
    parser = argparse.ArgumentParser(description='Generate visualizations for phosphorylation data')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file (aligned data)')
    parser.add_argument('-s', '--species', required=True, choices=['mouse', 'rat'], help='Species name')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for plots')
    parser.add_argument('--comparison_input', help='Optional: Second CSV file for species comparison (must also specify --comparison_species)')
    parser.add_argument('--comparison_species', choices=['mouse', 'rat'], help='Species name for comparison file')
    
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
    
    # Read comparison file if provided
    comparison_df = None
    if args.comparison_input:
        if not args.comparison_species:
            print("Error: --comparison_species must be specified when using --comparison_input")
            sys.exit(1)
        print(f"Reading comparison file: {args.comparison_input}")
        try:
            comparison_df = pd.read_csv(args.comparison_input)
            print(f"Loaded {len(comparison_df)} rows for comparison")
        except Exception as e:
            print(f"Error reading comparison CSV file: {e}")
            sys.exit(1)
    
    # Get tissue columns
    tissue_columns = get_tissue_columns(args.species)
    available_tissues = [col for col in tissue_columns if col in df.columns]
    
    if len(available_tissues) == 0:
        print("Warning: No tissue columns found in data.")
    
    # Check for required columns
    if 'amino_acid' not in df.columns:
        print("Warning: 'amino_acid' column not found. Some plots may be skipped.")
    
    print("\nGenerating visualizations...")
    
    # 0. S/T/Y distribution plot (matching Cell 3 style)
    if 'amino_acid' in df.columns:
        print("\n0. Creating S/T/Y distribution plot...")
        create_species_comparison_sty_plot(df, args.species, args.output_dir, 
                                          comparison_df=comparison_df, 
                                          comparison_species=args.comparison_species)
    
    # 1. Tissue phosphosite percentages (line plot with dots) - Cell 6
    if 'amino_acid' in df.columns:
        print("\n1. Creating tissue phosphosite percentages plot...")
        tissue_results = calculate_tissue_phosphosite_percentages(df, args.species)
        if tissue_results:
            create_tissue_phosphosite_plot(tissue_results, args.species, args.output_dir, 
                                          figsize=(15, 8) if args.species == 'rat' else (12, 8))
    
    # 2. Total phosphosites per tissue (bar chart with average line) - Cell 8
    if available_tissues:
        print("\n2. Creating total phosphosites per tissue plot...")
        tissue_counts = compute_tissue_counts(df, available_tissues)
        if tissue_counts:
            plot_tissue_counts(tissue_counts, args.species, args.output_dir,
                              figsize=(15, 6) if args.species == 'rat' else (12, 6))
    
    # 3. Tissue overlap analysis - Cell 10
    if available_tissues:
        print("\n3. Creating tissue overlap plot...")
        overlap_counts = calculate_tissue_overlap_data(df, available_tissues)
        plot_tissue_overlap(overlap_counts, args.species, args.output_dir,
                           figsize=(12, 6) if args.species == 'rat' else (10, 6))
    
    # 4. Tissue specificity analysis - Cell 11
    if available_tissues:
        print("\n4. Creating tissue specificity plot...")
        specificity_results = calculate_tissue_specific_percentage(df, available_tissues, args.species)
        plot_tissue_specificity(specificity_results, args.species, args.output_dir)
        
        # Print summary
        print(f"\n=== TISSUE SPECIFICITY ANALYSIS ===")
        print(f"Total phosphosites: {specificity_results['total_phosphosites']:,}")
        print(f"Tissue-specific (1 tissue): {specificity_results['tissue_specific_count']:,} ({specificity_results['tissue_specific_percentage']:.1f}%)")
        print(f"Multi-tissue (2+ tissues): {specificity_results['multi_tissue_count']:,} ({specificity_results['multi_tissue_percentage']:.1f}%)")
    
    print("\nVisualization generation complete!")


if __name__ == '__main__':
    main()
