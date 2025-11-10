#!/usr/bin/env python3
"""
Generate visualizations for phosphorylation data analysis pipeline.

This script creates visualizations showing:
1. Percentage of sites (total, aligned, etc.)
2. Number of proteins
3. S, T, Y distributions
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


def calculate_summary_stats(df, species):
    """Calculate summary statistics from the dataset."""
    stats = {}
    
    # Total phosphosites
    stats['total_phosphosites'] = len(df)
    
    # Successfully aligned phosphosites
    if 'alignment_success' in df.columns:
        stats['aligned_phosphosites'] = df['alignment_success'].sum()
        stats['alignment_percentage'] = (stats['aligned_phosphosites'] / stats['total_phosphosites']) * 100
    else:
        stats['aligned_phosphosites'] = 0
        stats['alignment_percentage'] = 0
    
    # Number of unique proteins
    if 'Protein' in df.columns:
        stats['unique_proteins'] = df['Protein'].nunique()
    elif 'protein_id' in df.columns:
        stats['unique_proteins'] = df['protein_id'].nunique()
    else:
        stats['unique_proteins'] = 0
    
    # S, T, Y counts and percentages
    if 'amino_acid' in df.columns:
        stats['ser_count'] = (df['amino_acid'] == 'S').sum()
        stats['thr_count'] = (df['amino_acid'] == 'T').sum()
        stats['tyr_count'] = (df['amino_acid'] == 'Y').sum()
        
        total_sty = stats['ser_count'] + stats['thr_count'] + stats['tyr_count']
        if total_sty > 0:
            stats['ser_percentage'] = (stats['ser_count'] / total_sty) * 100
            stats['thr_percentage'] = (stats['thr_count'] / total_sty) * 100
            stats['tyr_percentage'] = (stats['tyr_count'] / total_sty) * 100
        else:
            stats['ser_percentage'] = 0
            stats['thr_percentage'] = 0
            stats['tyr_percentage'] = 0
    else:
        stats['ser_count'] = 0
        stats['thr_count'] = 0
        stats['tyr_count'] = 0
        stats['ser_percentage'] = 0
        stats['thr_percentage'] = 0
        stats['tyr_percentage'] = 0
    
    return stats


def create_summary_plot(stats, species, output_dir):
    """Create a summary plot showing key statistics."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Phosphorylation Data Summary - {species.title()}', fontsize=16, fontweight='bold')
    
    # Plot 1: Total sites and alignment percentage
    ax1 = axes[0, 0]
    categories = ['Total\nSites', 'Aligned\nSites', 'Unique\nProteins']
    values = [stats['total_phosphosites'], stats['aligned_phosphosites'], stats['unique_proteins']]
    colors = ['#4C78A8', '#54A24B', '#E45756']
    bars = ax1.bar(categories, values, color=colors, alpha=0.8)
    ax1.set_ylabel('Count')
    ax1.set_title('Total Sites, Aligned Sites, and Unique Proteins')
    ax1.grid(True, axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax1.annotate(f'{int(height):,}',
                    xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=10)
    
    # Add alignment percentage as text
    ax1.text(0.5, 0.95, f"Alignment: {stats['alignment_percentage']:.1f}%",
            transform=ax1.transAxes, ha='center', va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
            fontsize=11, fontweight='bold')
    
    # Plot 2: S, T, Y distribution (pie chart)
    ax2 = axes[0, 1]
    if stats['ser_count'] + stats['thr_count'] + stats['tyr_count'] > 0:
        sizes = [stats['ser_percentage'], stats['thr_percentage'], stats['tyr_percentage']]
        labels = ['Ser', 'Thr', 'Tyr']
        colors_pie = ['#2E8B57', '#FF6347', '#4169E1']
        explode = (0.05, 0.05, 0.05)
        
        wedges, texts, autotexts = ax2.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%',
                                           startangle=90, explode=explode, shadow=True)
        ax2.set_title('S, T, Y Distribution')
        
        # Make percentage text bold
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
    else:
        ax2.text(0.5, 0.5, 'No S/T/Y data available', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_title('S, T, Y Distribution')
    
    # Plot 3: S, T, Y counts (bar chart)
    ax3 = axes[1, 0]
    if stats['ser_count'] + stats['thr_count'] + stats['tyr_count'] > 0:
        residues = ['Ser', 'Thr', 'Tyr']
        counts = [stats['ser_count'], stats['thr_count'], stats['tyr_count']]
        colors_bar = ['#2E8B57', '#FF6347', '#4169E1']
        bars = ax3.bar(residues, counts, color=colors_bar, alpha=0.8)
        ax3.set_ylabel('Count')
        ax3.set_title('S, T, Y Counts')
        ax3.grid(True, axis='y', alpha=0.3)
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax3.annotate(f'{int(height):,}',
                        xy=(bar.get_x() + bar.get_width()/2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=10)
    else:
        ax3.text(0.5, 0.5, 'No S/T/Y data available', ha='center', va='center', transform=ax3.transAxes)
        ax3.set_title('S, T, Y Counts')
    
    # Plot 4: Alignment success rate
    ax4 = axes[1, 1]
    if stats['total_phosphosites'] > 0:
        aligned = stats['aligned_phosphosites']
        failed = stats['total_phosphosites'] - aligned
        sizes = [aligned, failed]
        labels = ['Aligned', 'Failed']
        colors_pie = ['#54A24B', '#E45756']
        
        wedges, texts, autotexts = ax4.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%',
                                           startangle=90, shadow=True)
        ax4.set_title('Alignment Success Rate')
        
        # Make percentage text bold
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
    else:
        ax4.text(0.5, 0.5, 'No alignment data available', ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title('Alignment Success Rate')
    
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(output_dir, f'{species}_summary_statistics.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved summary plot: {output_path}")
    
    output_path_pdf = os.path.join(output_dir, f'{species}_summary_statistics.pdf')
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved summary plot (PDF): {output_path_pdf}")
    
    plt.close()


def create_sty_distribution_plot(df, species, output_dir):
    """Create a detailed S, T, Y distribution plot."""
    if 'amino_acid' not in df.columns:
        print("Warning: 'amino_acid' column not found. Skipping S/T/Y distribution plot.")
        return
    
    # Filter for S, T, Y only
    df_sty = df[df['amino_acid'].isin(['S', 'T', 'Y'])].copy()
    
    if len(df_sty) == 0:
        print("Warning: No S/T/Y data found. Skipping S/T/Y distribution plot.")
        return
    
    # Calculate percentages
    total = len(df_sty)
    ser_pct = (df_sty['amino_acid'] == 'S').sum() / total * 100
    thr_pct = (df_sty['amino_acid'] == 'T').sum() / total * 100
    tyr_pct = (df_sty['amino_acid'] == 'Y').sum() / total * 100
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Create bar chart
    residues = ['Ser', 'Thr', 'Tyr']
    percentages = [ser_pct, thr_pct, tyr_pct]
    colors = ['#2E8B57', '#FF6347', '#4169E1']
    
    bars = ax.bar(residues, percentages, color=colors, alpha=0.8, width=0.6)
    ax.set_ylabel('Percentage of Total Phosphosites (%)', fontsize=12)
    ax.set_xlabel('Residue Type', fontsize=12)
    ax.set_title(f'S, T, Y Distribution - {species.title()}', fontsize=14, fontweight='bold')
    ax.set_ylim([0, max(percentages) * 1.15])
    ax.grid(True, axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar, pct in zip(bars, percentages):
        height = bar.get_height()
        ax.annotate(f'{pct:.1f}%',
                   xy=(bar.get_x() + bar.get_width()/2, height),
                   xytext=(0, 3),
                   textcoords="offset points",
                   ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    # Add total count text
    ax.text(0.02, 0.98, f'Total S/T/Y sites: {total:,}',
           transform=ax.transAxes, ha='left', va='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
           fontsize=10)
    
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(output_dir, f'{species}_sty_distribution.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved S/T/Y distribution plot: {output_path}")
    
    output_path_pdf = os.path.join(output_dir, f'{species}_sty_distribution.pdf')
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved S/T/Y distribution plot (PDF): {output_path_pdf}")
    
    plt.close()


def create_tissue_sty_plot(df, species, output_dir):
    """Create S, T, Y distribution by tissue type."""
    if 'amino_acid' not in df.columns:
        print("Warning: 'amino_acid' column not found. Skipping tissue S/T/Y plot.")
        return
    
    tissue_columns = get_tissue_columns(species)
    available_tissues = [col for col in tissue_columns if col in df.columns]
    
    if len(available_tissues) == 0:
        print("Warning: No tissue columns found. Skipping tissue S/T/Y plot.")
        return
    
    # Filter for S, T, Y only
    df_sty = df[df['amino_acid'].isin(['S', 'T', 'Y'])].copy()
    
    if len(df_sty) == 0:
        print("Warning: No S/T/Y data found. Skipping tissue S/T/Y plot.")
        return
    
    # Calculate percentages for each tissue
    results = []
    for tissue in available_tissues:
        tissue_data = df_sty[df_sty[tissue] > 0]
        if len(tissue_data) > 0:
            total_in_tissue = len(tissue_data)
            ser_pct = (tissue_data['amino_acid'] == 'S').sum() / total_in_tissue * 100
            thr_pct = (tissue_data['amino_acid'] == 'T').sum() / total_in_tissue * 100
            tyr_pct = (tissue_data['amino_acid'] == 'Y').sum() / total_in_tissue * 100
            
            results.append({
                'tissue': tissue,
                'ser': ser_pct,
                'thr': thr_pct,
                'tyr': tyr_pct,
                'total': total_in_tissue
            })
    
    if len(results) == 0:
        print("Warning: No tissue data found. Skipping tissue S/T/Y plot.")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    tissues = [r['tissue'] for r in results]
    ser_percentages = [r['ser'] for r in results]
    thr_percentages = [r['thr'] for r in results]
    tyr_percentages = [r['tyr'] for r in results]
    
    x = np.arange(len(tissues))
    width = 0.25
    
    # Create grouped bar chart
    bars1 = ax.bar(x - width, ser_percentages, width, label='Ser', color='#2E8B57', alpha=0.8)
    bars2 = ax.bar(x, thr_percentages, width, label='Thr', color='#FF6347', alpha=0.8)
    bars3 = ax.bar(x + width, tyr_percentages, width, label='Tyr', color='#4169E1', alpha=0.8)
    
    ax.set_ylabel('Percentage of Phosphosites (%)', fontsize=12)
    ax.set_xlabel('Tissue Type', fontsize=12)
    ax.set_title(f'S, T, Y Distribution by Tissue - {species.title()}', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(tissues, rotation=45, ha='right')
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    
    # Save the plot
    output_path = os.path.join(output_dir, f'{species}_tissue_sty_distribution.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved tissue S/T/Y distribution plot: {output_path}")
    
    output_path_pdf = os.path.join(output_dir, f'{species}_tissue_sty_distribution.pdf')
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved tissue S/T/Y distribution plot (PDF): {output_path_pdf}")
    
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
    
    # Calculate summary statistics
    print("\nCalculating summary statistics...")
    stats = calculate_summary_stats(df, args.species)
    
    # Print summary statistics
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"Total phosphosites: {stats['total_phosphosites']:,}")
    print(f"Successfully aligned: {stats['aligned_phosphosites']:,} ({stats['alignment_percentage']:.2f}%)")
    print(f"Unique proteins: {stats['unique_proteins']:,}")
    print(f"\nS/T/Y Distribution:")
    print(f"  Ser: {stats['ser_count']:,} ({stats['ser_percentage']:.2f}%)")
    print(f"  Thr: {stats['thr_count']:,} ({stats['thr_percentage']:.2f}%)")
    print(f"  Tyr: {stats['tyr_count']:,} ({stats['tyr_percentage']:.2f}%)")
    print("="*60)
    
    # Create visualizations
    print("\nGenerating visualizations...")
    
    # 1. Summary statistics plot
    create_summary_plot(stats, args.species, args.output_dir)
    
    # 2. S, T, Y distribution plot
    create_sty_distribution_plot(df, args.species, args.output_dir)
    
    # 3. Tissue-specific S, T, Y distribution plot
    create_tissue_sty_plot(df, args.species, args.output_dir)
    
    print("\nVisualization generation complete!")


if __name__ == '__main__':
    main()

