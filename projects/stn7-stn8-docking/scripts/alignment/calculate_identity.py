#!/usr/bin/env python3
"""
STN7 and STN8 Sequence Identity Matrix Calculator
Reads the combined multiple sequence alignment (all_msa.fasta),
calculates the pairwise sequence identity matrix, saves it as a CSV,
and generates a publication-quality Seaborn heatmap.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

def calculate_identity(seq1, seq2):
    """
    Calculates sequence identity percentage between two aligned sequences.
    Identity % = (identical residues / alignment length) * 100
    We exclude columns where both sequences have gaps.
    """
    identical = 0
    total_len = 0
    
    for c1, c2 in zip(seq1, seq2):
        if c1 == '-' and c2 == '-':
            continue # Exclude columns with gaps in both sequences
        
        total_len += 1
        if c1 == c2:
            identical += 1
            
    if total_len == 0:
        return 0.0
    return (identical / total_len) * 100

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    msa_path = os.path.join(base_dir, "data/msa/all_msa.fasta")
    figures_dir = os.path.join(base_dir, "figures")
    results_dir = os.path.join(base_dir, "results")
    
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    print("=== Calculating Pairwise Sequence Identity ===")
    records = list(SeqIO.parse(msa_path, "fasta"))
    num_seqs = len(records)
    
    if num_seqs == 0:
        print("Error: No aligned sequences found.")
        return
        
    names = [rec.id.replace("STN7_", "STN7:").replace("STN8_", "STN8:") for rec in records]
    
    # Calculate matrix
    matrix = np.zeros((num_seqs, num_seqs))
    for i in range(num_seqs):
        for j in range(num_seqs):
            matrix[i, j] = calculate_identity(str(records[i].seq), str(records[j].seq))
            
    df = pd.DataFrame(matrix, index=names, columns=names)
    
    # Save CSV
    csv_path = os.path.join(results_dir, "sequence_identity_matrix.csv")
    df.to_csv(csv_path)
    print(f"Saved sequence identity matrix to: {csv_path}")
    
    # Generate Heatmap
    plt.figure(figsize=(10, 8))
    sns.set_theme(style="white")
    
    # Harmonies palette: sleek dark-teal theme
    cmap = sns.diverging_palette(220, 170, as_cmap=True)
    
    sns.heatmap(
        df, 
        annot=True, 
        fmt=".1f", 
        cmap="viridis", 
        linewidths=.5, 
        square=True, 
        cbar_kws={"shrink": .75, "label": "Identity Percentage (%)"}
    )
    
    plt.title("STN7 & STN8 Comparative Sequence Identity Matrix", fontsize=14, fontweight="bold", pad=15)
    plt.tight_layout()
    
    heatmap_path = os.path.join(figures_dir, "sequence_identity_heatmap.png")
    plt.savefig(heatmap_path, dpi=300)
    plt.close()
    print(f"Saved identity heatmap to: {heatmap_path}")
    
    # Print basic summary
    print("\n--- Summary of Sequence Identities ---")
    print(f"Average pairwise identity across all homologs: {matrix.mean():.2f}%")
    print(f"Identity between Arabidopsis STN7 and STN8: {df.loc[names[0], names[4]]:.2f}%")
    print(f"Identity between Arabidopsis STN7 and Rice STN7: {df.loc[names[0], names[1]]:.2f}%")
    print(f"Identity between Arabidopsis STN8 and Rice STN8: {df.loc[names[4], names[5]]:.2f}%")
    print("====================================================")

if __name__ == "__main__":
    main()
