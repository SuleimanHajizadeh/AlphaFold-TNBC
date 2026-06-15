#!/usr/bin/env python3
"""
STN7 and STN8 Residue-by-Species Conservation Heatmap Generator
Parses the multiple sequence alignment (MSA) and extracts alignment columns
corresponding to the key kinase functional motifs:
  1. G-loop (P-loop)
  2. Catalytic base loop (HRD region)
  3. Activation loop/metal binding (DFG region)
  4. Activation loop/substrate binding (APE region)
Visualizes conservation patterns across 8 homologous sequences.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from Bio import SeqIO

# Physicochemical property classification and color mapping for amino acids
AA_CATEGORIES = {
    'Acidic (D, E)': {'residues': 'DE', 'color': '#ef4444'},       # Red
    'Basic (K, R, H)': {'residues': 'KRH', 'color': '#3b82f6'},     # Blue
    'Polar (S, T, N, Q, Y, C)': {'residues': 'STNQYC', 'color': '#10b981'}, # Green
    'Hydrophobic (A, V, I, L, M, F, W, P)': {'residues': 'AVILMFWP', 'color': '#f59e0b'}, # Orange/Amber
    'Glycine (G)': {'residues': 'G', 'color': '#8b5cf6'},           # Purple
    'Gap (-)': {'residues': '-', 'color': '#e2e8f0'}               # Light Gray
}

def get_aa_color(aa):
    """
    Returns the hex color code for a given amino acid.
    """
    aa = aa.upper()
    for cat_name, cat_data in AA_CATEGORIES.items():
        if aa in cat_data['residues']:
            return cat_data['color']
    return '#ffffff' # Default white

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    msa_path = os.path.join(base_dir, "data/msa/all_msa.fasta")
    figures_dir = os.path.join(base_dir, "figures")
    results_dir = os.path.join(base_dir, "results")
    
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    print("=== Generating Kinase Motif Conservation Heatmap ===")
    
    # 1. Parse MSA records
    records = list(SeqIO.parse(msa_path, "fasta"))
    print(f"Loaded {len(records)} sequences from alignment.")
    
    # Alignments positions (0-indexed columns in all_msa.fasta)
    motifs = {
        "G-loop": (344, 349, "ATP-Binding G-loop"),
        "HRD region": (608, 613, "Catalytic HRD Region"),
        "DFG region": (632, 637, "Mg2+-binding DFG Region"),
        "APE region": (658, 663, "Substrate-binding APE Region")
    }
    
    # Extract headers/species names
    species_labels = []
    for rec in records:
        # Simplify ID for plot labels
        parts = rec.id.split('_')
        # E.g. STN7_Arabidopsis_thaliana_Q9S713 -> STN7 (Arath)
        gene = parts[0]
        species = parts[1][:5]  # Arath, Oryza, Zeama, Chlam
        uniprot = parts[-1]
        label = f"{gene} ({species} - {uniprot})"
        species_labels.append(label)
        
    # Write alignment data to text file
    align_report_path = os.path.join(results_dir, "kinase_motifs_alignment.txt")
    with open(align_report_path, "w") as f:
        f.write("=================================================================\n")
        f.write("STN7/STN8 Kinase Motifs Comparative Sequence Alignment\n")
        f.write("=================================================================\n\n")
        for m_name, (start, end, desc) in motifs.items():
            f.write(f"--- {m_name} ({desc}) columns {start+1}-{end+1} ---\n")
            for i, rec in enumerate(records):
                seq_seg = str(rec.seq)[start:end+1]
                f.write(f"  {species_labels[i]:<30} : {seq_seg}\n")
            f.write("\n")
            
    print(f"Saved motif alignment report to: {align_report_path}")
    
    # ---- Plotting Heatmap ----
    # Create a 4-panel horizontal subplot figure
    fig, axes = plt.subplots(1, 4, figsize=(16, 7), sharey=True)
    fig.suptitle("STN7 & STN8 Kinase Motif Sequence Conservation Analysis", fontsize=15, fontweight='bold', y=0.97)
    
    for idx, (m_name, (start, end, desc)) in enumerate(motifs.items()):
        ax = axes[idx]
        
        # Build numerical matrix for heatmap color indices
        width = end - start + 1
        height = len(records)
        color_matrix = np.zeros((height, width, 3)) # RGB grid
        
        # Extract segments and populate RGB matrix
        for row in range(height):
            seq_seg = str(records[row].seq)[start:end+1]
            for col in range(width):
                aa = seq_seg[col]
                hex_color = get_aa_color(aa)
                # Convert hex to RGB values between 0 and 1
                rgb = [int(hex_color[i:i+2], 16)/255.0 for i in (1, 3, 5)]
                color_matrix[row, col] = rgb
                
        # Draw heatmap using imshow
        ax.imshow(color_matrix, aspect='equal')
        
        # Add labels inside each cell
        for r in range(height):
            seq_seg = str(records[r].seq)[start:end+1]
            for c in range(width):
                aa = seq_seg[c]
                # Determine text color (white for colored cells, dark gray for gaps)
                t_color = '#ffffff' if aa != '-' else '#475569'
                ax.text(c, r, aa, ha='center', va='center', fontsize=11, fontweight='bold', color=t_color)
                
        # Subplot configurations
        ax.set_title(m_name, fontsize=12, fontweight='bold', pad=10)
        ax.set_xticks(np.arange(width))
        ax.set_xticklabels(np.arange(start + 1, end + 2), fontsize=8, rotation=45)
        ax.set_yticks(np.arange(height))
        if idx == 0:
            ax.set_yticklabels(species_labels, fontsize=10, fontweight='bold')
        else:
            ax.set_yticklabels([])
            
        ax.set_xticks(np.arange(width) - 0.5, minor=True)
        ax.set_yticks(np.arange(height) - 0.5, minor=True)
        ax.grid(which="minor", color="white", linestyle='-', linewidth=1.5)
        ax.tick_params(which="minor", size=0)
        
        # Add description below each panel
        ax.set_xlabel(desc, fontsize=9.5, fontstyle='italic', labelpad=8)
        
        # Remove borders
        for spine in ax.spines.values():
            spine.set_visible(False)
            
    # Add category legend
    legend_patches = []
    for cat_name, cat_data in AA_CATEGORIES.items():
        legend_patches.append(mpatches.Patch(color=cat_data['color'], label=cat_name))
        
    fig.legend(handles=legend_patches, loc='lower center', ncol=6, 
               fontsize=10.5, framealpha=0.9, bbox_to_anchor=(0.5, -0.06))
               
    plt.tight_layout()
    plot_path = os.path.join(figures_dir, "kinase_motifs_conservation_heatmap.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved motif conservation heatmap to: {plot_path}")
    print("=====================================================")

if __name__ == "__main__":
    main()
