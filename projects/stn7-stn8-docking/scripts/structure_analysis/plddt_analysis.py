#!/usr/bin/env python3
"""
STN7 and STN8 AlphaFold pLDDT Confidence Profile Plotter
Extracts B-factors (representing pLDDT scores) from the PDB files,
computes average confidence scores for full structures and kinase domains,
and plots residue-level confidence profiles categorized by standard AlphaFold thresholds.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB

def extract_plddt(pdb_file):
    """
    Extracts residue numbers and pLDDT scores (B-factors) from a PDB file.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]
    chain = model['A']
    
    res_nums = []
    plddts = []
    
    for residue in chain:
        if PDB.Polypeptide.is_aa(residue):
            res_num = residue.get_id()[1]
            # B-factor of any atom in the residue represents pLDDT in AF models
            # We take the B-factor of the CA atom
            try:
                plddt = residue['CA'].get_bfactor()
                res_nums.append(res_num)
                plddts.append(plddt)
            except KeyError:
                pass
                
    return np.array(res_nums), np.array(plddts)

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    stn7_pdb_path = os.path.join(base_dir, "data/structures/Q9S713_AlphaFold.pdb")
    stn8_pdb_path = os.path.join(base_dir, "data/structures/Q9LZV4_AlphaFold.pdb")
    
    figures_dir = os.path.join(base_dir, "figures")
    results_dir = os.path.join(base_dir, "results")
    
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    print("=== Analyzing AlphaFold pLDDT Confidence Scores ===")
    
    # 1. Extract pLDDT scores
    stn7_res, stn7_plddt = extract_plddt(stn7_pdb_path)
    stn8_res, stn8_plddt = extract_plddt(stn8_pdb_path)
    
    # 2. Segment Kinase domains
    # Kinase domains: STN7 (134-452), STN8 (133-477)
    stn7_kin_mask = (stn7_res >= 134) & (stn7_res <= 452)
    stn8_kin_mask = (stn8_res >= 133) & (stn8_res <= 477)
    
    stn7_kin_plddt = stn7_plddt[stn7_kin_mask]
    stn8_kin_plddt = stn8_plddt[stn8_kin_mask]
    
    # 3. Calculate statistics
    stats = {
        "Protein": ["STN7 (Q9S713)", "STN8 (Q9LZV4)"],
        "Total_Residues": [len(stn7_res), len(stn8_res)],
        "Overall_Mean_pLDDT": [np.mean(stn7_plddt), np.mean(stn8_plddt)],
        "Overall_Median_pLDDT": [np.median(stn7_plddt), np.median(stn8_plddt)],
        "Kinase_Mean_pLDDT": [np.mean(stn7_kin_plddt), np.mean(stn8_kin_plddt)],
        "Kinase_Median_pLDDT": [np.median(stn7_kin_plddt), np.median(stn8_kin_plddt)],
        "Very_High_Confidence_Pct (>90)": [np.mean(stn7_plddt > 90) * 100, np.mean(stn8_plddt > 90) * 100],
        "Confident_Pct (70-90)": [np.mean((stn7_plddt >= 70) & (stn7_plddt <= 90)) * 100, np.mean((stn8_plddt >= 70) & (stn8_plddt <= 90)) * 100],
        "Low_Confidence_Pct (50-70)": [np.mean((stn7_plddt >= 50) & (stn7_plddt < 70)) * 100, np.mean((stn8_plddt >= 50) & (stn8_plddt < 70)) * 100],
        "Very_Low_Confidence_Pct (<50)": [np.mean(stn7_plddt < 50) * 100, np.mean(stn8_plddt < 50) * 100]
    }
    
    df_stats = pd.DataFrame(stats)
    stats_csv = os.path.join(results_dir, "plddt_summary.csv")
    df_stats.to_csv(stats_csv, index=False)
    print(f"Saved pLDDT statistics summary to: {stats_csv}")
    print(df_stats.to_string(index=False))
    
    # 4. Generate Confidence Profile Plot
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharey=True)
    
    colors = {
        'very_high': '#1e3a8a',  # Dark Blue (>90)
        'confident': '#3b82f6',  # Light Blue (70-90)
        'low': '#eab308',        # Yellow (50-70)
        'very_low': '#ef4444'    # Red (<50)
    }
    
    # Plot STN7
    axes[0].plot(stn7_res, stn7_plddt, color='#374151', linewidth=1, alpha=0.5)
    # Color code by kinase domain
    axes[0].axvspan(134, 452, color='#14b8a6', alpha=0.15, label='Stromal Kinase Domain (134-452)')
    # Add horizontal threshold lines
    axes[0].axhline(y=90, color='#10b981', linestyle=':', alpha=0.7)
    axes[0].axhline(y=70, color='#3b82f6', linestyle=':', alpha=0.7)
    axes[0].axhline(y=50, color='#ef4444', linestyle=':', alpha=0.7)
    
    # Scatter plot residues by confidence level for visual clarity
    for threshold, color, label in [
        (stn7_plddt > 90, colors['very_high'], 'Very High (>90)'),
        ((stn7_plddt >= 70) & (stn7_plddt <= 90), colors['confident'], 'Confident (70-90)'),
        ((stn7_plddt >= 50) & (stn7_plddt < 70), colors['low'], 'Low (50-70)'),
        (stn7_plddt < 50, colors['very_low'], 'Very Low (<50)')
    ]:
        axes[0].scatter(stn7_res[threshold], stn7_plddt[threshold], color=color, s=4, alpha=0.8)
        
    axes[0].set_title("Arabidopsis thaliana STN7 AlphaFold Confidence Profile", fontsize=12, fontweight='bold')
    axes[0].set_ylabel("pLDDT Score", fontsize=10)
    axes[0].grid(True, linestyle=':', alpha=0.5)
    axes[0].legend(loc='lower left', fontsize=9)
    
    # Plot STN8
    axes[1].plot(stn8_res, stn8_plddt, color='#374151', linewidth=1, alpha=0.5)
    # Color code by kinase domain
    axes[1].axvspan(133, 477, color='#f43f5e', alpha=0.15, label='Stromal Kinase Domain (133-477)')
    # Add horizontal threshold lines
    axes[1].axhline(y=90, color='#10b981', linestyle=':', alpha=0.7)
    axes[1].axhline(y=70, color='#3b82f6', linestyle=':', alpha=0.7)
    axes[1].axhline(y=50, color='#ef4444', linestyle=':', alpha=0.7)
    
    # Scatter residues for STN8
    for threshold, color, label in [
        (stn8_plddt > 90, colors['very_high'], 'Very High (>90)'),
        ((stn8_plddt >= 70) & (stn8_plddt <= 90), colors['confident'], 'Confident (70-90)'),
        ((stn8_plddt >= 50) & (stn8_plddt < 70), colors['low'], 'Low (50-70)'),
        (stn8_plddt < 50, colors['very_low'], 'Very Low (<50)')
    ]:
        axes[1].scatter(stn8_res[threshold], stn8_plddt[threshold], color=color, s=4, alpha=0.8, label=label)
        
    axes[1].set_title("Arabidopsis thaliana STN8 AlphaFold Confidence Profile", fontsize=12, fontweight='bold')
    axes[1].set_xlabel("Amino Acid Position", fontsize=10)
    axes[1].set_ylabel("pLDDT Score", fontsize=10)
    axes[1].grid(True, linestyle=':', alpha=0.5)
    axes[1].legend(loc='lower left', fontsize=9)
    
    plt.tight_layout()
    plddt_plot_path = os.path.join(figures_dir, "plddt_confidence_profile.png")
    plt.savefig(plddt_plot_path, dpi=300)
    plt.close()
    
    print(f"\nSaved dual pLDDT confidence plot to: {plddt_plot_path}")
    print("====================================================")

if __name__ == "__main__":
    main()
