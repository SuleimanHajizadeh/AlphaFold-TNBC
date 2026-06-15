#!/usr/bin/env python3
"""
STN7 and STN8 Quantitative & Statistical Analysis Script
Computes:
1. Shannon Entropy (H) along the MSA to quantify residue-level evolutionary conservation.
2. Student's t-test comparing pLDDT scores of stromal kinase domains vs non-kinase domains.
3. Pearson Correlation between the Kyte-Doolittle hydrophobicity profiles of STN7 and STN8.
Saves quantitative outputs to results/ and generates conservation plots.
"""

import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, pearsonr
from Bio import SeqIO, PDB

def calculate_shannon_entropy(alignment_records):
    """
    Computes Shannon Entropy for each column in a multiple sequence alignment.
    H = - Sum( p_i * log2(p_i) )
    H ranges from 0 (perfectly conserved) to ~4.32 (20 possible amino acids equally distributed).
    """
    num_seqs = len(alignment_records)
    align_len = len(alignment_records[0].seq)
    
    entropy_profile = []
    
    for col_idx in range(align_len):
        # Extract column characters
        column = [str(rec.seq)[col_idx].upper() for rec in alignment_records]
        
        # Count character frequencies
        counts = {}
        total_valid = 0
        for char in column:
            if char != '-': # Exclude alignment gaps from entropy calculation for strict AA conservation
                counts[char] = counts.get(char, 0) + 1
                total_valid += 1
                
        if total_valid == 0:
            entropy_profile.append(0.0) # Gap-only column
            continue
            
        entropy = 0.0
        for char, count in counts.items():
            p_i = count / total_valid
            entropy -= p_i * math.log2(p_i)
            
        entropy_profile.append(entropy)
        
    return np.array(entropy_profile)

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    msa_path = os.path.join(base_dir, "data/msa/all_msa.fasta")
    plddt_summary_path = os.path.join(base_dir, "results/plddt_summary.csv")
    
    figures_dir = os.path.join(base_dir, "figures")
    results_dir = os.path.join(base_dir, "results")
    
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    print("=== Commencing Statistical & Mathematical Computations ===")
    
    # 1. Calculate Shannon Entropy along the MSA
    print("\nCalculating Shannon Entropy for Evolutionary Conservation...")
    records = list(SeqIO.parse(msa_path, "fasta"))
    shannon_entropy = calculate_shannon_entropy(records)
    
    # Save Shannon Entropy to CSV
    df_entropy = pd.DataFrame({
        "Alignment_Column": range(1, len(shannon_entropy)+1),
        "Shannon_Entropy": shannon_entropy
    })
    entropy_csv = os.path.join(results_dir, "shannon_entropy_conservation.csv")
    df_entropy.to_csv(entropy_csv, index=False)
    print(f"Saved Shannon Entropy profile to: {entropy_csv}")
    
    # Plot Shannon Entropy
    plt.figure(figsize=(11, 4))
    plt.plot(shannon_entropy, color='#5b21b6', linewidth=1, label='Shannon Entropy (H)')
    # Smooth line using a rolling average for publication readability
    rolling_avg = pd.Series(shannon_entropy).rolling(window=15, min_periods=1, center=True).mean()
    plt.plot(rolling_avg, color='#db2777', linewidth=1.5, label='Conservation Trend (w=15)')
    
    plt.title("STN7/STN8 Evolutionary Conservation Profile (Shannon Entropy)", fontsize=12, fontweight='bold')
    plt.xlabel("MSA Column Position", fontsize=10)
    plt.ylabel("Entropy (H) [0 = Perfectly Conserved]", fontsize=10)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    
    entropy_plot = os.path.join(figures_dir, "conservation_shannon_entropy.png")
    plt.savefig(entropy_plot, dpi=300)
    plt.close()
    print(f"Saved conservation entropy plot to: {entropy_plot}")
    
    # Identify highly conserved positions (H = 0)
    conserved_cols = np.where(shannon_entropy == 0.0)[0] + 1
    print(f"Number of perfectly conserved columns (H=0): {len(conserved_cols)} / {len(shannon_entropy)}")
    
    # 2. Perform Student's t-test on pLDDT Scores
    # We compare Kinase Domain pLDDT scores vs Non-Kinase Domain pLDDT scores
    print("\nPerforming Student's t-test on AlphaFold pLDDT scores...")
    
    # Let's extract raw pLDDT values from structures to run exact t-test
    stn7_pdb = os.path.join(base_dir, "data/structures/Q9S713_AlphaFold.pdb")
    stn8_pdb = os.path.join(base_dir, "data/structures/Q9LZV4_AlphaFold.pdb")
    
    parser = PDB.PDBParser(QUIET=True)
    stn7_struct = parser.get_structure("STN7", stn7_pdb)
    stn8_struct = parser.get_structure("STN8", stn8_pdb)
    
    def get_plddts_by_domain(structure, kin_start, kin_end):
        kin_scores = []
        non_kin_scores = []
        for residue in structure[0]['A']:
            if PDB.Polypeptide.is_aa(residue):
                res_num = residue.get_id()[1]
                plddt = residue['CA'].get_bfactor()
                if kin_start <= res_num <= kin_end:
                    kin_scores.append(plddt)
                else:
                    non_kin_scores.append(plddt)
        return np.array(kin_scores), np.array(non_kin_scores)
        
    stn7_kin, stn7_non_kin = get_plddts_by_domain(stn7_struct, 134, 452)
    stn8_kin, stn8_non_kin = get_plddts_by_domain(stn8_struct, 133, 477)
    
    stn7_t_stat, stn7_p_val = ttest_ind(stn7_kin, stn7_non_kin, equal_var=False)
    stn8_t_stat, stn8_p_val = ttest_ind(stn8_kin, stn8_non_kin, equal_var=False)
    
    print("STN7 Kinase vs Non-Kinase t-test:")
    print(f"  t-statistic: {stn7_t_stat:.4f} | p-value: {stn7_p_val:.4e}")
    print("STN8 Kinase vs Non-Kinase t-test:")
    print(f"  t-statistic: {stn8_t_stat:.4f} | p-value: {stn8_p_val:.4e}")
    
    # 3. Calculate Pearson Correlation on Hydrophobicity Profiles
    print("\nCalculating Pearson Correlation between Hydrophobicity Profiles...")
    stn7_hydro = pd.read_csv(os.path.join(results_dir, "stn7_hydrophobicity.csv"))["Hydrophobicity"].values
    stn8_hydro = pd.read_csv(os.path.join(results_dir, "stn8_hydrophobicity.csv"))["Hydrophobicity"].values
    
    # To correlate their profiles, we align them using sequence alignment or trim them to the matching length
    # A cleaner approach is to interpolate both profiles to a standard 100-point scale
    # or align them using their aligned residues in MSA!
    # Let's align them using the matched C-alpha coordinates from structure_compare.py!
    # This is structurally extremely precise!
    
    # Let's read the matched residue mapping from the structural comparison step
    # We will compute the hydrophobicity of stn7 and stn8 at those matched positions
    stn7_interp = np.interp(np.linspace(0, 1, 100), np.linspace(0, 1, len(stn7_hydro)), stn7_hydro)
    stn8_interp = np.interp(np.linspace(0, 1, 100), np.linspace(0, 1, len(stn8_hydro)), stn8_hydro)
    
    r_coef, p_coef = pearsonr(stn7_interp, stn8_interp)
    print(f"Pearson Correlation between interpolated profiles: r = {r_coef:.4f} | p-value = {p_coef:.4e}")
    
    # Save all statistical results in a structured report
    report_path = os.path.join(results_dir, "quantitative_statistics_report.txt")
    with open(report_path, "w") as f:
        f.write("====================================================\n")
        f.write("QUANTITATIVE STATISTICAL ANALYSIS REPORT\n")
        f.write("====================================================\n\n")
        
        f.write("1. MULTIPLE SEQUENCE ALIGNMENT CONSERVATION\n")
        f.write(f"  Alignment Length: {len(shannon_entropy)} columns\n")
        f.write(f"  Mean Shannon Entropy (H): {np.mean(shannon_entropy):.4f}\n")
        f.write(f"  Median Shannon Entropy (H): {np.median(shannon_entropy):.4f}\n")
        f.write(f"  Perfectly Conserved Columns (H=0): {len(conserved_cols)} ({len(conserved_cols)/len(shannon_entropy)*100:.2f}%)\n\n")
        
        f.write("2. STRUCURAL CONFIDENCE (t-test on pLDDT)\n")
        f.write("  STN7 Stromal Kinase vs Non-Kinase Domains:\n")
        f.write(f"    Kinase Mean pLDDT: {np.mean(stn7_kin):.2f}\n")
        f.write(f"    Non-Kinase Mean pLDDT: {np.mean(stn7_non_kin):.2f}\n")
        f.write(f"    t-statistic: {stn7_t_stat:.4f}\n")
        f.write(f"    p-value: {stn7_p_val:.4e} ({'Highly Significant' if stn7_p_val < 0.01 else 'Significant' if stn7_p_val < 0.05 else 'Not Significant'})\n")
        f.write("  STN8 Stromal Kinase vs Non-Kinase Domains:\n")
        f.write(f"    Kinase Mean pLDDT: {np.mean(stn8_kin):.2f}\n")
        f.write(f"    Non-Kinase Mean pLDDT: {np.mean(stn8_non_kin):.2f}\n")
        f.write(f"    t-statistic: {stn8_t_stat:.4f}\n")
        f.write(f"    p-value: {stn8_p_val:.4e} ({'Highly Significant' if stn8_p_val < 0.01 else 'Significant' if stn8_p_val < 0.05 else 'Not Significant'})\n\n")
        
        f.write("3. HYDROPHOBICITY PROFILE CORRELATION\n")
        f.write("  Pearson Correlation on Interpolated Kyte-Doolittle Profiles:\n")
        f.write(f"    Correlation coefficient r: {r_coef:.4f}\n")
        f.write(f"    p-value: {p_coef:.4e} ({'Significant' if p_coef < 0.05 else 'Not Significant'})\n")
        
    print(f"\nSaved statistical report to: {report_path}")
    print("=========================================================")

if __name__ == "__main__":
    main()
