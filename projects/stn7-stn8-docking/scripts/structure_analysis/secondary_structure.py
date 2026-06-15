#!/usr/bin/env python3
"""
STN7 vs STN8 Secondary Structure Analysis and Comparison
Calculates the secondary structure composition (Helix, Sheet, Coil)
of full-length sequences and kinase domains of STN7 and STN8
using dihedral angles (phi/psi) from PDB coordinates.
Generates a publication-quality stacked bar plot.
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB

def classify_secondary_structure(phi, psi):
    """
    Classifies secondary structure using Ramachandran dihedral angles.
    """
    if phi is None or psi is None:
        return 'Coil'
    
    phi_deg = math.degrees(phi)
    psi_deg = math.degrees(psi)
    
    # Helix region: phi [-120, -30], psi [-90, 0]
    if -120 <= phi_deg <= -30 and -90 <= psi_deg <= 0:
        return 'Helix'
    # Sheet region: phi [-180, -45], psi [90, 180] or [-180, -150]
    elif -180 <= phi_deg <= -45 and (90 <= psi_deg <= 180 or -180 <= psi_deg <= -150):
        return 'Sheet'
    else:
        return 'Coil'

def analyze_pdb_secondary_structure(pdb_path, domain_range=None):
    """
    Parses a PDB file and calculates secondary structure composition.
    Returns counts of Helix, Sheet, Coil, and total residues.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = structure[0]
    
    ppbuilder = PDB.PPBuilder()
    peptides = ppbuilder.build_peptides(model['A'])
    
    counts = {'Helix': 0, 'Sheet': 0, 'Coil': 0}
    total = 0
    
    for poly in peptides:
        phi_psi = poly.get_phi_psi_list()
        for i, res in enumerate(poly):
            res_num = res.get_id()[1]
            if domain_range is not None and res_num not in domain_range:
                continue
            phi, psi = phi_psi[i]
            ss = classify_secondary_structure(phi, psi)
            counts[ss] += 1
            total += 1
            
    return counts, total

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    pdb_7 = os.path.join(base_dir, "data/structures/Q9S713_AlphaFold.pdb")
    pdb_8 = os.path.join(base_dir, "data/structures/Q9LZV4_AlphaFold.pdb")
    
    figures_dir = os.path.join(base_dir, "figures")
    results_dir = os.path.join(base_dir, "results")
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    print("=== Starting Secondary Structure Analysis ===")
    
    # Define targets and domains
    # Kinase domains: STN7 (134-452), STN8 (133-477)
    targets = [
        {"name": "STN7 Full-Length", "pdb": pdb_7, "range": None},
        {"name": "STN7 Kinase Domain", "pdb": pdb_7, "range": range(134, 453)},
        {"name": "STN8 Full-Length", "pdb": pdb_8, "range": None},
        {"name": "STN8 Kinase Domain", "pdb": pdb_8, "range": range(133, 478)}
    ]
    
    results = {}
    
    # Save text report
    report_path = os.path.join(results_dir, "secondary_structure_comparison.txt")
    with open(report_path, "w") as f:
        f.write("=================================================================\n")
        f.write("STN7 vs STN8 Secondary Structure Composition Report\n")
        f.write("Calculated using Ramachandran dihedral angles (phi/psi)\n")
        f.write("=================================================================\n\n")
        
        for t in targets:
            counts, total = analyze_pdb_secondary_structure(t["pdb"], t["range"])
            results[t["name"]] = {
                "counts": counts,
                "total": total,
                "pct": {ss: (count / total * 100) if total > 0 else 0 for ss, count in counts.items()}
            }
            
            f.write(f"--- {t['name']} ---\n")
            f.write(f"  Total residues analyzed: {total}\n")
            for ss in ['Helix', 'Sheet', 'Coil']:
                count = counts[ss]
                pct = results[t["name"]]["pct"][ss]
                f.write(f"  {ss:<8}: {count:<4} ({pct:.2f}%)\n")
            f.write("\n")
            
            print(f"{t['name']}: Helix={pct_h:.1f}%, Sheet={pct_e:.1f}%, Coil={pct_c:.1f}%" 
                  if False else f"{t['name']}: Total={total}, H={counts['Helix']}, E={counts['Sheet']}, C={counts['Coil']}")
            
    print(f"Saved secondary structure text report to: {report_path}")
    
    # ---- Draw publication-quality stacked bar plot ----
    labels = list(results.keys())
    helix_pcts = [results[lbl]["pct"]["Helix"] for lbl in labels]
    sheet_pcts = [results[lbl]["pct"]["Sheet"] for lbl in labels]
    coil_pcts = [results[lbl]["pct"]["Coil"] for lbl in labels]
    
    # Use modern, rich colors
    color_helix = '#0ea5e9'  # Modern Sky Blue
    color_sheet = '#f43f5e'  # Modern Rose/Red
    color_coil = '#cbd5e1'   # Slate Gray
    
    fig, ax = plt.subplots(figsize=(10, 5))
    
    # Draw horizontal stacked bars
    y_pos = np.arange(len(labels))
    bar_width = 0.55
    
    # Plot Helix
    ax.barh(y_pos, helix_pcts, bar_width, label='Alpha-Helix', color=color_helix, edgecolor='white', linewidth=0.5)
    # Plot Sheet (stacked on Helix)
    ax.barh(y_pos, sheet_pcts, bar_width, left=helix_pcts, label='Beta-Sheet', color=color_sheet, edgecolor='white', linewidth=0.5)
    # Plot Coil (stacked on Helix + Sheet)
    left_coil = np.array(helix_pcts) + np.array(sheet_pcts)
    ax.barh(y_pos, coil_pcts, bar_width, left=left_coil, label='Random Coil / Loop', color=color_coil, edgecolor='white', linewidth=0.5)
    
    # Add values on top of bars
    for i in range(len(labels)):
        # Helix label
        h_val = helix_pcts[i]
        if h_val > 5:
            ax.text(h_val / 2, i, f"{h_val:.1f}%", ha='center', va='center', color='white', fontweight='bold', fontsize=9)
        # Sheet label
        s_val = sheet_pcts[i]
        if s_val > 5:
            ax.text(h_val + s_val / 2, i, f"{s_val:.1f}%", ha='center', va='center', color='white', fontweight='bold', fontsize=9)
        # Coil label
        c_val = coil_pcts[i]
        if c_val > 5:
            ax.text(h_val + s_val + c_val / 2, i, f"{c_val:.1f}%", ha='center', va='center', color='#334155', fontweight='bold', fontsize=9)
            
    # Styling
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10, fontweight='bold')
    ax.set_xlabel('Percentage of Resides (%)', fontsize=11, labelpad=8)
    ax.set_xlim(0, 100)
    ax.set_title('Secondary Structure Composition Comparison\nSTN7 vs STN8 (Full-Length & Kinase Domains)', fontsize=12, fontweight='bold', pad=15)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#cbd5e1')
    ax.spines['bottom'].set_color('#cbd5e1')
    
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.22), ncol=3, framealpha=0.9, fontsize=9.5)
    plt.tight_layout()
    
    plot_path = os.path.join(figures_dir, "secondary_structure_comparison.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved secondary structure comparison plot to: {plot_path}")
    print("=============================================")

if __name__ == "__main__":
    main()
