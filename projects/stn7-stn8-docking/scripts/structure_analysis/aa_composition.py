#!/usr/bin/env python3
"""
STN7 and STN8 Amino Acid Composition Comparative Analysis
Calculates and compares the amino acid composition (mole %)
of Arabidopsis STN7 and STN8 full-length sequences and their
stromal kinase domains. Generates a publication-quality grouped barplot.
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from Bio import SeqIO
from collections import Counter

AMINO_ACIDS_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

# Physicochemical groups for coloring
AA_GROUPS = {
    "Nonpolar": ["A", "V", "I", "L", "M", "F", "W", "P", "G"],
    "Polar": ["S", "T", "C", "Y", "N", "Q"],
    "Positive": ["K", "R", "H"],
    "Negative": ["D", "E"]
}

GROUP_COLORS = {
    "Nonpolar": "#475569",
    "Polar": "#3b82f6",
    "Positive": "#f43f5e",
    "Negative": "#f59e0b"
}

def get_aa_group(aa):
    for group, aas in AA_GROUPS.items():
        if aa in aas:
            return group
    return "Other"

def compute_composition(sequence):
    """
    Returns amino acid composition as percentage (mole %).
    """
    clean_seq = sequence.upper().replace("-", "")
    counts = Counter(clean_seq)
    total = len(clean_seq)
    composition = {aa: (counts.get(aa, 0) / total) * 100 for aa in AMINO_ACIDS_ORDER}
    return composition, total

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    fasta_path = os.path.join(base_dir, "data/fasta/all_homologs.fasta")
    figures_dir = os.path.join(base_dir, "figures")
    results_dir = os.path.join(base_dir, "results")

    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    print("=== Amino Acid Composition Comparative Analysis ===\n")

    records = list(SeqIO.parse(fasta_path, "fasta"))
    stn7_rec = next(r for r in records if "STN7_Arabidopsis" in r.id)
    stn8_rec = next(r for r in records if "STN8_Arabidopsis" in r.id)

    stn7_seq = str(stn7_rec.seq).upper().replace("-", "")
    stn8_seq = str(stn8_rec.seq).upper().replace("-", "")

    # Full-length compositions
    stn7_comp, stn7_len = compute_composition(stn7_seq)
    stn8_comp, stn8_len = compute_composition(stn8_seq)

    # Kinase domain compositions
    stn7_kin_seq = stn7_seq[133:452]   # 0-indexed (residues 134-452)
    stn8_kin_seq = stn8_seq[132:477]   # 0-indexed (residues 133-477)

    stn7_kin_comp, _ = compute_composition(stn7_kin_seq)
    stn8_kin_comp, _ = compute_composition(stn8_kin_seq)

    # Build comparison DataFrame
    df = pd.DataFrame({
        "STN7 (Full)": stn7_comp,
        "STN8 (Full)": stn8_comp,
        "STN7 (Kinase)": stn7_kin_comp,
        "STN8 (Kinase)": stn8_kin_comp
    }, index=AMINO_ACIDS_ORDER)

    csv_path = os.path.join(results_dir, "aa_composition_comparison.csv")
    df.to_csv(csv_path)
    print(f"Saved amino acid composition table to: {csv_path}")
    print(df.to_string())

    # ---- Plot ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
    fig.suptitle("Amino Acid Composition: STN7 vs STN8 (Arabidopsis thaliana)",
                 fontsize=13, fontweight='bold', y=1.02)

    x = np.arange(len(AMINO_ACIDS_ORDER))
    width = 0.4

    # Assign bar colors by physicochemical group
    bar_colors = [GROUP_COLORS[get_aa_group(aa)] for aa in AMINO_ACIDS_ORDER]

    for ax_idx, (label1, label2, title) in enumerate([
        ("STN7 (Full)", "STN8 (Full)", "Full-length Sequences"),
        ("STN7 (Kinase)", "STN8 (Kinase)", "Stromal Kinase Domains Only")
    ]):
        ax = axes[ax_idx]
        v1 = [df.loc[aa, label1] for aa in AMINO_ACIDS_ORDER]
        v2 = [df.loc[aa, label2] for aa in AMINO_ACIDS_ORDER]

        bars1 = ax.bar(x - width/2, v1, width, label="STN7", alpha=0.85,
                       color=[c + "cc" for c in bar_colors],
                       edgecolor="white", linewidth=0.4)
        bars2 = ax.bar(x + width/2, v2, width, label="STN8", alpha=0.65,
                       color=bar_colors,
                       edgecolor="white", linewidth=0.4)

        ax.set_xticks(x)
        ax.set_xticklabels(AMINO_ACIDS_ORDER, fontsize=9)
        ax.set_xlabel("Amino Acid (Single-letter Code)", fontsize=10)
        ax.set_ylabel("Mole %", fontsize=10)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, axis='y', linestyle=':', alpha=0.5)
        ax.set_ylim(0, None)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add physicochemical group annotations
        group_positions = {"Nonpolar": [], "Polar": [], "Positive": [], "Negative": []}
        for i, aa in enumerate(AMINO_ACIDS_ORDER):
            group = get_aa_group(aa)
            group_positions[group].append(i)

        for group, positions in group_positions.items():
            if positions:
                mid = np.mean(positions)
                ax.text(mid, ax.get_ylim()[1] * 0.95, group,
                        ha='center', va='top', fontsize=7,
                        color=GROUP_COLORS[group], fontstyle='italic')

    # Legend for physicochemical groups
    group_patches = [plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.8, label=g)
                     for g, c in GROUP_COLORS.items()]
    fig.legend(handles=group_patches, title="AA Property", loc='lower center',
               ncol=4, fontsize=9, framealpha=0.9, bbox_to_anchor=(0.5, -0.07))

    plt.tight_layout()
    fig_path = os.path.join(figures_dir, "aa_composition_comparison.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"\nSaved amino acid composition figure to: {fig_path}")
    print("===================================================")

if __name__ == "__main__":
    main()
