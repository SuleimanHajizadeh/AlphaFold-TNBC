#!/usr/bin/env python3
"""
AlphaFold2 pLDDT Score Analysis for AKT1 (Human)
Author: Suleiman Hajizadeh | IMBB, Azerbaijan
Description:
    Parses AlphaFold2-predicted PDB file, extracts per-residue pLDDT
    confidence scores from B-factor column, and generates a publication-
    ready plot with colour-coded confidence zones.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import PDB

# ─────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────
PDB_FILE   = "data/AKT1_ranked_0.pdb"   # ColabFold output name
OUTPUT_FIG = "figures/AKT1_pLDDT.png"
PROTEIN    = "AKT1 (P31749) — Homo sapiens"

# pLDDT confidence thresholds (DeepMind/EBI standard)
THRESHOLDS = {
    "Very high (≥90)":  (90, 100, "#0053D6"),   # blue
    "Confident (70–90)": (70,  90, "#65CBF3"),   # light blue
    "Low (50–70)":       (50,  70, "#FFDB13"),   # yellow
    "Very low (<50)":    ( 0,  50, "#FF7D45"),   # orange
}

# ─────────────────────────────────────────────
# Parse pLDDT from PDB B-factor column
# ─────────────────────────────────────────────
def extract_plddt(pdb_path: str):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    residue_ids = []
    plddt_scores = []

    for model in structure:
        for chain in model:
            for residue in chain:
                # Take CA atom bfactor as pLDDT (AlphaFold convention)
                if "CA" in residue:
                    ca = residue["CA"]
                    residue_ids.append(residue.get_id()[1])
                    plddt_scores.append(ca.get_bfactor())
    return np.array(residue_ids), np.array(plddt_scores)


# ─────────────────────────────────────────────
# Colour each residue by confidence band
# ─────────────────────────────────────────────
def colour_by_band(plddt):
    colours = []
    for v in plddt:
        if v >= 90:
            colours.append("#0053D6")
        elif v >= 70:
            colours.append("#65CBF3")
        elif v >= 50:
            colours.append("#FFDB13")
        else:
            colours.append("#FF7D45")
    return colours


# ─────────────────────────────────────────────
# Plot
# ─────────────────────────────────────────────
def plot_plddt(res_ids, plddt, output_path):
    fig, ax = plt.subplots(figsize=(14, 5))

    colours = colour_by_band(plddt)
    ax.bar(res_ids, plddt, color=colours, width=1.0, linewidth=0)

    # Threshold reference lines
    for thr, colour in zip([90, 70, 50], ["#0053D6", "#65CBF3", "#FFDB13"]):
        ax.axhline(thr, color=colour, linewidth=0.8, linestyle="--", alpha=0.6)

    # Legend
    patches = [
        mpatches.Patch(color=c, label=l)
        for l, (lo, hi, c) in THRESHOLDS.items()
    ]
    ax.legend(handles=patches, loc="lower right", fontsize=9, framealpha=0.9)

    ax.set_xlim(res_ids.min() - 1, res_ids.max() + 1)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Residue Position", fontsize=12)
    ax.set_ylabel("pLDDT Score", fontsize=12)
    ax.set_title(f"AlphaFold2 — Per-residue pLDDT Confidence | {PROTEIN}", fontsize=13)
    ax.tick_params(labelsize=10)
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)
    print(f"[✓] Figure saved → {output_path}")


# ─────────────────────────────────────────────
# Summary statistics
# ─────────────────────────────────────────────
def print_summary(plddt):
    print("\n===== pLDDT Summary Statistics =====")
    print(f"  Total residues  : {len(plddt)}")
    print(f"  Mean pLDDT      : {plddt.mean():.2f}")
    print(f"  Median pLDDT    : {np.median(plddt):.2f}")
    print(f"  Min / Max       : {plddt.min():.2f} / {plddt.max():.2f}")
    print(f"  Very high (≥90) : {(plddt >= 90).sum()} residues  ({(plddt >= 90).mean()*100:.1f}%)")
    print(f"  Confident (70+) : {(plddt >= 70).sum()} residues  ({(plddt >= 70).mean()*100:.1f}%)")
    print(f"  Low (<50)       : {(plddt <  50).sum()} residues  ({(plddt <  50).mean()*100:.1f}%)")
    print("=====================================\n")


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
if __name__ == "__main__":
    if not os.path.exists(PDB_FILE):
        print(f"[!] PDB file not found: {PDB_FILE}")
        print("    Please place the ColabFold output PDB (ranked_0.pdb) in data/")
        sys.exit(1)

    res_ids, plddt = extract_plddt(PDB_FILE)
    print_summary(plddt)
    plot_plddt(res_ids, plddt, OUTPUT_FIG)
