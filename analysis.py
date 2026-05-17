#!/usr/bin/env python3
"""
AKT1 Structural Analysis — B-factor Profile & Secondary Structure
Author : Suleiman Hajizadeh | IMBB, Azerbaijan
Dataset: AKT1 (PDB: 4EJN) — Crystal structure of autoinhibited AKT1
         in complex with kinase inhibitor. Resolution: 2.20 Å
Description:
    Parses the experimental AKT1 PDB file, extracts per-residue B-factors
    (crystallographic temperature factors — proxy for structural flexibility
    and intrinsic disorder), and generates a publication-ready figure.
    High B-factor regions correlate with structural flexibility or disorder,
    analogous to low pLDDT in AlphaFold2 predictions.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict
from Bio import PDB
from Bio.PDB import DSSP

# ──────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────
PDB_FILE   = "data/AKT1_ranked_0.pdb"
PDB_ID     = "4EJN"
PROTEIN    = f"AKT1 (PDB: {PDB_ID}) — Homo sapiens"
CHAIN_ID   = "A"                          # Primary AKT1 chain
OUTPUT_DIR = "figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ──────────────────────────────────────────────
# 1. Parse B-factors from Cα atoms
# ──────────────────────────────────────────────
def extract_bfactors(pdb_path: str, chain: str = "A"):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("AKT1", pdb_path)

    res_ids, bfactors, res_names = [], [], []
    for model in structure:
        if chain not in [c.id for c in model]:
            chain = list(model.get_chains())[0].id  # fallback to first chain
        for residue in model[chain]:
            if "CA" in residue and residue.get_id()[0] == " ":  # ATOM records only
                res_ids.append(residue.get_id()[1])
                bfactors.append(residue["CA"].get_bfactor())
                res_names.append(residue.get_resname())
    return np.array(res_ids), np.array(bfactors), res_names

# ──────────────────────────────────────────────
# 2. Classify flexibility zones
# ──────────────────────────────────────────────
def classify_bfactor(bfactors):
    """
    Classify residues by B-factor flexibility:
      Low  (<20)  → rigid, well-ordered → dark blue
      Med  (20-40) → moderate motion    → light blue
      High (40-60) → flexible loop      → yellow
      VHigh(>60)   → disordered region  → orange
    """
    colours = []
    for v in bfactors:
        if v < 20:
            colours.append("#0053D6")
        elif v < 40:
            colours.append("#65CBF3")
        elif v < 60:
            colours.append("#FFDB13")
        else:
            colours.append("#FF7D45")
    return colours

# ──────────────────────────────────────────────
# 3. Summary statistics
# ──────────────────────────────────────────────
def print_summary(res_ids, bfactors):
    print(f"\n{'='*45}")
    print(f"  AKT1 B-factor Structural Analysis — {PDB_ID}")
    print(f"{'='*45}")
    print(f"  Residues analysed  : {len(bfactors)}")
    print(f"  Mean B-factor      : {bfactors.mean():.2f} Å²")
    print(f"  Median B-factor    : {np.median(bfactors):.2f} Å²")
    print(f"  Min / Max          : {bfactors.min():.2f} / {bfactors.max():.2f} Å²")
    print(f"  Rigid (<20 Å²)     : {(bfactors < 20).sum()} res  ({(bfactors < 20).mean()*100:.1f}%)")
    print(f"  Flexible (>40 Å²)  : {(bfactors > 40).sum()} res  ({(bfactors > 40).mean()*100:.1f}%)")
    print(f"{'='*45}\n")

# ──────────────────────────────────────────────
# 4. Publication-ready plot
# ──────────────────────────────────────────────
def plot_bfactor(res_ids, bfactors, out_path):
    colours = classify_bfactor(bfactors)

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.bar(res_ids, bfactors, color=colours, width=1.0, linewidth=0, zorder=2)

    # Shaded domain annotations (AKT1 domain boundaries)
    domains = [
        (1,   107, "#E8F4FD", "PH Domain\n(lipid binding)"),
        (108, 152, "#FFF9E6", "Linker"),
        (153, 408, "#E8FDE8", "Kinase Domain\n(catalytic)"),
        (409, 480, "#FDE8F4", "Regulatory\nDomain (HM)"),
    ]
    for start, end, colour, label in domains:
        ax.axvspan(start, end, alpha=0.25, color=colour, zorder=1)
        mid = (start + end) / 2
        ax.text(mid, ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] > 0 else 60,
                label, ha="center", va="top", fontsize=7.5, color="#444",
                style="italic")

    # Reference lines
    for thr, col in [(20, "#0053D6"), (40, "#FFDB13"), (60, "#FF7D45")]:
        ax.axhline(thr, color=col, linewidth=0.8, linestyle="--", alpha=0.6)

    # Legend
    legend_labels = {
        "Rigid (<20 Å²)":      "#0053D6",
        "Moderate (20–40 Å²)": "#65CBF3",
        "Flexible (40–60 Å²)": "#FFDB13",
        "Disordered (>60 Å²)": "#FF7D45",
    }
    patches = [mpatches.Patch(color=c, label=l) for l, c in legend_labels.items()]
    ax.legend(handles=patches, loc="upper right", fontsize=9, framealpha=0.9)

    ax.set_xlim(res_ids.min() - 2, res_ids.max() + 2)
    ax.set_xlabel("Residue Position", fontsize=12)
    ax.set_ylabel("B-factor (Å²)", fontsize=12)
    ax.set_title(
        f"Crystallographic B-factor Profile | {PROTEIN}\n"
        f"Structural flexibility mapped across functional domains",
        fontsize=12,
    )
    ax.tick_params(labelsize=10)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"[✓] Figure saved → {out_path}")

# ──────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────
if __name__ == "__main__":
    if not os.path.exists(PDB_FILE):
        print(f"[!] PDB file not found: {PDB_FILE}")
        sys.exit(1)

    res_ids, bfactors, res_names = extract_bfactors(PDB_FILE, chain=CHAIN_ID)
    print_summary(res_ids, bfactors)
    plot_bfactor(res_ids, bfactors, f"{OUTPUT_DIR}/AKT1_Bfactor.png")
