#!/usr/bin/env python3
"""
AKT1 Ramachandran Plot
Author : Suleiman Hajizadeh | IMBB, Azerbaijan
Input  : data/AKT1_ranked_0.pdb  (PDB 4EJN — 2.20 Å crystal structure)
Output : figures/AKT1_Ramachandran.png

Description:
    Extracts backbone torsion angles (φ phi, ψ psi) from Cα coordinates
    using Biopython and plots a Ramachandran diagram to characterise the
    secondary structure composition of AKT1. Favoured regions (α-helix,
    β-sheet) are highlighted using standard Ramachandran boundaries.
"""

import os
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder

warnings.filterwarnings("ignore")

# ── Config ─────────────────────────────────────────────────────────────────
PDB_FILE   = "data/AKT1_ranked_0.pdb"
OUTPUT_FIG = "figures/AKT1_Ramachandran.png"
CHAIN_ID   = "A"

os.makedirs("figures", exist_ok=True)

# ── 1. Extract φ/ψ angles ──────────────────────────────────────────────────
def get_phi_psi(pdb_path: str):
    parser  = PDB.PDBParser(QUIET=True)
    struct  = parser.get_structure("AKT1", pdb_path)
    builder = PPBuilder()
    phi_list, psi_list, aa_list = [], [], []

    for pp in builder.build_peptides(struct):
        angles = pp.get_phi_psi_list()
        for res, (phi, psi) in zip(pp, angles):
            if phi is not None and psi is not None:
                phi_list.append(np.degrees(phi))
                psi_list.append(np.degrees(psi))
                aa_list.append(res.get_resname())

    return np.array(phi_list), np.array(psi_list), aa_list

# ── 2. Colour by secondary structure region ────────────────────────────────
def classify_region(phi, psi):
    """Returns colour code for each (φ, ψ) point."""
    colours = []
    for ph, ps in zip(phi, psi):
        # α-helix favoured: φ ≈ -60±30, ψ ≈ -45±30
        if -90 <= ph <= -30 and -70 <= ps <= -10:
            colours.append("#2166AC")   # blue
        # β-strand favoured: φ ≈ -120±40, ψ ≈ +130±40
        elif -160 <= ph <= -80 and 80 <= ps <= 180:
            colours.append("#1A9641")   # green
        # Left-handed helix (Gly): φ>0, ψ>0
        elif 0 <= ph <= 90 and 0 <= ps <= 90:
            colours.append("#F46D43")   # orange
        # All other allowed
        else:
            colours.append("#B2ABD2")   # grey-purple
    return colours

# ── 3. Plot ────────────────────────────────────────────────────────────────
def plot_ramachandran(phi, psi, colours, out_path):
    fig, ax = plt.subplots(figsize=(7, 7))

    # Shaded favoured regions (approximate)
    ax.axvspan(-90, -30, ymin=0, ymax=1, alpha=0.04, color="#2166AC")
    ax.add_patch(plt.Rectangle((-90, -70), 60, 60, linewidth=0,
                                edgecolor="none", facecolor="#2166AC", alpha=0.08,
                                label="__nolegend__"))
    ax.add_patch(plt.Rectangle((-160, 80), 80, 100, linewidth=0,
                                edgecolor="none", facecolor="#1A9641", alpha=0.08,
                                label="__nolegend__"))

    # Reference axes
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")

    # Data points
    ax.scatter(phi, psi, c=colours, s=16, alpha=0.75, linewidths=0, zorder=3)

    # Legend
    patches = [
        mpatches.Patch(color="#2166AC", label=f"α-helix region  (n={sum(1 for c in colours if c=='#2166AC')})"),
        mpatches.Patch(color="#1A9641", label=f"β-strand region (n={sum(1 for c in colours if c=='#1A9641')})"),
        mpatches.Patch(color="#F46D43", label=f"Left-handed Gly (n={sum(1 for c in colours if c=='#F46D43')})"),
        mpatches.Patch(color="#B2ABD2", label=f"Other allowed   (n={sum(1 for c in colours if c=='#B2ABD2')})"),
    ]
    ax.legend(handles=patches, loc="lower right", fontsize=9, framealpha=0.9)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel("φ (phi) angle  [degrees]", fontsize=12)
    ax.set_ylabel("ψ (psi) angle  [degrees]", fontsize=12)
    ax.set_title("Ramachandran Plot — AKT1 (PDB: 4EJN, 2.20 Å)\n"
                 "Backbone torsion angles coloured by secondary-structure region",
                 fontsize=11)
    ax.tick_params(labelsize=10)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"[✓] Ramachandran plot saved → {out_path}")

# ── 4. Summary ─────────────────────────────────────────────────────────────
def print_summary(phi, psi):
    total = len(phi)
    alpha = sum(1 for ph, ps in zip(phi, psi) if -90<=ph<=-30 and -70<=ps<=-10)
    beta  = sum(1 for ph, ps in zip(phi, psi) if -160<=ph<=-80 and 80<=ps<=180)
    left  = sum(1 for ph, ps in zip(phi, psi) if 0<=ph<=90 and 0<=ps<=90)
    other = total - alpha - beta - left
    print(f"\n{'='*46}")
    print(f"  AKT1 Ramachandran Analysis — PDB 4EJN")
    print(f"{'='*46}")
    print(f"  Total residues with φ/ψ : {total}")
    print(f"  α-helix favoured        : {alpha:3d}  ({alpha/total*100:.1f}%)")
    print(f"  β-strand favoured       : {beta:3d}  ({beta/total*100:.1f}%)")
    print(f"  Left-handed (Gly)       : {left:3d}  ({left/total*100:.1f}%)")
    print(f"  Other allowed           : {other:3d}  ({other/total*100:.1f}%)")
    print(f"{'='*46}\n")

# ── Main ───────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if not os.path.exists(PDB_FILE):
        print(f"[!] PDB file not found: {PDB_FILE}")
        sys.exit(1)
    phi, psi, aa = get_phi_psi(PDB_FILE)
    print_summary(phi, psi)
    colours = classify_region(phi, psi)
    plot_ramachandran(phi, psi, colours, OUTPUT_FIG)
