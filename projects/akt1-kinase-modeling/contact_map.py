#!/usr/bin/env python3
"""
AKT1 Cα Contact Map
Author : Suleiman Hajizadeh | IMBB, Azerbaijan
Input  : AKT1_TNBC_42642_0/*rank_001*.pdb  (ColabFold v1.6.1 / AlphaFold2)
Output : figures/AKT1_ContactMap.png

Description:
    Computes the pairwise Euclidean distance matrix between all Cα atoms
    in the top-ranked AlphaFold2 predicted structure and renders a contact
    map — residues within an 8 Å threshold are marked as "in contact".
    Contact maps reveal secondary structure topology and long-range tertiary
    interactions within the predicted AKT1 structure.
"""

import os
import sys
import glob
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import PDB

warnings.filterwarnings("ignore")

# ── Config ─────────────────────────────────────────────────────────────────
_pdb_hits   = glob.glob("AKT1_TNBC_42642_0/*rank_001*.pdb")
PDB_FILE    = _pdb_hits[0] if _pdb_hits else "data/AKT1_ranked_0.pdb"
OUTPUT_FIG  = "figures/AKT1_ContactMap.png"
CHAIN_ID    = "A"
THRESHOLD   = 8.0   # Å — standard contact-map cutoff

os.makedirs("figures", exist_ok=True)

# ── 1. Extract Cα coordinates ───────────────────────────────────────────────
def get_ca_coords(pdb_path: str, chain: str = "A"):
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("AKT1", pdb_path)
    coords, res_ids = [], []

    for model in struct:
        # Fallback if requested chain absent
        available = [c.id for c in model]
        if chain not in available:
            chain = available[0]

        for residue in model[chain]:
            if "CA" in residue and residue.get_id()[0] == " ":
                coords.append(residue["CA"].get_vector().get_array())
                res_ids.append(residue.get_id()[1])

    return np.array(coords), np.array(res_ids)

# ── 2. Distance matrix & contact map ───────────────────────────────────────
def compute_contacts(coords: np.ndarray, threshold: float):
    n = len(coords)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            dist[i, j] = np.linalg.norm(coords[i] - coords[j])
    contact = (dist < threshold).astype(int)
    # Exclude self-contacts and trivially bonded neighbours (|i-j| <= 2)
    for i in range(n):
        for j in range(max(0, i-2), min(n, i+3)):
            contact[i, j] = 0
    return dist, contact

# ── 3. Summary ──────────────────────────────────────────────────────────────
def print_summary(dist, contact, res_ids):
    n = len(res_ids)
    total_pairs   = n * (n - 1) // 2
    contact_pairs = int(contact.sum() / 2)
    print(f"\n{'='*46}")
    print(f"  AKT1 Contact Map Analysis — PDB 4EJN")
    print(f"{'='*46}")
    print(f"  Residues (Cα)         : {n}")
    print(f"  Residue range         : {res_ids[0]}–{res_ids[-1]}")
    print(f"  Contact threshold     : {THRESHOLD} Å")
    print(f"  Total residue pairs   : {total_pairs}")
    print(f"  Contacts detected     : {contact_pairs}  ({contact_pairs/total_pairs*100:.1f}%)")
    print(f"  Min Cα–Cα distance    : {dist[dist>0].min():.2f} Å")
    print(f"  Max Cα–Cα distance    : {dist.max():.2f} Å")
    print(f"{'='*46}\n")

# ── 4. Plot ─────────────────────────────────────────────────────────────────
def plot_contact_map(contact, res_ids, out_path):
    fig, ax = plt.subplots(figsize=(9, 8))

    ax.imshow(contact, cmap="Blues", origin="lower",
              extent=[res_ids[0], res_ids[-1],
                      res_ids[0], res_ids[-1]],
              aspect="auto", interpolation="none", vmin=0, vmax=1)

    # Domain boundary lines
    boundaries = [
        (107, "PH | Linker",   "#E74C3C"),
        (152, "Linker | Kinase","#E67E22"),
        (408, "Kinase | HM",   "#8E44AD"),
    ]
    for pos, label, col in boundaries:
        ax.axhline(pos, color=col, linewidth=0.9, linestyle="--", alpha=0.7)
        ax.axvline(pos, color=col, linewidth=0.9, linestyle="--", alpha=0.7)
        ax.text(pos + 3, res_ids[-1] - 10, label,
                color=col, fontsize=7, va="top", rotation=90, alpha=0.85)

    # Domain region labels
    regions = [
        (55,  55,  "PH\nDomain",   "#2980B9"),
        (280, 280, "Kinase\nDomain","#27AE60"),
        (440, 440, "HM\nDomain",   "#8E44AD"),
    ]
    for rx, ry, lbl, col in regions:
        if rx <= res_ids[-1] and ry <= res_ids[-1]:
            ax.text(rx, ry, lbl, color=col, fontsize=7.5,
                    ha="center", va="center", style="italic", alpha=0.6)

    ax.set_xlabel("Residue Index", fontsize=11)
    ax.set_ylabel("Residue Index", fontsize=11)
    ax.set_title(
        f"Cα Contact Map — AKT1 (PDB: 4EJN, 2.20 Å)\n"
        f"Threshold: {THRESHOLD} Å | Blue = contact; White = no contact",
        fontsize=11,
    )
    ax.tick_params(labelsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"[✓] Contact map saved → {out_path}")

# ── Main ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if not os.path.exists(PDB_FILE):
        print(f"[!] PDB file not found: {PDB_FILE}")
        sys.exit(1)

    coords, res_ids = get_ca_coords(PDB_FILE, CHAIN_ID)
    dist, contact   = compute_contacts(coords, THRESHOLD)
    print_summary(dist, contact, res_ids)
    plot_contact_map(contact, res_ids, OUTPUT_FIG)
