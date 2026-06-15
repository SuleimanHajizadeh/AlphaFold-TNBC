#!/usr/bin/env python3
"""
AlphaFold2 pLDDT Analysis — AKT1 (TNBC Hub Kinase)
Author : Suleiman Hajizadeh | IMBB, Azerbaijan
Input  : AKT1_TNBC_42642_0/ (ColabFold v1.6.1 output)
Output : figures/AKT1_pLDDT.png

Description:
    Extracts per-residue pLDDT (Predicted Local Distance Difference Test)
    confidence scores from the AlphaFold2 top-ranked model JSON file and
    renders a publication-ready plot with colour-coded confidence zones
    following DeepMind/EBI standard thresholds.

    pLDDT is the primary quality metric for AlphaFold2 predictions:
      ≥ 90  Very high confidence  — accurate to atomic resolution
      70-90 Confident            — accurate backbone, loops may vary
      50-70 Low confidence       — correct fold likely, details uncertain
      < 50  Very low             — intrinsically disordered region (IDR)
"""

import os
import sys
import json
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

os.makedirs("figures", exist_ok=True)

# ── Locate best-ranked JSON ────────────────────────────────────────────────
ALPHAFOLD_DIR = "AKT1_TNBC_42642_0"
PROTEIN       = "AKT1 (ColabFold v1.6.1 | AlphaFold2-PTM)"
OUTPUT_FIG    = "figures/AKT1_pLDDT.png"

def find_rank001_json(directory: str) -> str:
    pattern = os.path.join(directory, "*scores_rank_001*.json")
    files   = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"No rank_001 scores JSON found in {directory}")
    return files[0]

# ── Load pLDDT from JSON ───────────────────────────────────────────────────
def load_plddt(json_path: str):
    with open(json_path) as f:
        data = json.load(f)
    plddt = np.array(data["plddt"])
    ptm   = data.get("ptm", None)
    return plddt, ptm

# ── Colour by confidence band ──────────────────────────────────────────────
BANDS = [
    (90, 100, "#0053D6", "Very high (≥90)"),
    (70,  90, "#65CBF3", "Confident (70–90)"),
    (50,  70, "#FFDB13", "Low (50–70)"),
    ( 0,  50, "#FF7D45", "Very low (<50) — IDR"),
]

def colour_by_band(plddt):
    colours = []
    for v in plddt:
        for lo, hi, col, _ in BANDS:
            if lo <= v < hi or (hi == 100 and v >= lo):
                colours.append(col)
                break
    return colours

# ── Summary statistics ─────────────────────────────────────────────────────
def print_summary(plddt, ptm, json_path):
    print(f"\n{'='*50}")
    print(f"  AlphaFold2 pLDDT Analysis — AKT1")
    print(f"  Model: {os.path.basename(json_path)}")
    print(f"{'='*50}")
    print(f"  Total residues       : {len(plddt)}")
    print(f"  Mean pLDDT           : {plddt.mean():.2f}")
    print(f"  Median pLDDT         : {np.median(plddt):.2f}")
    print(f"  Min / Max            : {plddt.min():.2f} / {plddt.max():.2f}")
    print(f"  pTM score            : {ptm:.3f}" if ptm else "  pTM score           : N/A")
    print(f"  Very high (≥90)      : {(plddt>=90).sum():3d}  ({(plddt>=90).mean()*100:.1f}%)")
    print(f"  Confident  (≥70)     : {(plddt>=70).sum():3d}  ({(plddt>=70).mean()*100:.1f}%)")
    print(f"  Low        (50–70)   : {((plddt>=50)&(plddt<70)).sum():3d}  ({((plddt>=50)&(plddt<70)).mean()*100:.1f}%)")
    print(f"  IDR        (<50)     : {(plddt<50).sum():3d}  ({(plddt<50).mean()*100:.1f}%)")
    print(f"{'='*50}\n")

# ── Plot ───────────────────────────────────────────────────────────────────
def plot_plddt(plddt, output_path):
    residues = np.arange(1, len(plddt) + 1)
    colours  = colour_by_band(plddt)

    fig, ax = plt.subplots(figsize=(14, 5))

    ax.bar(residues, plddt, color=colours, width=1.0, linewidth=0, zorder=2)

    # Threshold reference lines
    for thr, _, col, _ in BANDS[:-1]:
        ax.axhline(thr, color=col, linewidth=0.8, linestyle="--", alpha=0.7)

    # Shaded background bands
    shade = [
        (90, 100, "#EEF4FF"), (70, 90, "#F0FAFF"),
        (50,  70, "#FFFDE7"), ( 0, 50, "#FFF3E0"),
    ]
    for lo, hi, bg in shade:
        ax.axhspan(lo, hi, alpha=0.25, color=bg, zorder=0)

    # Band labels on right axis
    for lo, hi, col, lbl in BANDS:
        mid = (lo + hi) / 2
        ax.text(len(plddt) * 1.005, mid, lbl,
                color=col, fontsize=7.5, va="center", ha="left")

    # Legend
    patches = [mpatches.Patch(color=c, label=l) for _, _, c, l in BANDS]
    ax.legend(handles=patches, loc="lower left", fontsize=9, framealpha=0.9)

    ax.set_xlim(0, len(plddt) + len(plddt) * 0.12)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Residue Position", fontsize=12)
    ax.set_ylabel("pLDDT Score", fontsize=12)
    ax.set_title(
        f"AlphaFold2 — Per-residue pLDDT Confidence | {PROTEIN}\n"
        f"Mean pLDDT: {plddt.mean():.1f} | pTM: N/A",
        fontsize=12,
    )
    ax.tick_params(labelsize=10)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"[✓] pLDDT figure saved → {output_path}")

# ── Main ───────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if not os.path.isdir(ALPHAFOLD_DIR):
        print(f"[!] AlphaFold output directory not found: {ALPHAFOLD_DIR}")
        sys.exit(1)

    json_path    = find_rank001_json(ALPHAFOLD_DIR)
    plddt, ptm   = load_plddt(json_path)
    print_summary(plddt, ptm, json_path)
    plot_plddt(plddt, OUTPUT_FIG)
