#!/usr/bin/env python3
"""
STN7 and STN8 Publication-Quality Multi-Panel Summary Figure
Combines the 6 key results into one integrated figure panel:
  Panel A: Sequence Identity Heatmap
  Panel B: pLDDT Confidence Profile (Arabidopsis)
  Panel C: Hydrophobicity Profile (Arabidopsis)
  Panel D: Shannon Entropy Conservation Profile
  Panel E: Kinase CA Distance Profile
  Panel F: Amino Acid Composition (Kinase domains)
"""

import os
import re
import json
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from Bio import SeqIO, PDB
from collections import Counter

# ─── Constants ───────────────────────────────────────────────
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

AMINO_ACIDS_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
KD_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5,
    'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
    'Y': -1.3, 'V': 4.2
}

D3_TO_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def extract_plddt(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("p", pdb_file)
    res_nums, plddts = [], []
    for residue in structure[0]['A']:
        if PDB.Polypeptide.is_aa(residue):
            try:
                res_nums.append(residue.get_id()[1])
                plddts.append(residue['CA'].get_bfactor())
            except KeyError:
                pass
    return np.array(res_nums), np.array(plddts)

def kd_profile(seq, w=19):
    L = len(seq)
    profile = np.zeros(L)
    hw = w // 2
    for i in range(L):
        window = seq[max(0, i-hw): min(L, i+hw+1)]
        vals = [KD_SCALE[a] for a in window if a in KD_SCALE]
        profile[i] = np.mean(vals) if vals else 0.0
    return profile

def shannon_entropy_per_col(records):
    align_len = len(records[0].seq)
    entropy = []
    for col in range(align_len):
        column = [str(r.seq)[col].upper() for r in records]
        counts = {}
        total = 0
        for c in column:
            if c != '-':
                counts[c] = counts.get(c, 0) + 1
                total += 1
        if total == 0:
            entropy.append(0.0)
            continue
        H = 0.0
        for cnt in counts.values():
            p = cnt / total
            H -= p * math.log2(p)
        entropy.append(H)
    return np.array(entropy)

def aa_composition(seq):
    clean = seq.upper().replace('-', '')
    counts = Counter(clean)
    total = len(clean)
    return {aa: (counts.get(aa, 0) / total) * 100 for aa in AMINO_ACIDS_ORDER}

def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    print("=== Generating Multi-Panel Publication Summary Figure ===\n")

    # ── Load Data ──────────────────────────────────────────────
    records = list(SeqIO.parse(os.path.join(BASE_DIR, "data/fasta/all_homologs.fasta"), "fasta"))
    msa_records = list(SeqIO.parse(os.path.join(BASE_DIR, "data/msa/all_msa.fasta"), "fasta"))

    stn7_rec = next(r for r in records if "STN7_Arabidopsis" in r.id)
    stn8_rec = next(r for r in records if "STN8_Arabidopsis" in r.id)
    stn7_seq = str(stn7_rec.seq).upper().replace("-", "")
    stn8_seq = str(stn8_rec.seq).upper().replace("-", "")

    stn7_pdb = os.path.join(BASE_DIR, "data/structures/Q9S713_AlphaFold.pdb")
    stn8_pdb = os.path.join(BASE_DIR, "data/structures/Q9LZV4_AlphaFold.pdb")

    stn7_res, stn7_plddt = extract_plddt(stn7_pdb)
    stn8_res, stn8_plddt = extract_plddt(stn8_pdb)

    stn7_kd = kd_profile(stn7_seq)
    stn8_kd = kd_profile(stn8_seq)

    entropy = shannon_entropy_per_col(msa_records)
    rolling_avg = pd.Series(entropy).rolling(window=15, min_periods=1, center=True).mean()

    identity_df = pd.read_csv(os.path.join(RESULTS_DIR, "sequence_identity_matrix.csv"), index_col=0)
    identity_df.index = [x.split()[0].replace("STN7:", "STN7\n").replace("STN8:", "STN8\n")
                          for x in identity_df.index]
    identity_df.columns = identity_df.index

    stn7_kin_comp = aa_composition(stn7_seq[133:452])
    stn8_kin_comp = aa_composition(stn8_seq[132:477])

    # ── Build Figure ──────────────────────────────────────────
    fig = plt.figure(figsize=(20, 22))
    fig.patch.set_facecolor('#f8fafc')

    gs = gridspec.GridSpec(4, 2, figure=fig, hspace=0.45, wspace=0.3,
                           top=0.94, bottom=0.05, left=0.07, right=0.97)

    ax_A = fig.add_subplot(gs[0, 0])   # Sequence Identity Heatmap
    ax_B = fig.add_subplot(gs[0, 1])   # Shannon Entropy
    ax_C1 = fig.add_subplot(gs[1, 0])  # pLDDT STN7
    ax_C2 = fig.add_subplot(gs[1, 1])  # pLDDT STN8
    ax_D = fig.add_subplot(gs[2, 0])   # Hydrophobicity
    ax_E = fig.add_subplot(gs[2, 1])   # Kinase CA Distance
    ax_F = fig.add_subplot(gs[3, :])   # AA Composition full width

    panel_labels = {'A': ax_A, 'B': ax_B, 'C': ax_C1, 'D': ax_D, 'E': ax_E, 'F': ax_F}
    for label, ax in panel_labels.items():
        ax.text(-0.08, 1.08, label, transform=ax.transAxes,
                fontsize=16, fontweight='bold', va='top', color='#1e293b')

    # ── Panel A: Sequence Identity Heatmap ────────────────────
    mask = np.zeros_like(identity_df.values, dtype=bool)
    cmap_id = LinearSegmentedColormap.from_list("id_cmap",
        ["#f8fafc", "#bfdbfe", "#3b82f6", "#1d4ed8", "#1e3a8a"])
    sns.heatmap(identity_df, ax=ax_A, annot=True, fmt=".0f",
                cmap=cmap_id, linewidths=0.4, square=True,
                cbar_kws={"shrink": 0.7, "label": "Identity %"},
                annot_kws={"size": 6.5})
    ax_A.set_title("Pairwise Sequence Identity Matrix (%)", fontsize=11, fontweight='bold', pad=10)
    ax_A.tick_params(axis='both', labelsize=7)
    plt.setp(ax_A.get_xticklabels(), rotation=40, ha='right')
    plt.setp(ax_A.get_yticklabels(), rotation=0)

    # ── Panel B: Shannon Entropy ─────────────────────────────
    ax_B.fill_between(range(len(entropy)), entropy, alpha=0.25, color='#8b5cf6')
    ax_B.plot(entropy, color='#8b5cf6', linewidth=0.6, alpha=0.7)
    ax_B.plot(rolling_avg, color='#db2777', linewidth=1.8, label='Trend (w=15)')
    n_conserved = np.sum(entropy == 0.0)
    ax_B.axhline(y=0, color='#10b981', linestyle='--', linewidth=1.0, label=f'H=0 ({n_conserved} cols)')
    ax_B.set_title("Evolutionary Conservation\n(Shannon Entropy, H)", fontsize=11, fontweight='bold')
    ax_B.set_xlabel("MSA Column Position", fontsize=9)
    ax_B.set_ylabel("H [0 = Conserved]", fontsize=9)
    ax_B.legend(fontsize=8, loc='upper right')
    ax_B.grid(True, linestyle=':', alpha=0.4)
    ax_B.set_facecolor('#fafafa')

    # ── Panel C: pLDDT (shared y) ────────────────────────────
    conf_colors = {'>90': '#1d4ed8', '70-90': '#3b82f6', '50-70': '#eab308', '<50': '#ef4444'}
    for (ax_c, res, plddt, name, kin_s, kin_e) in [
        (ax_C1, stn7_res, stn7_plddt, "STN7", 134, 452),
        (ax_C2, stn8_res, stn8_plddt, "STN8", 133, 477)
    ]:
        ax_c.axvspan(kin_s, kin_e, color='#e0f2fe', alpha=0.6, label='Kinase domain')
        for mask_fn, col, lbl in [
            (plddt > 90, '#1d4ed8', '>90'), ((plddt >= 70) & (plddt <= 90), '#3b82f6', '70-90'),
            ((plddt >= 50) & (plddt < 70), '#eab308', '50-70'), (plddt < 50, '#ef4444', '<50')
        ]:
            ax_c.scatter(res[mask_fn], plddt[mask_fn], s=2.5, color=col, alpha=0.85)
        ax_c.plot(res, plddt, color='#374151', linewidth=0.5, alpha=0.4)
        for y, col in [(90, '#10b981'), (70, '#3b82f6'), (50, '#ef4444')]:
            ax_c.axhline(y=y, color=col, linestyle=':', alpha=0.6, linewidth=0.8)
        ax_c.set_title(f"{name} AlphaFold pLDDT Profile", fontsize=10, fontweight='bold')
        ax_c.set_xlabel("Residue Position", fontsize=8)
        ax_c.set_ylabel("pLDDT", fontsize=8)
        ax_c.set_ylim(0, 105)
        ax_c.grid(True, linestyle=':', alpha=0.35)
        ax_c.set_facecolor('#fafafa')
        if ax_c is ax_C1:
            ax_c.legend(fontsize=7, loc='lower left')

    # ── Panel D: Hydrophobicity ──────────────────────────────
    ax_D.plot(range(1, len(stn7_kd)+1), stn7_kd, color='#0f766e',
              linewidth=1.5, label='STN7', alpha=0.85)
    ax_D.plot(range(1, len(stn8_kd)+1), stn8_kd, color='#f43f5e',
              linewidth=1.5, label='STN8', alpha=0.75)
    ax_D.axhline(y=1.6, color='#dc2626', linestyle='--', linewidth=1.0,
                 alpha=0.8, label='TM threshold (1.6)')
    ax_D.axhline(y=0, color='#9ca3af', linestyle=':', alpha=0.5, linewidth=0.8)
    ax_D.set_title("Kyte-Doolittle Hydrophobicity\n(Sliding Window w=19)", fontsize=10, fontweight='bold')
    ax_D.set_xlabel("Amino Acid Position", fontsize=9)
    ax_D.set_ylabel("Hydrophobicity Score", fontsize=9)
    ax_D.legend(fontsize=8)
    ax_D.grid(True, linestyle=':', alpha=0.35)
    ax_D.set_facecolor('#fafafa')

    # ── Panel E: CA Distance Profile ────────────────────────
    dist_path = os.path.join(RESULTS_DIR, "kinase_structural_rmsd.txt")
    # Re-read the RMSD from the saved file
    rmsd_value = None
    with open(dist_path) as f:
        for line in f:
            if "RMSD:" in line:
                rmsd_value = float(line.strip().split()[-2])
                break

    # Read aligned structures and compute distances again for the plot
    parser = PDB.PDBParser(QUIET=True)
    from Bio import Align
    def get_kin_seq_resnums(pdb_path, kin_start, kin_end):
        s = parser.get_structure("p", pdb_path)
        seq, res = [], []
        for r in s[0]['A']:
            if PDB.Polypeptide.is_aa(r):
                rn = r.get_resname().upper()
                if rn in D3_TO_1:
                    seq.append(D3_TO_1[rn])
                    res.append(r.get_id()[1])
        kin_seq = [seq[res.index(i)] for i in range(kin_start, kin_end+1) if i in res]
        return "".join(kin_seq), res

    stn7_kin_seq_str, stn7_all_res = get_kin_seq_resnums(stn7_pdb, 134, 452)
    stn8_kin_seq_str, stn8_all_res = get_kin_seq_resnums(stn8_pdb, 133, 477)

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    best_aln = aligner.align(stn7_kin_seq_str, stn8_kin_seq_str)[0]

    stn7_struct = parser.get_structure("s7", stn7_pdb)
    stn8_struct = parser.get_structure("s8", stn8_pdb)

    ref_coords_aln, query_coords_aln = best_aln.aligned
    mapped_7, mapped_8 = [], []
    for r_range, q_range in zip(ref_coords_aln, query_coords_aln):
        for idx in range(r_range[1] - r_range[0]):
            mapped_7.append(int(134 + r_range[0] + idx))
            mapped_8.append(int(133 + q_range[0] + idx))

    s7_atoms, s8_atoms = [], []
    for r7, r8 in zip(mapped_7, mapped_8):
        try:
            s7_atoms.append(stn7_struct[0]['A'][r7]['CA'])
            s8_atoms.append(stn8_struct[0]['A'][r8]['CA'])
        except KeyError:
            pass

    sup = PDB.Superimposer()
    sup.set_atoms(s7_atoms, s8_atoms)
    distances = [np.linalg.norm(a.get_coord() - b.get_coord())
                 for a, b in zip(s7_atoms, s8_atoms)]

    ax_E.fill_between(range(len(distances)), distances, alpha=0.2, color='#0f766e')
    ax_E.plot(distances, color='#0f766e', linewidth=1.2, label='C-α Distance')
    ax_E.axhline(y=np.mean(distances), color='#be123c', linestyle='--',
                 linewidth=1.2, label=f'Mean: {np.mean(distances):.2f} Å')
    ax_E.set_title(f"Kinase Domain Structural Distance Profile\n(RMSD = {sup.rms:.4f} Å over {len(distances)} C-α pairs)",
                   fontsize=10, fontweight='bold')
    ax_E.set_xlabel("Aligned Residue Pair", fontsize=9)
    ax_E.set_ylabel("C-α Distance (Å)", fontsize=9)
    ax_E.legend(fontsize=8)
    ax_E.grid(True, linestyle=':', alpha=0.35)
    ax_E.set_facecolor('#fafafa')

    # ── Panel F: AA Composition (Kinase Domain) ──────────────
    x = np.arange(len(AMINO_ACIDS_ORDER))
    width = 0.35
    AA_GROUPS = {
        "Nonpolar": list("AVILMFWPG"),
        "Polar":    list("STCYNQ"),
        "Positive": list("KRH"),
        "Negative": list("DE")
    }
    GROUP_COLORS = {"Nonpolar": "#64748b", "Polar": "#3b82f6",
                    "Positive": "#f43f5e", "Negative": "#f59e0b"}

    def get_group_color(aa):
        for g, aas in AA_GROUPS.items():
            if aa in aas:
                return GROUP_COLORS[g]
        return "#94a3b8"

    bar_colors = [get_group_color(aa) for aa in AMINO_ACIDS_ORDER]
    v7 = [stn7_kin_comp[aa] for aa in AMINO_ACIDS_ORDER]
    v8 = [stn8_kin_comp[aa] for aa in AMINO_ACIDS_ORDER]

    ax_F.bar(x - width/2, v7, width, label="STN7 Kinase Domain", alpha=0.85,
             color=bar_colors, edgecolor='white', linewidth=0.4)
    ax_F.bar(x + width/2, v8, width, label="STN8 Kinase Domain", alpha=0.55,
             color=bar_colors, edgecolor='white', linewidth=0.4, hatch='//')

    ax_F.set_xticks(x)
    ax_F.set_xticklabels(AMINO_ACIDS_ORDER, fontsize=10)
    ax_F.set_xlabel("Amino Acid (Single-letter Code)", fontsize=10)
    ax_F.set_ylabel("Mole %", fontsize=10)
    ax_F.set_title("Amino Acid Composition — Stromal Kinase Domains",
                   fontsize=11, fontweight='bold')
    ax_F.legend(fontsize=9)
    ax_F.grid(True, axis='y', linestyle=':', alpha=0.4)
    ax_F.spines['top'].set_visible(False)
    ax_F.spines['right'].set_visible(False)
    ax_F.set_facecolor('#fafafa')

    # Physicochemical group labels
    for group, aas in AA_GROUPS.items():
        positions = [AMINO_ACIDS_ORDER.index(a) for a in aas if a in AMINO_ACIDS_ORDER]
        if positions:
            ax_F.annotate('', xy=(max(positions) + 0.5, -1.6),
                          xytext=(min(positions) - 0.5, -1.6),
                          xycoords='data',
                          arrowprops=dict(arrowstyle='-', color=GROUP_COLORS[group], lw=2.0))
            ax_F.text(np.mean(positions), -2.4, group, ha='center', va='top',
                      fontsize=8, color=GROUP_COLORS[group], fontstyle='italic')

    # ── Super-title ──────────────────────────────────────────
    fig.text(0.5, 0.975,
             "Comparative Structural Bioinformatics of STN7 and STN8 Chloroplast Thylakoid Kinases",
             ha='center', va='top', fontsize=14, fontweight='bold', color='#0f172a')
    fig.text(0.5, 0.963,
             "Arabidopsis thaliana | AlphaFold v6 Models | UniProt Q9S713 & Q9LZV4",
             ha='center', va='top', fontsize=9.5, color='#475569', style='italic')

    # ── Save ─────────────────────────────────────────────────
    out_path = os.path.join(FIGURES_DIR, "summary_multipanel.png")
    plt.savefig(out_path, dpi=300, facecolor=fig.get_facecolor())
    plt.close()
    print(f"Saved multi-panel summary figure to: {out_path}")
    print("=======================================================")

if __name__ == "__main__":
    main()
