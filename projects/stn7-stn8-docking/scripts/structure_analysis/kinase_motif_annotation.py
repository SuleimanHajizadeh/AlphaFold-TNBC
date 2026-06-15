#!/usr/bin/env python3
"""
STN7 and STN8 Kinase Functional Motif Annotation Script
Scans the protein sequences of STN7 and STN8 for conserved kinase motifs:
  - G-loop (P-loop): GxGxxG - ATP binding
  - C-loop (catalytic HRD): HRD/YRD - catalytic base
  - DFG motif (activation loop start): DFG - Mg2+ coordination
  - APE motif (activation loop end): APE
Annotates exact positions and generates a visual motif map.
"""

import os
import re
import json
import requests
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO

# ============================================================
# Kinase motif definitions (regex patterns)
# ============================================================
KINASE_MOTIFS = {
    "G-loop (P-loop)": {
        "pattern": r"G.G..G",
        "color": "#3b82f6",
        "description": "ATP-binding phosphate clamp"
    },
    "HRD motif (C-loop)": {
        "pattern": r"[HY]RD",
        "color": "#ef4444",
        "description": "Catalytic aspartate base"
    },
    "DFG motif": {
        "pattern": r"DF[GS]",
        "color": "#f59e0b",
        "description": "Activation loop / Mg2+ coordination"
    },
    "APE motif": {
        "pattern": r"APE",
        "color": "#10b981",
        "description": "Activation loop end / substrate binding"
    },
    "GXGXXG (alternative)": {
        "pattern": r"G.G..[AG]",
        "color": "#8b5cf6",
        "description": "Alternative ATP-binding variation"
    }
}

def find_motifs(sequence, motifs):
    """
    Finds all occurrences of each kinase motif in the sequence.
    Returns a dict mapping motif_name -> list of (start, end, matched_sequence)
    All positions are 1-indexed for biological convention.
    """
    results = {}
    for motif_name, motif_data in motifs.items():
        pattern = motif_data["pattern"]
        matches = []
        for m in re.finditer(pattern, sequence):
            start = m.start() + 1  # 1-indexed
            end = m.end()           # 1-indexed inclusive
            matches.append({
                "start": start,
                "end": end,
                "sequence": m.group()
            })
        results[motif_name] = matches
    return results

def draw_motif_map(seq_name, sequence, motif_results, ax, seq_len_max):
    """
    Draws a linear motif diagram for a given sequence on a matplotlib Axes.
    """
    seq_len = len(sequence)
    
    # Draw backbone
    ax.barh(0, seq_len, left=0, height=0.3, color="#e2e8f0", edgecolor="#94a3b8", linewidth=0.5)
    
    # Mark kinase domain region
    if "STN7" in seq_name:
        kin_start, kin_end = 134, 452
    else:
        kin_start, kin_end = 133, 477
    
    ax.barh(0, kin_end - kin_start, left=kin_start - 1, height=0.3,
            color="#dbeafe", edgecolor="#3b82f6", linewidth=0.8, alpha=0.6,
            label="Kinase Domain")
    
    # Draw motifs as blocks
    y_levels = [0.25, 0.5, 0.75, 1.0, 1.25]  # Stack overlapping motifs
    motif_occupancy = {}
    
    for motif_name, matches in motif_results.items():
        color = KINASE_MOTIFS[motif_name]["color"]
        for match in matches:
            start = match["start"] - 1
            width = match["end"] - match["start"] + 1
            seq_match = match["sequence"]
            pos = match["start"]
            
            # Find an available level to avoid overlap
            level = 0.45
            for y in y_levels:
                occupied = False
                for prev_start, prev_end in motif_occupancy.get(y, []):
                    if not (start > prev_end or (start + width) < prev_start):
                        occupied = True
                        break
                if not occupied:
                    level = y
                    if y not in motif_occupancy:
                        motif_occupancy[y] = []
                    motif_occupancy[y].append((start, start + width))
                    break
            
            ax.barh(level, width, left=start, height=0.22,
                    color=color, edgecolor="white", linewidth=0.5, alpha=0.9)
            
            # Annotate with matched sequence if there's room
            if width > 5:
                ax.text(start + width / 2, level, seq_match,
                        ha='center', va='center', fontsize=6.5,
                        fontweight='bold', color='white')
    
    ax.set_xlim(0, seq_len_max)
    ax.set_ylim(-0.3, 1.6)
    ax.set_title(seq_name, fontsize=11, fontweight='bold', pad=8)
    ax.set_xlabel("Amino Acid Position", fontsize=9)
    ax.set_yticks([])
    ax.grid(True, axis='x', linestyle=':', alpha=0.4)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    fasta_path = os.path.join(base_dir, "data/fasta/all_homologs.fasta")
    figures_dir = os.path.join(base_dir, "figures")
    results_dir = os.path.join(base_dir, "results")
    
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    print("=== Kinase Functional Motif Annotation ===\n")
    
    # Read all sequences
    records = list(SeqIO.parse(fasta_path, "fasta"))
    
    # Focus on Arabidopsis STN7 and STN8 for detailed annotation
    stn7_rec = next(r for r in records if "STN7_Arabidopsis" in r.id)
    stn8_rec = next(r for r in records if "STN8_Arabidopsis" in r.id)
    
    stn7_seq = str(stn7_rec.seq).upper().replace("-", "")
    stn8_seq = str(stn8_rec.seq).upper().replace("-", "")
    
    print(f"Scanning STN7 ({len(stn7_seq)} aa) for kinase motifs...")
    stn7_motifs = find_motifs(stn7_seq, KINASE_MOTIFS)
    
    print(f"Scanning STN8 ({len(stn8_seq)} aa) for kinase motifs...")
    stn8_motifs = find_motifs(stn8_seq, KINASE_MOTIFS)
    
    # ---- Print results ----
    print("\n--- STN7 Kinase Motif Locations ---")
    for motif_name, matches in stn7_motifs.items():
        if matches:
            for m in matches:
                print(f"  {motif_name}: {m['sequence']} @ residues {m['start']}-{m['end']}")
        else:
            print(f"  {motif_name}: Not found")
    
    print("\n--- STN8 Kinase Motif Locations ---")
    for motif_name, matches in stn8_motifs.items():
        if matches:
            for m in matches:
                print(f"  {motif_name}: {m['sequence']} @ residues {m['start']}-{m['end']}")
        else:
            print(f"  {motif_name}: Not found")
    
    # Save motif annotation to JSON
    motif_report = {
        "STN7_Q9S713": stn7_motifs,
        "STN8_Q9LZV4": stn8_motifs
    }
    json_path = os.path.join(results_dir, "kinase_motif_annotation.json")
    with open(json_path, "w") as f:
        json.dump(motif_report, f, indent=2)
    print(f"\nSaved motif annotation JSON to: {json_path}")
    
    # ---- Draw motif map figure ----
    fig, axes = plt.subplots(2, 1, figsize=(12, 6))
    fig.suptitle("STN7 and STN8 Kinase Functional Motif Map\n(Arabidopsis thaliana)",
                 fontsize=13, fontweight='bold', y=1.01)
    
    max_len = max(len(stn7_seq), len(stn8_seq))
    
    draw_motif_map("STN7 (Q9S713) — 562 aa", stn7_seq, stn7_motifs, axes[0], max_len)
    draw_motif_map("STN8 (Q9LZV4) — 495 aa", stn8_seq, stn8_motifs, axes[1], max_len)
    
    # ---- Create legend ----
    legend_patches = [
        mpatches.Patch(color="#dbeafe", edgecolor="#3b82f6", label="Stromal Kinase Domain")
    ]
    for motif_name, motif_data in KINASE_MOTIFS.items():
        legend_patches.append(
            mpatches.Patch(color=motif_data["color"], label=f"{motif_name} ({motif_data['description']})")
        )
    fig.legend(handles=legend_patches, loc='lower center', ncol=3,
               fontsize=8, framealpha=0.9, bbox_to_anchor=(0.5, -0.08))
    
    plt.tight_layout()
    motif_fig_path = os.path.join(figures_dir, "kinase_motif_map.png")
    plt.savefig(motif_fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved kinase motif map to: {motif_fig_path}")
    print("==========================================")

if __name__ == "__main__":
    main()
