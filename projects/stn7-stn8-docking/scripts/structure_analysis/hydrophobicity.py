#!/usr/bin/env python3
"""
STN7 and STN8 Kyte-Doolittle Hydrophobicity Profiler
Computes the Kyte-Doolittle hydrophobicity profile using a sliding window of 19 residues
to identify and compare potential transmembrane (TM) or membrane-association helices.
Saves profiles to results/ and generates a publication-quality comparative plot.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

# Standard Kyte-Doolittle Hydrophobicity Scale for single-letter Amino Acids
KD_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

def calculate_kd_profile(sequence, window_size=19):
    """
    Computes Kyte-Doolittle hydrophobicity score using a sliding window.
    """
    L = len(sequence)
    profile = np.zeros(L)
    half_w = window_size // 2
    
    # Calculate window averages
    for i in range(L):
        # Determine window boundaries
        start = max(0, i - half_w)
        end = min(L, i + half_w + 1)
        window_seq = sequence[start:end]
        
        # Calculate average of known amino acids in the window
        scores = [KD_SCALE[aa] for aa in window_seq if aa in KD_SCALE]
        profile[i] = np.mean(scores) if scores else 0.0
        
    return profile

def identify_tm_segments(profile, threshold=1.6, min_len=15):
    """
    Identifies contiguous residues crossing the hydrophobicity threshold
    indicating potential transmembrane segments.
    """
    segments = []
    in_segment = False
    start = 0
    
    for i, val in enumerate(profile):
        if val >= threshold:
            if not in_segment:
                start = i + 1 # 1-based index
                in_segment = True
        else:
            if in_segment:
                end = i
                if (end - start + 1) >= min_len:
                    segments.append((start, end, np.mean(profile[start-1:end])))
                in_segment = False
                
    # Handle end of sequence
    if in_segment:
        end = len(profile)
        if (end - start + 1) >= min_len:
            segments.append((start, end, np.mean(profile[start-1:end])))
            
    return segments

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    fasta_path = os.path.join(base_dir, "data/fasta/all_homologs.fasta")
    
    figures_dir = os.path.join(base_dir, "figures")
    results_dir = os.path.join(base_dir, "results")
    
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    print("=== Kyte-Doolittle Hydrophobicity & Transmembrane Analysis ===")
    
    # Read Arabidopsis sequences
    records = list(SeqIO.parse(fasta_path, "fasta"))
    stn7_rec = next(rec for rec in records if "STN7_Arabidopsis" in rec.id)
    stn8_rec = next(rec for rec in records if "STN8_Arabidopsis" in rec.id)
    
    stn7_seq = str(stn7_rec.seq).upper().replace('-', '') # Remove any alignment gaps
    stn8_seq = str(stn8_rec.seq).upper().replace('-', '')
    
    print(f"STN7 pure sequence length: {len(stn7_seq)} aa")
    print(f"STN8 pure sequence length: {len(stn8_seq)} aa")
    
    # 1. Compute profiles
    window_size = 19
    stn7_kd = calculate_kd_profile(stn7_seq, window_size)
    stn8_kd = calculate_kd_profile(stn8_seq, window_size)
    
    # 2. Identify potential TM regions
    tm_threshold = 1.6
    stn7_tms = identify_tm_segments(stn7_kd, threshold=tm_threshold)
    stn8_tms = identify_tm_segments(stn8_kd, threshold=tm_threshold)
    
    print("\n--- Predicted Transmembrane / Hydrophobic Segments (Threshold > 1.6) ---")
    print("STN7 Predicted TM Regions:")
    for tm in stn7_tms:
        print(f"  Residues {tm[0]}-{tm[1]} | Length: {tm[1]-tm[0]+1} aa | Mean Hydrophobicity: {tm[2]:.2f}")
    if not stn7_tms:
        # Fallback to absolute maximum to see highest hydrophobic region
        max_idx = np.argmax(stn7_kd)
        print(f"  No segments > 1.6. Highest hydrophobic peak at residues {max_idx-9}-{max_idx+9} (Score: {stn7_kd[max_idx]:.2f})")
        
    print("STN8 Predicted TM Regions:")
    for tm in stn8_tms:
        print(f"  Residues {tm[0]}-{tm[1]} | Length: {tm[1]-tm[0]+1} aa | Mean Hydrophobicity: {tm[2]:.2f}")
    if not stn8_tms:
        max_idx = np.argmax(stn8_kd)
        print(f"  No segments > 1.6. Highest hydrophobic peak at residues {max_idx-9}-{max_idx+9} (Score: {stn8_kd[max_idx]:.2f})")
        
    # Save profiles to CSV
    # Since they have different lengths, we pad or write separately. Let's write separately.
    pd.DataFrame({"Residue": range(1, len(stn7_kd)+1), "Hydrophobicity": stn7_kd}).to_csv(os.path.join(results_dir, "stn7_hydrophobicity.csv"), index=False)
    pd.DataFrame({"Residue": range(1, len(stn8_kd)+1), "Hydrophobicity": stn8_kd}).to_csv(os.path.join(results_dir, "stn8_hydrophobicity.csv"), index=False)
    
    # 3. Plot profiles
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharey=True)
    
    # STN7 Plot
    axes[0].plot(range(1, len(stn7_kd)+1), stn7_kd, color='#0f766e', linewidth=1.5, label='Hydrophobicity (w=19)')
    axes[0].axhline(y=tm_threshold, color='#be123c', linestyle='--', linewidth=1.2, label=f'TM Threshold ({tm_threshold})')
    axes[0].axhline(y=0, color='#9ca3af', linestyle=':', alpha=0.7)
    
    # Highlight predicted TMs
    for tm in stn7_tms:
        axes[0].axvspan(tm[0], tm[1], color='#f59e0b', alpha=0.3, label='Predicted TMD')
    # If no TMD, highlight the peak
    if not stn7_tms:
        max_idx = np.argmax(stn7_kd)
        axes[0].axvspan(max_idx-9, max_idx+9, color='#ef4444', alpha=0.15, label='Hydrophobic Peak (Non-TM)')
        
    axes[0].set_title("STN7 Chloroplast Kinase Kyte-Doolittle Hydrophobicity Profile", fontsize=12, fontweight='bold')
    axes[0].set_ylabel("Hydrophobicity Score", fontsize=10)
    axes[0].grid(True, linestyle=':', alpha=0.5)
    axes[0].legend(loc='lower left', fontsize=9)
    
    # STN8 Plot
    axes[1].plot(range(1, len(stn8_kd)+1), stn8_kd, color='#f43f5e', linewidth=1.5, label='Hydrophobicity (w=19)')
    axes[1].axhline(y=tm_threshold, color='#be123c', linestyle='--', linewidth=1.2, label=f'TM Threshold ({tm_threshold})')
    axes[1].axhline(y=0, color='#9ca3af', linestyle=':', alpha=0.7)
    
    # Highlight predicted TMs
    for tm in stn8_tms:
        axes[1].axvspan(tm[0], tm[1], color='#f59e0b', alpha=0.3, label='Predicted TMD')
    if not stn8_tms:
        max_idx = np.argmax(stn8_kd)
        axes[1].axvspan(max_idx-9, max_idx+9, color='#ef4444', alpha=0.15, label='Hydrophobic Peak (Non-TM)')
        
    axes[1].set_title("STN8 Chloroplast Kinase Kyte-Doolittle Hydrophobicity Profile", fontsize=12, fontweight='bold')
    axes[1].set_xlabel("Amino Acid Position", fontsize=10)
    axes[1].set_ylabel("Hydrophobicity Score", fontsize=10)
    axes[1].grid(True, linestyle=':', alpha=0.5)
    axes[1].legend(loc='lower left', fontsize=9)
    
    plt.tight_layout()
    plot_path = os.path.join(figures_dir, "hydrophobicity_profile.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    print(f"\nSaved hydrophobicity profile plot to: {plot_path}")
    print("====================================================")

if __name__ == "__main__":
    main()
