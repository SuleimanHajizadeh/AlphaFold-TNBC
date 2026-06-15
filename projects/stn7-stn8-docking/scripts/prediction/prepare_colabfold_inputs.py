#!/usr/bin/env python3
"""
STN7 and STN8 ColabFold Query Input Preparer
Reads the downloaded Arabidopsis STN7 and STN8 sequences from alignment / fasta files,
cleans them of any gaps, and writes ready-to-run FASTA query files for:
  - STN7 Monomer
  - STN8 Monomer
  - STN7-STN8 Heterodimer
  - STN7 Homodimer
  - STN8 Homodimer
Query files are stored under data/prediction_queries/.
"""

import os
from Bio import SeqIO

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    fasta_path = os.path.join(base_dir, "data/fasta/all_homologs.fasta")
    output_dir = os.path.join(base_dir, "data/prediction_queries")
    os.makedirs(output_dir, exist_ok=True)
    
    print("=== Preparing ColabFold Query Inputs ===")
    
    # 1. Parse alignment FASTA to get sequences
    if not os.path.exists(fasta_path):
        print(f"Error: Base FASTA not found at {fasta_path}")
        return
        
    records = list(SeqIO.parse(fasta_path, "fasta"))
    
    # Extract Arabidopsis sequences and remove alignment gaps (-)
    stn7_rec = next((r for r in records if "STN7_Arabidopsis" in r.id), None)
    stn8_rec = next((r for r in records if "STN8_Arabidopsis" in r.id), None)
    
    if not stn7_rec or not stn8_rec:
        print("Error: Could not find Arabidopsis sequences in all_homologs.fasta")
        return
        
    stn7_seq = str(stn7_rec.seq).upper().replace("-", "")
    stn8_seq = str(stn8_rec.seq).upper().replace("-", "")
    
    print(f"Parsed STN7 (Q9S713) sequence length: {len(stn7_seq)} residues")
    print(f"Parsed STN8 (Q9LZV4) sequence length: {len(stn8_seq)} residues")
    
    # 2. Write Query Files
    # A. STN7 Monomer
    stn7_mono_path = os.path.join(output_dir, "stn7_monomer.fasta")
    with open(stn7_mono_path, "w") as f:
        f.write(f">STN7_Arabidopsis_Q9S713\n{stn7_seq}\n")
    print(f"Saved STN7 Monomer query → {stn7_mono_path}")
    
    # B. STN8 Monomer
    stn8_mono_path = os.path.join(output_dir, "stn8_monomer.fasta")
    with open(stn8_mono_path, "w") as f:
        f.write(f">STN8_Arabidopsis_Q9LZV4\n{stn8_seq}\n")
    print(f"Saved STN8 Monomer query → {stn8_mono_path}")
    
    # C. STN7-STN8 Heterodimer (Multiple entry format & colon-separated format)
    stn7_stn8_hetero_path = os.path.join(output_dir, "stn7_stn8_heterodimer.fasta")
    with open(stn7_stn8_hetero_path, "w") as f:
        # ColabFold multi-chain format (colon separated in a single entry or multi-entry)
        f.write(f">STN7_STN8_heterodimer\n{stn7_seq}:{stn8_seq}\n")
    print(f"Saved STN7-STN8 Heterodimer query → {stn7_stn8_hetero_path}")
    
    # D. STN7 Homodimer
    stn7_homo_path = os.path.join(output_dir, "stn7_homodimer.fasta")
    with open(stn7_homo_path, "w") as f:
        f.write(f">STN7_homodimer\n{stn7_seq}:{stn7_seq}\n")
    print(f"Saved STN7 Homodimer query → {stn7_homo_path}")
    
    # E. STN8 Homodimer
    stn8_homo_path = os.path.join(output_dir, "stn8_homodimer.fasta")
    with open(stn8_homo_path, "w") as f:
        f.write(f">STN8_homodimer\n{stn8_seq}:{stn8_seq}\n")
    print(f"Saved STN8 Homodimer query → {stn8_homo_path}")
    
    print("\n=== ColabFold Query Input Preparation Completed ===")
    
if __name__ == "__main__":
    main()
