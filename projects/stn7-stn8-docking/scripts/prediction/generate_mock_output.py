#!/usr/bin/env python3
"""
STN7-STN8 Dimer Mock ColabFold Output Generator
Creates a simulated ColabFold output directory to test and verify
the post-prediction dimer interface analysis script.
Merges the aligned kinase domains of STN7 (Chain A) and STN8 (Chain B) into a single PDB
and writes a mock scores/PAE JSON file.
"""

import os
import json
import numpy as np
from Bio import PDB

def merge_pdbs(pdb7_path, pdb8_path, out_pdb_path):
    """
    Reads pdb7 (fixed) and pdb8 (moving) and merges them into a single PDB.
    Assigns Chain A to STN7 and Chain B to STN8.
    """
    parser = PDB.PDBParser(QUIET=True)
    struct7 = parser.get_structure("STN7", pdb7_path)
    struct8 = parser.get_structure("STN8", pdb8_path)
    
    # Create new PDB structure
    new_struct = PDB.Structure.Structure("mock_dimer")
    new_model = PDB.Model.Model(0)
    new_struct.add(new_model)
    
    # Extract chains
    model7 = struct7[0]
    model8 = struct8[0]
    
    # We assume first chain in both is the target
    chain7 = list(model7.get_chains())[0]
    chain8 = list(model8.get_chains())[0]
    
    # Clone and rename chain7 to A
    chainA = PDB.Chain.Chain("A")
    for residue in chain7:
        if PDB.Polypeptide.is_aa(residue):
            # Clone residue and atoms to avoid modification side effects
            res_clone = residue.copy()
            chainA.add(res_clone)
            
    # Clone and rename chain8 to B
    chainB = PDB.Chain.Chain("B")
    for residue in chain8:
        if PDB.Polypeptide.is_aa(residue):
            res_clone = residue.copy()
            chainB.add(res_clone)
            
    new_model.add(chainA)
    new_model.add(chainB)
    
    # Save the merged PDB
    io = PDB.PDBIO()
    io.set_structure(new_struct)
    io.save(out_pdb_path)
    print(f"Merged PDB saved → {out_pdb_path}")
    return len(chainA), len(chainB)

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    pdb7_path = os.path.join(base_dir, "data/structures/STN7_aligned.pdb")
    pdb8_path = os.path.join(base_dir, "data/structures/STN8_aligned.pdb")
    
    output_dir = os.path.join(base_dir, "results/prediction_outputs/stn7_stn8_heterodimer")
    os.makedirs(output_dir, exist_ok=True)
    
    print("=== Generating Mock ColabFold Output ===")
    
    if not os.path.exists(pdb7_path) or not os.path.exists(pdb8_path):
        print(f"Error: Aligned PDB files not found under data/structures/.")
        print("Please ensure Stage 4 structure_compare.py has been run.")
        return
        
    pdb_filename = "STN7_STN8_predicted_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb"
    json_filename = "STN7_STN8_predicted_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json"
    
    out_pdb_path = os.path.join(output_dir, pdb_filename)
    out_json_path = os.path.join(output_dir, json_filename)
    
    # 1. Merge structures and get chain lengths
    len_a, len_b = merge_pdbs(pdb7_path, pdb8_path, out_pdb_path)
    total_len = len_a + len_b
    print(f"Chain A (STN7 Kinase Domain) length: {len_a} residues")
    print(f"Chain B (STN8 Kinase Domain) length: {len_b} residues")
    print(f"Total complex length: {total_len} residues")
    
    # 2. Generate simulated PAE and pLDDT data
    # Create mock pLDDT scores (high in domain cores, lower in loops)
    plddt = np.random.uniform(85.0, 95.0, total_len).tolist()
    
    # Create mock PAE matrix
    # Intra-chain PAE is low (e.g. 1.0 to 4.0)
    # Inter-chain PAE is higher (e.g., 20.0 to 30.0 for non-contacting, 5.0 to 12.0 for contacting interface residues)
    pae = np.zeros((total_len, total_len))
    
    # Populate Intra-chain A
    for i in range(len_a):
        for j in range(len_a):
            pae[i, j] = abs(i - j) * 0.05 + np.random.uniform(1.0, 3.0)
            
    # Populate Intra-chain B
    for i in range(len_b):
        for j in range(len_b):
            idx_i = len_a + i
            idx_j = len_a + j
            pae[idx_i, idx_j] = abs(i - j) * 0.05 + np.random.uniform(1.0, 3.0)
            
    # Populate Inter-chain PAE (A to B and B to A)
    # Set background inter-chain PAE to high (no confidence in alignment)
    pae[0:len_a, len_a:total_len] = np.random.uniform(22.0, 28.0, (len_a, len_b))
    pae[len_a:total_len, 0:len_a] = np.random.uniform(22.0, 28.0, (len_b, len_a))
    
    # Simulate a contacting interface region with high confidence (low PAE)
    # Let's say residues 100-130 in Chain A interact with residues 120-150 in Chain B
    for i in range(100, 131):
        for j in range(120, 151):
            idx_j = len_a + j
            score = np.random.uniform(4.0, 9.0)
            pae[i, idx_j] = score
            pae[idx_j, i] = score
            
    # Ensure PAE is bounded (AlphaFold standard is capped at 31.75 Å)
    pae = np.clip(pae, 0.0, 31.75).tolist()
    
    # Save scores to JSON file matching ColabFold v1.5/v1.6 format
    scores_data = {
        "plddt": plddt,
        "max_pae": 31.75,
        "pae": pae,
        "ptm": 0.82,
        "iptm": 0.75
    }
    
    with open(out_json_path, "w") as f:
        json.dump(scores_data, f, indent=2)
        
    print(f"Simulated scores JSON saved → {out_json_path}")
    print("=== Dimer Mock Generation Completed ===\n")

if __name__ == "__main__":
    main()
