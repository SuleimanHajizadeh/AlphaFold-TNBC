#!/usr/bin/env python3
"""
STN7 vs STN8 Structural Alignment and RMSD Calculator
Aligns the kinase domains of STN7 (residues 134-452) and STN8 (residues 133-477).
Maps matching residues using sequence alignment, superimposes their CA atoms,
calculates exact structural RMSD, and exports aligned PDB coordinates and distance profile.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB, Align

def get_sequence_from_pdb(structure_file, model_id=0, chain_id='A'):
    """
    Extracts the amino acid sequence and residue numbers from a PDB file.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", structure_file)
    model = structure[model_id]
    chain = model[chain_id]
    
    seq = []
    res_nums = []
    
    d3_to_1 = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    for residue in chain:
        if PDB.Polypeptide.is_aa(residue):
            res_name = residue.get_resname().upper()
            if res_name in d3_to_1:
                seq.append(d3_to_1[res_name])
                res_nums.append(residue.get_id()[1])
                
    return "".join(seq), res_nums

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    stn7_pdb_path = os.path.join(base_dir, "data/structures/Q9S713_AlphaFold.pdb")
    stn8_pdb_path = os.path.join(base_dir, "data/structures/Q9LZV4_AlphaFold.pdb")
    
    figures_dir = os.path.join(base_dir, "figures")
    results_dir = os.path.join(base_dir, "results")
    structures_dir = os.path.join(base_dir, "data/structures")
    
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    print("=== Aligning STN7 & STN8 Kinase Domain Structures ===")
    
    # 1. Extract sequences and residue numbers
    stn7_seq, stn7_res_nums = get_sequence_from_pdb(stn7_pdb_path)
    stn8_seq, stn8_res_nums = get_sequence_from_pdb(stn8_pdb_path)
    
    print(f"STN7 PDB Sequence Length: {len(stn7_seq)} residues")
    print(f"STN8 PDB Sequence Length: {len(stn8_seq)} residues")
    
    # 2. Extract Kinase domain sequence fragments
    # Kinase domains: STN7 (134-452), STN8 (133-477)
    # Map UniProt residue number to PDB structure index
    stn7_kin_range = range(134, 453)
    stn8_kin_range = range(133, 478)
    
    stn7_kin_seq_list = [stn7_seq[stn7_res_nums.index(r)] for r in stn7_kin_range if r in stn7_res_nums]
    stn8_kin_seq_list = [stn8_seq[stn8_res_nums.index(r)] for r in stn8_kin_range if r in stn8_res_nums]
    
    stn7_kin_seq = "".join(stn7_kin_seq_list)
    stn8_kin_seq = "".join(stn8_kin_seq_list)
    
    print(f"STN7 Kinase Domain Fragment: {len(stn7_kin_seq)} residues")
    print(f"STN8 Kinase Domain Fragment: {len(stn8_kin_seq)} residues")
    
    # 3. Perform pairwise alignment of the kinase domain sequences
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    alignments = aligner.align(stn7_kin_seq, stn8_kin_seq)
    best_aln = alignments[0]
    
    print("\nPairwise Sequence Alignment of Kinase Domains computed:")
    print(f"Alignment Score: {best_aln.score:.2f}")
    
    # Map matched columns to matching residues
    ref_aligned_coords, query_aligned_coords = best_aln.aligned
    
    # Construct indices maps
    stn7_kin_mapped_indices = []
    stn8_kin_mapped_indices = []
    
    for r_range, q_range in zip(ref_aligned_coords, query_aligned_coords):
        r_start, r_end = r_range
        q_start, q_end = q_range
        # Add corresponding index mapping (with explicit cast to native python int to avoid numpy.int64 KeyError in Bio.PDB)
        for idx in range(r_end - r_start):
            stn7_kin_mapped_indices.append(int(134 + r_start + idx))
            stn8_kin_mapped_indices.append(int(133 + q_start + idx))
            
    print(f"Mapped {len(stn7_kin_mapped_indices)} matching structural residues for superposition.")
    
    # 4. Superimpose structural C-alpha atoms
    parser = PDB.PDBParser(QUIET=True)
    stn7_struct = parser.get_structure("STN7", stn7_pdb_path)
    stn8_struct = parser.get_structure("STN8", stn8_pdb_path)
    
    stn7_atoms = []
    stn8_atoms = []
    
    # Select CA atoms corresponding to mapped indices
    for stn7_res_num, stn8_res_num in zip(stn7_kin_mapped_indices, stn8_kin_mapped_indices):
        try:
            stn7_atom = stn7_struct[0]['A'][stn7_res_num]['CA']
            stn8_atom = stn8_struct[0]['A'][stn8_res_num]['CA']
            stn7_atoms.append(stn7_atom)
            stn8_atoms.append(stn8_atom)
        except KeyError:
            # Skip if atom or residue missing
            continue
            
    print(f"Retrieved {len(stn7_atoms)} matching C-alpha atoms.")
    
    superimposer = PDB.Superimposer()
    # Align STN8 (moving) to STN7 (fixed)
    superimposer.set_atoms(stn7_atoms, stn8_atoms)
    superimposer.apply(stn8_struct.get_atoms())
    
    print(f"\nStructural Superposition Successful!")
    print(f"Calculated RMSD: {superimposer.rms:.4f} Å")
    
    # 5. Export aligned PDB structures
    io = PDB.PDBIO()
    stn7_out = os.path.join(structures_dir, "STN7_aligned.pdb")
    stn8_out = os.path.join(structures_dir, "STN8_aligned.pdb")
    
    io.set_structure(stn7_struct)
    io.save(stn7_out)
    io.set_structure(stn8_struct)
    io.save(stn8_out)
    
    print(f"Saved aligned STN7 PDB coordinates to: {stn7_out}")
    print(f"Saved aligned STN8 PDB coordinates to: {stn8_out}")
    
    # 6. Calculate CA distance profile
    distances = []
    for a7, a8 in zip(stn7_atoms, stn8_atoms):
        coord7 = a7.get_coord()
        coord8 = a8.get_coord()
        dist = np.linalg.norm(coord7 - coord8)
        distances.append(dist)
        
    # Save distance report
    with open(os.path.join(results_dir, "kinase_structural_rmsd.txt"), "w") as f:
        f.write(f"STN7 vs STN8 Structural Kinase Domain Superposition\n")
        f.write(f"Number of C-alpha atoms used: {len(distances)}\n")
        f.write(f"RMSD: {superimposer.rms:.4f} Angstroms\n")
        f.write(f"Average pairwise C-alpha distance: {np.mean(distances):.4f} Angstroms\n")
        f.write(f"Max C-alpha distance: {np.max(distances):.4f} Angstroms\n")
        
    # Plot distances
    plt.figure(figsize=(10, 5))
    plt.plot(distances, color='#0f766e', linewidth=1.5, label='CA Distance')
    plt.axhline(y=np.mean(distances), color='#be123c', linestyle='--', linewidth=1.2, label=f'Mean Distance ({np.mean(distances):.2f} Å)')
    
    plt.title("STN7 vs STN8 Kinase Domain Local Structural Distance Profile", fontsize=12, fontweight='bold')
    plt.xlabel("Aligned Amino Acid Position", fontsize=10)
    plt.ylabel("C-alpha Distance (Å)", fontsize=10)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    
    dist_plot_path = os.path.join(figures_dir, "kinase_ca_distance_profile.png")
    plt.savefig(dist_plot_path, dpi=300)
    plt.close()
    print(f"Saved distance profile plot to: {dist_plot_path}")
    print("=====================================================")

if __name__ == "__main__":
    main()
