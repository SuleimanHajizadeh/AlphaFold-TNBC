#!/usr/bin/env python3
"""
STN7-STN8 ColabFold Predicted Complex Interface Analyzer
Parses ColabFold dimer prediction outputs (PDB coordinates and scores JSON).
Calculates inter-chain residue-residue distances to map the interaction interface,
extracts and plots Predicted Aligned Error (PAE) matrices and pLDDT confidence scores,
and generates publication-quality visualization figures and a text report.
"""

import os
import sys
import glob
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB

def get_pdb_and_json(input_dir):
    """
    Search for predicted PDB (rank 001) and scores JSON files in the input directory.
    """
    # Look for rank 001 PDB
    pdb_patterns = [
        "*rank_001*.pdb",
        "*ranked_0.pdb",
        "*.pdb"
    ]
    pdb_file = None
    for pat in pdb_patterns:
        hits = glob.glob(os.path.join(input_dir, pat))
        if hits:
            # Sort to get top rank
            hits.sort()
            pdb_file = hits[0]
            break
            
    # Look for scores or PAE JSON
    json_patterns = [
        "*predicted_aligned_error_v1.json",
        "*scores*.json",
        "*.json"
    ]
    json_file = None
    for pat in json_patterns:
        hits = glob.glob(os.path.join(input_dir, pat))
        if hits:
            hits.sort()
            json_file = hits[0]
            break
            
    return pdb_file, json_file

def parse_dimer_structure(pdb_path):
    """
    Parses the dimer structure, extracts coordinates and residue names for Chain A and B.
    """
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("dimer", pdb_path)
    model = struct[0]
    
    chains = list(model.get_chains())
    if len(chains) < 2:
        print(f"Error: Expected at least 2 chains in PDB, found {len(chains)}")
        sys.exit(1)
        
    chain_a = chains[0]
    chain_b = chains[1]
    
    print(f"Analyzing Chain A: ID={chain_a.id}, residues={len(list(chain_a.get_residues()))}")
    print(f"Analyzing Chain B: ID={chain_b.id}, residues={len(list(chain_b.get_residues()))}")
    
    data = {}
    for c_idx, chain in [('A', chain_a), ('B', chain_b)]:
        coords = []
        res_info = [] # List of (res_num, res_name, 1-letter)
        d3_to_1 = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        for residue in chain:
            if PDB.Polypeptide.is_aa(residue) and "CA" in residue:
                coords.append(residue["CA"].get_coord())
                res_name = residue.get_resname().upper()
                one_letter = d3_to_1.get(res_name, 'X')
                res_info.append((residue.get_id()[1], res_name, one_letter))
                
        data[c_idx] = {
            'coords': np.array(coords),
            'res_info': res_info
        }
        
    return data

def analyze_interface(data, threshold):
    """
    Calculates distance matrix between Chain A and Chain B C-alpha atoms.
    Maps interface residues and contact pairs.
    """
    coords_a = data['A']['coords']
    coords_b = data['B']['coords']
    res_a = data['A']['res_info']
    res_b = data['B']['res_info']
    
    n_a = len(coords_a)
    n_b = len(coords_b)
    
    # Calculate distance matrix (n_a x n_b)
    dist_matrix = np.zeros((n_a, n_b))
    for i in range(n_a):
        for j in range(n_b):
            dist_matrix[i, j] = np.linalg.norm(coords_a[i] - coords_b[j])
            
    # Find contact pairs
    contacts = []
    interface_a = set()
    interface_b = set()
    
    for i in range(n_a):
        for j in range(n_b):
            d = dist_matrix[i, j]
            if d <= threshold:
                contacts.append({
                    'idx_a': i,
                    'idx_b': j,
                    'res_num_a': res_a[i][0],
                    'res_name_a': res_a[i][1],
                    'one_letter_a': res_a[i][2],
                    'res_num_b': res_b[j][0],
                    'res_name_b': res_b[j][1],
                    'one_letter_b': res_b[j][2],
                    'distance': d
                })
                interface_a.add(i)
                interface_b.add(j)
                
    return dist_matrix, contacts, list(interface_a), list(interface_b)

def parse_scores_json(json_path):
    """
    Parses scores/PAE JSON file. Supports both ColabFold format structures.
    """
    with open(json_path, 'r') as f:
        scores = json.load(f)
        
    plddt = None
    pae = None
    
    # ColabFold stores scores as a list of dicts or a single dict
    if isinstance(scores, list) and len(scores) > 0:
        data = scores[0]
    else:
        data = scores
        
    # Extract pLDDT
    if 'plddt' in data:
        plddt = np.array(data['plddt'])
        
    # Extract PAE
    if 'pae' in data:
        pae = np.array(data['pae'])
    elif 'predicted_aligned_error' in data:
        pae = np.array(data['predicted_aligned_error'])
        
    return plddt, pae

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    default_input = os.path.join(base_dir, "results/prediction_outputs/stn7_stn8_heterodimer")
    
    parser = argparse.ArgumentParser(description="STN7-STN8 ColabFold Dimer Interface Analyzer")
    parser.add_argument("--input", default=default_input, help="ColabFold output directory")
    parser.add_argument("--threshold", type=float, default=6.0, help="Interface contact distance cutoff in Å (default 6.0)")
    args = parser.parse_args()
    
    input_dir = args.input
    threshold = args.threshold
    
    figures_dir = os.path.join(base_dir, "figures")
    results_dir = os.path.join(base_dir, "results")
    
    print("=== Commencing ColabFold Predicted Complex Post-Analysis ===")
    print(f"Reading from: {input_dir}")
    print(f"Using contact distance threshold: {threshold} Å\n")
    
    if not os.path.exists(input_dir):
        print(f"Error: Input directory not found: {input_dir}")
        sys.exit(1)
        
    pdb_file, json_file = get_pdb_and_json(input_dir)
    
    if not pdb_file:
        print("Error: No predicted PDB files found.")
        sys.exit(1)
    print(f"Found predicted PDB structure file: {os.path.basename(pdb_file)}")
    
    if not json_file:
        print("Warning: No predicted scores/PAE JSON file found. PAE and pLDDT plots will be skipped.")
    else:
        print(f"Found predicted scores JSON file: {os.path.basename(json_file)}")
        
    # 1. Parse dimer coordinates
    data = parse_dimer_structure(pdb_file)
    res_a = data['A']['res_info']
    res_b = data['B']['res_info']
    
    # 2. Analyze interface contact residues
    dist_matrix, contacts, interface_a, interface_b = analyze_interface(data, threshold)
    print(f"\n--- Interface Contact Summary (d <= {threshold} Å) ---")
    print(f"  Interface residues detected in Chain A (STN7): {len(interface_a)} / {len(res_a)} ({len(interface_a)/len(res_a)*100:.1f}%)")
    print(f"  Interface residues detected in Chain B (STN8): {len(interface_b)} / {len(res_b)} ({len(interface_b)/len(res_b)*100:.1f}%)")
    print(f"  Total inter-chain residue contacting pairs: {len(contacts)}")
    
    # 3. Parse JSON scores if available
    plddt, pae = None, None
    if json_file:
        try:
            plddt, pae = parse_scores_json(json_file)
            print("\nSuccessfully parsed pLDDT scores and PAE matrix.")
        except Exception as e:
            print(f"Warning: Failed to parse JSON scores: {str(e)}")
            
    # 4. Save detailed contact report
    report_path = os.path.join(results_dir, "stn7_stn8_interface_report.txt")
    with open(report_path, "w") as f:
        f.write("=================================================================\n")
        f.write("STN7-STN8 Predicted Heterodimer Interface Report\n")
        f.write(f"Analyzed structure: {os.path.basename(pdb_file)}\n")
        f.write(f"Contact threshold: {threshold} Angstroms\n")
        f.write("=================================================================\n\n")
        
        f.write(f"Chain A (STN7) length: {len(res_a)} residues\n")
        f.write(f"Chain B (STN8) length: {len(res_b)} residues\n\n")
        
        f.write(f"Interface residues in Chain A (STN7): {len(interface_a)}\n")
        # List them
        a_list = [f"{res_a[idx][2]}{res_a[idx][0]}" for idx in sorted(interface_a)]
        f.write(f"  Residues: {', '.join(a_list)}\n\n")
        
        f.write(f"Interface residues in Chain B (STN8): {len(interface_b)}\n")
        b_list = [f"{res_b[idx][2]}{res_b[idx][0]}" for idx in sorted(interface_b)]
        f.write(f"  Residues: {', '.join(b_list)}\n\n")
        
        f.write(f"Total atom contact pairs at interface: {len(contacts)}\n")
        f.write(f"{'Chain A (STN7)':<18} | {'Chain B (STN8)':<18} | {'Distance (Å)':<12}\n")
        f.write("-" * 56 + "\n")
        
        # Sort contact pairs by distance
        sorted_contacts = sorted(contacts, key=lambda x: x['distance'])
        for c in sorted_contacts:
            res_str_a = f"{c['one_letter_a']}{c['res_num_a']} ({c['res_name_a']})"
            res_str_b = f"{c['one_letter_b']}{c['res_num_b']} ({c['res_name_b']})"
            f.write(f"{res_str_a:<18} | {res_str_b:<18} | {c['distance']:.2f} Å\n")
            
    print(f"Saved interface contact report to: {report_path}")
    
    # 5. Plot Figures
    # A. Inter-Chain Contact Map
    plt.figure(figsize=(8, 7))
    contact_mask = (dist_matrix <= threshold).astype(float)
    plt.imshow(contact_mask, cmap='Blues', aspect='auto', origin='lower')
    plt.colorbar(label=f"Contact Mask (d <= {threshold} Å)")
    
    # Label mapping (using actual residue index numbers)
    tick_step_a = max(1, len(res_a) // 10)
    tick_step_b = max(1, len(res_b) // 10)
    
    plt.xticks(np.arange(0, len(res_a), tick_step_a), [res_a[i][0] for i in range(0, len(res_a), tick_step_a)], rotation=45)
    plt.yticks(np.arange(0, len(res_b), tick_step_b), [res_b[i][0] for i in range(0, len(res_b), tick_step_b)])
    
    plt.xlabel("Chain A: STN7 Residue Position", fontsize=11)
    plt.ylabel("Chain B: STN8 Residue Position", fontsize=11)
    plt.title("STN7-STN8 Predicted Dimer Interface Contact Map", fontsize=12, fontweight='bold', pad=10)
    plt.grid(True, linestyle=':', alpha=0.3)
    plt.tight_layout()
    
    contact_map_path = os.path.join(figures_dir, "stn7_stn8_interface_contacts.png")
    plt.savefig(contact_map_path, dpi=300)
    plt.close()
    print(f"Saved interface contact map plot to: {contact_map_path}")
    
    # B. PAE Matrix Plot
    if pae is not None:
        plt.figure(figsize=(9, 8))
        plt.imshow(pae, cmap='bwr', aspect='equal', origin='upper', vmin=0, vmax=31.75)
        plt.colorbar(label="Predicted Aligned Error (Å)")
        
        # Draw chain boundaries
        plt.axhline(len(res_a) - 0.5, color='black', linestyle='--', linewidth=1.5)
        plt.axvline(len(res_a) - 0.5, color='black', linestyle='--', linewidth=1.5)
        
        # Add labels for domains
        plt.text(len(res_a)/2, len(res_a)/2, "STN7\n(Intra-chain)", ha='center', va='center', fontsize=10, fontweight='bold')
        plt.text(len(res_a) + len(res_b)/2, len(res_a) + len(res_b)/2, "STN8\n(Intra-chain)", ha='center', va='center', fontsize=10, fontweight='bold')
        plt.text(len(res_a) + len(res_b)/2, len(res_a)/2, "Interface\n(Inter-chain)", ha='center', va='center', fontsize=10, fontweight='bold', color='darkred')
        plt.text(len(res_a)/2, len(res_a) + len(res_b)/2, "Interface\n(Inter-chain)", ha='center', va='center', fontsize=10, fontweight='bold', color='darkred')
        
        plt.xlabel("Aligned Residue Index", fontsize=11)
        plt.ylabel("Aligned Residue Index", fontsize=11)
        plt.title("STN7-STN8 Predicted Dimer PAE Matrix", fontsize=12, fontweight='bold', pad=10)
        plt.tight_layout()
        
        pae_path = os.path.join(figures_dir, "stn7_stn8_pae_matrix.png")
        plt.savefig(pae_path, dpi=300)
        plt.close()
        print(f"Saved PAE matrix heatmap to: {pae_path}")
        
    # C. Complex pLDDT Plot
    if plddt is not None:
        plt.figure(figsize=(10, 5))
        plt.plot(np.arange(len(res_a)), plddt[0:len(res_a)], color='#0ea5e9', label='STN7 (Chain A)', linewidth=1.5)
        plt.plot(np.arange(len(res_a), len(plddt)), plddt[len(res_a):], color='#f43f5e', label='STN8 (Chain B)', linewidth=1.5)
        
        plt.axvline(len(res_a) - 0.5, color='black', linestyle=':', linewidth=1.2)
        
        plt.axhline(y=70, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
        plt.text(0, 71, "High Confidence (pLDDT > 70)", color='gray', fontsize=8)
        plt.axhline(y=90, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
        plt.text(0, 91, "Very High Confidence (pLDDT > 90)", color='gray', fontsize=8)
        
        plt.xlabel("Residue Position in Dimer Complex", fontsize=11)
        plt.ylabel("pLDDT Score", fontsize=11)
        plt.ylim(0, 100)
        plt.legend(loc='lower left')
        plt.title("STN7-STN8 Predicted Dimer Complex pLDDT Confidence Profile", fontsize=12, fontweight='bold', pad=10)
        plt.grid(True, linestyle=':', alpha=0.4)
        plt.tight_layout()
        
        plddt_path = os.path.join(figures_dir, "stn7_stn8_complex_plddt.png")
        plt.savefig(plddt_path, dpi=300)
        plt.close()
        print(f"Saved complex pLDDT profile plot to: {plddt_path}")
        
    print("\n=== ColabFold Dimer Analysis Completed Successfully ===")
    print("=========================================================")

if __name__ == "__main__":
    main()
