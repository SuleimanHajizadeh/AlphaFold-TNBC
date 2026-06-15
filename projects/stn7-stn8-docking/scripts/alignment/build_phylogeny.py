#!/usr/bin/env python3
"""
STN7 and STN8 Evolutionary Phylogeny Construction Script
Reads the combined multiple sequence alignment (all_msa.fasta),
calculates the pairwise amino acid distance matrix, constructs a
Neighbor-Joining (NJ) phylogenetic tree, saves it in standard Newick format,
and generates a publication-quality phylogenetic tree figure.
"""

import os
import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    msa_path = os.path.join(base_dir, "data/msa/all_msa.fasta")
    
    figures_dir = os.path.join(base_dir, "figures")
    results_dir = os.path.join(base_dir, "results")
    phylogeny_dir = os.path.join(base_dir, "data/phylogeny")
    
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(phylogeny_dir, exist_ok=True)
    
    print("=== Commencing Evolutionary Phylogenetic Tree Construction ===")
    
    # 1. Read Multiple Sequence Alignment
    print("Reading MSA alignment...")
    alignment = AlignIO.read(msa_path, "fasta")
    
    # Clean sequence IDs for neat tree labels
    for record in alignment:
        record.id = record.id.split()[0].replace("STN7_", "STN7:").replace("STN8_", "STN8:")
        record.name = record.id
        
    # 2. Compute Distance Matrix
    print("Calculating pairwise sequence distance matrix (blosum62 model)...")
    calculator = DistanceCalculator('blosum62')
    dm = calculator.get_distance(alignment)
    
    # Save distance matrix report
    dm_txt_path = os.path.join(results_dir, "phylogeny_distance_matrix.txt")
    with open(dm_txt_path, "w") as f:
        f.write(str(dm))
    print(f"Saved distance matrix to: {dm_txt_path}")
    
    # 3. Construct Neighbor-Joining Tree
    print("Constructing Neighbor-Joining (NJ) phylogenetic tree...")
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(dm)
    
    # Set tree root if needed (optional, NJ is unrooted, but let's root it to alga outgroup for publication style!)
    # Chlamydomonas reinhardtii is our evolutionary outgroup
    outgroup_id = None
    for term in nj_tree.get_terminals():
        if "Chlamydomonas" in term.name:
            outgroup_id = term.name
            break
            
    if outgroup_id:
        print(f"Rooting the tree with outgroup: {outgroup_id}")
        nj_tree.root_with_outgroup(outgroup_id)
        
    # Standardize branch lengths and names
    nj_tree.ladderize() # Sort branches for neat presentation
    
    # 4. Save Tree in Newick Format
    tree_path = os.path.join(phylogeny_dir, "kinase_nj_tree.nwk")
    Phylo.write(nj_tree, tree_path, "newick")
    print(f"Saved phylogenetic tree in Newick format to: {tree_path}")
    
    # 5. Visualize Tree in publication quality
    plt.figure(figsize=(10, 6))
    # Elegant styling for trees in Matplotlib
    # We use a custom color for tree branches and labels
    fig = plt.gcf()
    ax = fig.add_subplot(1, 1, 1)
    
    # Customize the drawing of branches
    Phylo.draw(
        nj_tree, 
        do_show=False, 
        axes=ax, 
        branch_labels=None,
        label_func=lambda x: str(x).replace("Inner", "") if x.name else ""
    )
    
    # Decorate plot
    plt.title("Neighbor-Joining Phylogenetic Tree of STN7 and STN8 Chloroplast Kinases", fontsize=12, fontweight='bold', pad=15)
    plt.xlabel("Evolutionary Distance", fontsize=10)
    plt.ylabel("Protein Lineages", fontsize=10)
    
    # Stylistic tweaks
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_visible(False) # Tree draw handles y-axis labels
    ax.grid(True, linestyle=':', alpha=0.4)
    
    plt.tight_layout()
    plot_path = os.path.join(figures_dir, "phylogenetic_tree.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    
    print(f"Saved phylogenetic tree plot to: {plot_path}")
    print("==========================================================")

if __name__ == "__main__":
    main()
