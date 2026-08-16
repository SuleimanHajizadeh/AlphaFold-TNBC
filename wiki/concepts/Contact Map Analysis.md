---
title: "Contact Map Analysis"
type: concept
tags:
  - structural-biology/topology
  - distance-matrix
  - contact-map
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Contact Map"
  - "Distance Matrix"
  - "Residue Pairwise Topology"
---

# Contact Map Analysis

A **Contact Map** is a 2D binary or distance-weighted symmetric matrix representation of a protein's 3D tertiary structure. It captures the global folding topology by mapping spatial proximities between all pairs of amino acid residues in the chain.

---

## 📐 Mathematical Formulation

Given $N$ residues with 3D Cartesian coordinates of their $C_\alpha$ atoms $(x_i, y_i, z_i)$, the $N \times N$ Euclidean distance matrix $D$ is populated by:

$$D_{ij} = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2}$$

### Contact Declaration Threshold
A true tertiary contact between residue $i$ and residue $j$ is declared when:

$$D_{ij} \leq d_{\text{cutoff}} \quad \text{and} \quad |i - j| \geq k_{\text{sep}}$$

* Standard distance cutoff: $d_{\text{cutoff}} = 8.0\text{ Å}$ (or $6.0\text{ Å}$ for heavy atoms).
* Sequence-separation filter: $k_{\text{sep}} = 6$ (eliminates trivial local backbone contacts and focuses on long-range folding architecture).

---

## 📊 Characteristic Visual Patterns

* **Main Diagonal**: Local sequence continuity ($i \approx j$).
* **Thick Parallel Off-diagonals**: $\alpha$-helices (contacts at $i, i+3$ and $i, i+4$).
* **Orthogonal Cross-lines**: Antiparallel $\beta$-sheets.
* **Parallel Off-diagonal Bands**: Parallel $\beta$-sheets.
* **Dispersed Points**: Tertiary contacts between distant loops and domains.

---

## 🔗 Key Cross-References
* **Application**: [[Source - AKT1 Kinase Modeling Project]] (573 long-range contacts mapped across [[AKT1]]).
* **Related Concepts**: [[Ramachandran Dihedral Angles]], [[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]].
