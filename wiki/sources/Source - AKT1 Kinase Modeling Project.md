---
title: "Source - AKT1 Kinase Modeling Project"
type: source
tags:
  - project-source/akt1
  - structural-biology/af2
  - oncology/tnbc
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Source - AKT1"
  - "AKT1 Modeling Project"
---

# Source - AKT1 Kinase Modeling Project

* **Path**: `projects/akt1-kinase-modeling/`
* **Target Protein**: Human AKT1 (*RAC-alpha serine/threonine-protein kinase*, UniProt: [P31749](https://www.uniprot.org/uniprot/P31749))
* **Pipeline Engine**: [[AlphaFold2]] (ColabFold v1.6.1 MMseqs2 MSA)

---

## 📌 Executive Summary

This computational structural biology project models the 3D structure of the central Triple-Negative Breast Cancer (TNBC) hub kinase **AKT1**. The pipeline extracts and evaluates prediction confidence ([[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]]), backbone dihedral torsion vectors ([[Ramachandran Dihedral Angles]]), $C_\alpha$ pairwise Euclidean distance matrices ([[Contact Map Analysis]]), and Predicted Aligned Error (PAE) matrices.

---

## 📊 Key Quantitative Findings

1. **Confidence Profile**:
   * Mean pLDDT: **65.61** | Predicted TM-score (pTM): **0.450**
   * Very High Confidence ($\text{pLDDT} \geq 90$): **24.7%** (72 residues)
   * Confident ($70\text{–}90$): **23.7%** (69 residues)
   * Low ($50\text{–}70$): **21.3%** (62 residues)
   * Very Low / Disordered ($\text{pLDDT} < 50$): **30.2%** (88 residues)
2. **Backbone Dihedral Conformations**:
   * Favored $\alpha$-helix: **24.6%**
   * Favored $\beta$-strand: **32.9%**
   * Loops / other allowed: **39.4%**
   * Glycine / Left-handed: **3.1%**
3. **Tertiary Distance Geometry**:
   * Mapped Residues: 291
   * Active long-range contacts ($D_{ij} \leq 8.0\text{ Å}, |i-j| \geq 6$): **573** (1.4% contact density)
   * Distance Range: $3.02\text{ Å} \text{ to } 76.69\text{ Å}$

---

## 🔬 Biological Insight
The structural modeling captures the intrinsically disordered Pleckstrin Homology (PH) domain linker and C-terminal regulatory tail (pLDDT < 50), which provide the physical plasticity required for membrane recruitment and [[Kinase Activation Loop & Allostery|allosteric activation]] in oncogenic environments.

---

## 🔗 Referenced Wiki Entities & Concepts
* Entities: [[AKT1]], [[AlphaFold2]].
* Concepts: [[Local Distance Difference Test (LDDT & pLDDT)]], [[Ramachandran Dihedral Angles]], [[Contact Map Analysis]], [[Kinase Activation Loop & Allostery]].
* Synthesis: [[Comparative Kinase Mechanics - AKT1 vs STN7-STN8]].
