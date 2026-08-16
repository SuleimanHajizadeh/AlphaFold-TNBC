---
title: "Structural Superposition & RMSD"
type: concept
tags:
  - structural-biology/alignment
  - superposition
  - rmsd
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "RMSD"
  - "Root-Mean-Square Deviation"
  - "Structural Alignment"
  - "Kabsch Algorithm"
---

# Structural Superposition & RMSD

**Root-Mean-Square Deviation (RMSD)** is the fundamental quantitative measure of the average distance between corresponding backbones (typically $C_\alpha$ atoms) of two superimposed protein 3D structures.

---

## 📐 Mathematical Formulation

Given two structures with $N$ paired atoms having Cartesian coordinates $\mathbf{x}_i$ and $\mathbf{y}_i$ after optimal rigid-body superposition (rotation matrix $\mathbf{R}$ and translation vector $\mathbf{t}$):

$$\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^N \| \mathbf{x}_i - (\mathbf{R} \mathbf{y}_i + \mathbf{t}) \|^2}$$

### The Kabsch Algorithm (Singular Value Decomposition - SVD)
1. Center both coordinate sets at their respective centers of mass.
2. Compute the cross-covariance matrix $\mathbf{C} = \mathbf{X}^T \mathbf{Y}$.
3. Decompose $\mathbf{C}$ via SVD: $\mathbf{C} = \mathbf{V} \mathbf{S} \mathbf{W}^T$.
4. The optimal rotation matrix is $\mathbf{R} = \mathbf{V} \mathbf{D} \mathbf{W}^T$, where $\mathbf{D} = \operatorname{diag}(1, 1, \det(\mathbf{V}\mathbf{W}^T))$ ensures a right-handed coordinate rotation without reflections.

---

## 🔬 Domain-Level Local Alignment vs Global Alignment

* **Global RMSD**: Sensitive to flexible loop excursions and rigid domain-domain hinge shifts. Two structurally identical kinase lobes can produce a high global RMSD if connected by a flexible linker.
* **Core Kinase Superposition**: Local domain alignment of catalytic cores (e.g., between [[STN7]] and [[STN8]]) isolates the catalytic triad and ATP-binding loop, yielding high-precision matching ($\text{RMSD} = 4.25\text{ Å}$ across 173 matched core atoms).

---

## 🔗 Key Cross-References
* **Superposition-Free Metric**: [[Local Distance Difference Test (LDDT & pLDDT)|LDDT & pLDDT]].
* **Kinase Paralog Analysis**: [[Source - STN7-STN8 Docking Project]], [[STN7]], [[STN8]].
* **Comparative Syntheses**: [[Comparative Kinase Mechanics - AKT1 vs STN7-STN8]].
