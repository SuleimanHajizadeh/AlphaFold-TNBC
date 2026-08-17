---
title: "Source - STN7-STN8 Docking Project"
type: source
tags:
  - project-source/stn7-stn8
  - structural-biology/plant
  - thylakoid/kinases
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Source - STN7-STN8"
  - "STN7-STN8 Research Pipeline"
---

# Source - STN7-STN8 Docking Project

* **Path**: `projects/stn7-stn8-docking/`
* **Target Kinases**: Plant chloroplast thylakoid kinases [[STN7]] and [[STN8]]
* **Organisms Analyzed**: *Arabidopsis thaliana* (dicot), *Oryza sativa*, *Zea mays* (monocots), and *Chlamydomonas reinhardtii* (green algae outgroup)

---

## 📌 Executive Summary

A publication-quality comparative structural, evolutionary, and molecular docking bioinformatics pipeline analyzing thylakoid kinases **STN7** (mediator of [[Thylakoid State Transitions]]) and **STN8** (mediator of Photosystem II core repair). The workflow integrates:
1. Automated **AlphaFold v6** structure retrieval from EBI API.
2. Progressive star-alignment MSA and Shannon entropy analysis.
3. 3D structural superposition with Singular Value Decomposition (SVD).
4. Kyte-Doolittle hydrophobicity mapping.
5. Statistical validation ($t$-test on domain pLDDT, Pearson correlation).
6. **ColabFold-Multimer** heterodimer complex modeling and PAE interface mapping.
7. **AutoDock Vina** automated virtual screening pipeline (transit peptide cleavage, receptor preparation, grid-box auto-parameterization, and pose scoring).

---

## 📊 Key Quantitative Findings (from Manuscript Draft)

1. **Evolutionary & Paralog Conservation**:
   * *A. thaliana* STN7 vs STN8 sequence identity: **36.12%**.
   * Dicot vs monocot ortholog conservation: **73.78%** (Arabidopsis vs rice STN7) vs **35.66%** (Arabidopsis vs rice STN8).
   * Multi-species alignment (1089 columns): **57.02%** perfectly conserved columns ($H = 0.0$).
2. **Asymmetric Kinase Motifs**:
   * **G-loop**: `GEGSFG` (conserved across all homolog groups).
   * **Catalytic Loop**: `HRD` (Asp-279 in STN7, Asp-308 in STN8).
   * **Activation Loop**: Canonical **`DFG`** in STN8 vs non-canonical **`DLG`** in STN7 ($F \rightarrow L$ substitution).
   * **Substrate-Binding Loop**: Canonical **`APE`** in STN7 vs proline-substituted **`PPE`** in STN8.
3. **3D Structural Superposition ([[Structural Superposition & RMSD]])**:
   * Superimposed stromal catalytic core: **173 $C_\alpha$ pairs matched**.
   * Core $\text{RMSD} = \mathbf{4.2507\text{ Å}}$ (Average $C_\alpha$ distance = $3.67\text{ Å}$) using SVD.
4. **Domain Confidence Metrics ([[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]])**:
   * **STN7**: Overall mean pLDDT = 69.98 vs Kinase domain = **84.46** ($t = 29.32, p = 1.96 \times 10^{-102}$).
   * **STN8**: Overall mean pLDDT = 75.46 vs Kinase domain = **85.45** ($t = 17.21, p = 1.65 \times 10^{-41}$).
5. **Membrane Association & Hydrophobicity**:
   * Peak Kyte-Doolittle hydrophobicity: **1.49** (STN7, res 83-101) and **1.59** (STN8, res 101-119) ($w=19$), both below the $>1.6$ threshold for typical transmembrane spans.
   * Statistically significant profile correlation: Pearson $r = 0.2039, p = 0.0419$.
6. **ColabFold-Multimer Heterodimer Interface**:
   * **2,190 inter-chain contact pairs** ($d \leq 6\text{ Å}$), with 66.7% of STN7 and 78.0% of STN8 participating.
   * Low inter-chain PAE ($< 10\text{ Å}$) across STN7 residues 100-130 and STN8 residues 120-150.
7. **Automated AutoDock Vina Virtual Screening**:
   * Automated transit peptide sequence truncation via `prepare_receptor.py`.
   * Dynamic box calculation, multi-ligand screening with Vina 1.2.5, and ensemble pose scoring.

---

## 🔗 Referenced Wiki Entities & Concepts
* Entities: [[STN7]], [[STN8]], [[AlphaFold2]].
* Concepts: [[Thylakoid State Transitions]], [[Structural Superposition & RMSD]], [[Local Distance Difference Test (LDDT & pLDDT)]], [[Kinase Activation Loop & Allostery]].
* Synthesis: [[Comparative Kinase Mechanics - AKT1 vs STN7-STN8]].
