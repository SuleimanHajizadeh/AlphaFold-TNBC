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

A publication-quality comparative structural and evolutionary bioinformatics pipeline analyzing thylakoid kinases **STN7** (mediator of [[Thylakoid State Transitions]]) and **STN8** (mediator of Photosystem II core repair). The workflow integrates automated [[AlphaFold2]] database retrieval, progressive star-alignment MSA, 3D structural superposition with SVD, Kyte-Doolittle hydrophobicity mapping, and statistical validation.

---

## 📊 Key Quantitative Findings

1. **Evolutionary & Paralog Conservation**:
   * *A. thaliana* STN7 vs STN8 sequence identity: **36.12%**.
   * Multi-species alignment: **57.02%** of alignment columns (621 / 1089) exhibit perfect conservation (Shannon Entropy $H = 0.0$).
2. **3D Structural Superposition ([[Structural Superposition & RMSD]])**:
   * Superimposed stromal catalytic core: **173 $C_\alpha$ pairs matched**.
   * Core $\text{RMSD} = \mathbf{4.2507\text{ Å}}$ using Singular Value Decomposition (SVD).
3. **Domain Confidence Metrics ([[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]])**:
   * Kinase catalytic domain: High definition (Mean pLDDT $\approx \mathbf{85}$).
   * Significantly higher confidence in catalytic core vs variable transit peptides/loops ($p < 0.001$, Student's t-test).
4. **Membrane Association & Hydrophobicity**:
   * Peak Kyte-Doolittle hydrophobicity: **1.49** (STN7) and **1.59** (STN8) ($w=19$), both below the $>1.6$ threshold for typical transmembrane spans.
   * Moderate but statistically significant profile correlation (Pearson $r = 0.2039, p < 0.05$).

---

## 🔗 Referenced Wiki Entities & Concepts
* Entities: [[STN7]], [[STN8]], [[AlphaFold2]].
* Concepts: [[Thylakoid State Transitions]], [[Structural Superposition & RMSD]], [[Local Distance Difference Test (LDDT & pLDDT)]], [[Kinase Activation Loop & Allostery]].
* Synthesis: [[Comparative Kinase Mechanics - AKT1 vs STN7-STN8]].
