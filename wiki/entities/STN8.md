---
title: "STN8 (State Transition Kinase 8)"
type: entity
tags:
  - kinase/plant
  - structural-biology/photosynthesis
  - chloroplast/photoprotection
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "STN8"
  - "State Transition Kinase 8"
  - "At5g01920"
  - "UniProt:Q9FM15"
---

# STN8 (State Transition Kinase 8)

**STN8** is a chloroplast thylakoid protein kinase in photosynthetic organisms dedicated to the phosphorylation of **Photosystem II (PSII) core reaction center subunits** (D1, D2, CP43, and PsbH) under high-light and photoinhibitory stress.

---

## ☀️ Biological Function & Photoprotection

Unlike its paralog [[STN7]] which modulates light harvesting antennae via [[Thylakoid State Transitions]], **STN8** is primarily responsible for the **PSII repair and turnover cycle**:

```
High Light & Photo-damage
           │
           ▼
    STN8 Activation
           │
           ▼
Phosphorylation of PSII Core (D1 / D2 / CP43)
           │
           ▼
Unstacking of Grana & Facilitated FtsH Protease Degradation of Damaged D1
```

* **Target Substrates**: PSII reaction center proteins (D1/PsbA, D2/PsbD, CP43/PsbC, and PsbH).
* **Physiological Role**: Phosphorylation of damaged PSII reaction center proteins promotes the disassembly of PSII-LHCII supercomplexes and allows damaged D1 to be degraded by the stromal FtsH protease and replaced by newly synthesized D1.

---

## 🔬 Structural Properties & Comparative Metrics

* **Paralog Comparison with [[STN7]]**:
  * Sequence Identity: **36.12%** (with *A. thaliana* STN7).
  * Structural Kinase Core Superposition: **173 $C_\alpha$ matched atoms**, **$\text{RMSD} = 4.2507\text{ Å}$** ([[Structural Superposition & RMSD]]).
  * Active Site Conservation: Catalytic residues and ATP-coordinating triad are strictly preserved.
* **Membrane Hydrophobicity Profile**:
  * Peak Kyte-Doolittle Hydrophobicity score = **1.59** (below the typical >1.6 cutoff for standard transmembrane helices).
  * Significantly correlated hydrophobicity profile with STN7 (Pearson $r = 0.2039, p < 0.05$).
* **Structural Model Confidence**:
  * High stromal kinase domain confidence (mean [[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]] $\approx 85$) contrasted with disordered transit peptide and loop regions ($p < 0.001$, Student's t-test).

---

## 🔗 Key Cross-References
* **Paralog Kinase**: [[STN7]].
* **Biological Pathway**: [[Thylakoid State Transitions]].
* **Comparative Syntheses**: [[Comparative Kinase Mechanics - AKT1 vs STN7-STN8]].
* **Primary Project Source**: [[Source - STN7-STN8 Docking Project]].
