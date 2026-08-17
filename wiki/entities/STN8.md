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
  * Sequence Identity: **36.12%** (with *A. thaliana* STN7). Higher divergence in monocot STN8 (rice vs Arabidopsis identity is **35.66%**).
  * Structural Kinase Core Superposition: **173 $C_\alpha$ matched atoms**, **$\text{RMSD} = 4.2507\text{ Å}$** ([[Structural Superposition & RMSD]]), average $C_\alpha$ distance = $3.67\text{ Å}$.
* **Conserved & Divergent Kinase Motifs**:
  * **ATP-Binding G-Loop (P-loop)**: Conserved `GEGSFG` (residues 140-145 in Arabidopsis).
  * **Catalytic Base Loop (HRD Motif)**: Perfectly conserved `HRD` (residues 306-308, with Asp-308 acting as the catalytic proton acceptor).
  * **Canonical Activation Loop**: Maintains the standard **`DFG`** motif (residues 326-328 in Arabidopsis) for coordinating catalytic $\text{Mg}^{2+}$ ions, contrasting with STN7's non-canonical `DLG`.
  * **Divergent Substrate-Binding Loop**: Contains a non-canonical **`PPE`** motif (residues 346-348 in Arabidopsis), representing a proline-substituted variant of the standard `APE` motif.
* **Structural Model Confidence & Secondary Structure (AlphaFold v6)**:
  * Overall Mean [[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]]: **75.46** vs Kinase Domain Mean pLDDT: **85.45** ($t = 17.21, p = 1.65 \times 10^{-41}$, Student's t-test).
  * Secondary structure (Kinase domain): **53.0%** $\alpha$-helix, **34.5%** $\beta$-strand, **12.5%** coil.
* **Membrane Hydrophobicity Profile**:
  * Peak Kyte-Doolittle Hydrophobicity score = **1.59** (residues 101-119, $w=19$, below the $>1.6$ threshold for standard transmembrane helices).
  * Significant correlation with STN7 hydrophobicity profile (Pearson $r = 0.2039, p = 0.0419$).
* **ColabFold Heterodimer Complex with [[STN7]]**:
  * Forms 2,190 inter-chain contacts ($d \leq 6\text{ Å}$), with 78.0% (386/495) of STN8 participating in the interface.
  * Inter-chain PAE is low ($< 10\text{ Å}$) across residues 120–150.

---

## 🔗 Key Cross-References
* **Paralog Kinase**: [[STN7]].
* **Biological Pathway**: [[Thylakoid State Transitions]].
* **Comparative Syntheses**: [[Comparative Kinase Mechanics - AKT1 vs STN7-STN8]].
* **Primary Project Source**: [[Source - STN7-STN8 Docking Project]].
