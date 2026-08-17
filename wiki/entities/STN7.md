---
title: "STN7 (State Transition Kinase 7)"
type: entity
tags:
  - kinase/plant
  - structural-biology/photosynthesis
  - chloroplast/thylakoid
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "STN7"
  - "State Transition Kinase 7"
  - "At1g68830"
  - "UniProt:Q949Z3"
---

# STN7 (State Transition Kinase 7)

**STN7** is a chloroplast thylakoid membrane-associated serine/threonine protein kinase found across green algae and land plants. It plays an indispensable role in **photosynthetic acclimation** by catalyzing the reversible phosphorylation of the Light-Harvesting Complex II (LHCII), driving [[Thylakoid State Transitions]].

---

## 🌿 Biological Role & Function

```
Plastoquinone Pool Reduced (High PSII excitation)
                 │
                 ▼
          STN7 Activation
                 │
                 ▼
       Phosphorylation of LHCII
                 │
                 ▼
LHCII detaches from PSII and migrates to PSI (State 2 Transition)
```

* **Target Substrate**: Mobile trimeric LHCII proteins (Lhcb1, Lhcb2).
* **Physiological Role**: Under preferential excitation of Photosystem II (PSII), the plastoquinone (PQ) pool becomes over-reduced. STN7 senses the redox state via interaction with the cytochrome $b_6f$ complex and phosphorylates LHCII. Phosphorylated LHCII dissociates from PSII and associates with PSI (State 2), rebalancing excitation energy distribution.
* **Regulation**: Rapidly inactivated under extreme high-light stress or oxidizing conditions to protect photosynthetic machinery.

---

## 🔬 Structural & Computational Characteristics

* **Evolutionary Conservation**: Orthologs are conserved across green algae (*Chlamydomonas reinhardtii* STT7: `Q84V18`), monocots (*Oryza sativa* `B9FLG7`, *Zea mays* `A0A3L6E9I0`), and dicots (*Arabidopsis thaliana* `Q9S713`).
  * High dicot-monocot ortholog conservation: **73.78%** identity between *A. thaliana* and rice.
* **Paralog Relationship with [[STN8]]**: Displays **36.12%** pairwise sequence identity with *A. thaliana* STN8, arising from an ancient gene duplication event.
* **Conserved & Divergent Kinase Motifs**:
  * **ATP-Binding G-Loop (P-loop)**: Canonical `GEGSFG` (residues 141-146 in Arabidopsis).
  * **Catalytic Base Loop (HRD Motif)**: Perfectly conserved `HRD` (residues 277-279, with Asp-279 acting as the proton acceptor).
  * **Non-Canonical Activation Loop**: Features a non-canonical **`DLG`** motif (residues 295-297). The key $\text{Phe} \rightarrow \text{Leu}$ ($F \rightarrow L$) substitution influences metal coordination and is an evolutionary signature of STN7.
  * **Substrate-Binding APE Motif**: Canonical `APE` motif (residues 324-326).
* **Domain Confidence & Secondary Structure (AlphaFold v6)**:
  * Overall Mean [[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]]: **69.98** vs Kinase Domain Mean pLDDT: **84.46** ($t = 29.32, p = 1.96 \times 10^{-102}$, Student's t-test).
  * Secondary structure (Kinase domain): **58.6%** $\alpha$-helix, **31.3%** $\beta$-strand, **10.0%** coil.
* **Membrane Anchor**: Lacks a canonical transmembrane helix (maximum Kyte-Doolittle hydrophobicity = **1.49** at residues 83-101, $w=19$), exhibiting an atypical peripheral/loop membrane association mode.
* **ColabFold Heterodimer Complex with [[STN8]]**:
  * Forms a putative heterodimer interface with 2,190 inter-chain contacts ($d \leq 6\text{ Å}$), with 66.7% (375/562) of STN7 residues participating in the contact zone.
  * Displays low inter-chain PAE ($< 10\text{ Å}$) across residues 100–130.

---

## 🔗 Key Cross-References
* **Paralog Kinase**: [[STN8]].
* **Biological Mechanism**: [[Thylakoid State Transitions]].
* **Comparative Syntheses**: [[Comparative Kinase Mechanics - AKT1 vs STN7-STN8]].
* **Methodologies**: [[Structural Superposition & RMSD]], [[AlphaFold2]].
* **Primary Project Source**: [[Source - STN7-STN8 Docking Project]].
