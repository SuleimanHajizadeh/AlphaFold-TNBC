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

* **Evolutionary Conservation**: Orthologs are conserved across green algae (*Chlamydomonas reinhardtii*), monocots (*Oryza sativa*, *Zea mays*), and dicots (*Arabidopsis thaliana*).
* **Paralog Relationship with [[STN8]]**: Displays **36.12%** pairwise sequence identity with paralog STN8 in *A. thaliana*, arising from an ancient gene duplication event.
* **Catalytic Domain**: Contains a highly structured stromal catalytic core (mean [[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]] $\approx 85$) superimposed with STN8 across 173 $C_\alpha$ positions with an RMSD of **4.25 Å**.
* **Membrane Anchor**: Lacks a classical deep transmembrane $\alpha$-helix (maximum Kyte-Doolittle hydrophobicity score = **1.49**, window $w=19$), exhibiting an atypical peripheral/loop membrane association mode.

---

## 🔗 Key Cross-References
* **Paralog Kinase**: [[STN8]].
* **Biological Mechanism**: [[Thylakoid State Transitions]].
* **Comparative Syntheses**: [[Comparative Kinase Mechanics - AKT1 vs STN7-STN8]].
* **Methodologies**: [[Structural Superposition & RMSD]], [[AlphaFold2]].
* **Primary Project Source**: [[Source - STN7-STN8 Docking Project]].
