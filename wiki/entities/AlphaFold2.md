---
title: "AlphaFold2"
type: entity
tags:
  - machine-learning/structural-biology
  - protein-structure-prediction
  - deepmind
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "AlphaFold2"
  - "AF2"
  - "ColabFold"
---

# AlphaFold2

**AlphaFold2** is a deep learning system developed by Google DeepMind (Jumper et al., 2021) that predicts 3D coordinates of proteins from their primary amino acid sequences with atomic accuracy.

In this repository, AlphaFold2 is deployed via **ColabFold v1.6.1** (using MMseqs2 fast homology searching) and the **AlphaFold Protein Structure Database (AlphaFold DB v4/v6)**.

---

## 🏗️ Model Architecture

```
Sequence / MSA (MMseqs2) ───► [ Evoformer: 48 blocks ] ───► [ Structure Module: 8 blocks ] ───► 3D Backbone + Sidechains
                                 (Pairwise + MSA reps)          (Invariant Point Attention - IPA)
```

1. **Multiple Sequence Alignment (MSA) & Pair Representation**: Captures co-evolutionary patterns between residue pairs.
2. **Evoformer Blocks**: Uses axial attention across sequence alignments and triangle updates across spatial pair representations.
3. **Structure Module**: Utilizes Invariant Point Attention (IPA) and 3D rigid-body transformations without bond geometry constraints to directly predict $C_\alpha$ positions and residue rotations ($SO(3)$).

---

## 📊 Core Confidence & Error Metrics

AlphaFold2 generates calibrated internal metrics:

### 1. Predicted Local Distance Difference Test (pLDDT)
* Per-residue confidence metric ($0\text{–}100$).
* Direct prediction of physical [[Local Distance Difference Test (LDDT & pLDDT)|LDDT]].
* **pLDDT $\geq 90$**: High-accuracy backbone and side chains.
* **pLDDT $< 50$**: Intrinsically Disordered Regions (IDRs).

### 2. Predicted Aligned Error (PAE)
* Evaluates relative spatial orientation error (in Å) between any two residues $i$ and $j$.
* Key for evaluating domain-domain packing versus flexible inter-domain linkers.

### 3. Predicted TM-Score (pTM)
* Measures global structural fold confidence on a scale of $0.0\text{–}1.0$.

---

## 🔗 Key Cross-References
* **Evaluation Metrics**: [[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]], [[Ramachandran Dihedral Angles]], [[Contact Map Analysis]].
* **Projects Modeled with AF2**: [[Source - AKT1 Kinase Modeling Project]], [[Source - STN7-STN8 Docking Project]].
* **Entities Analyzed**: [[AKT1]], [[STN7]], [[STN8]].
