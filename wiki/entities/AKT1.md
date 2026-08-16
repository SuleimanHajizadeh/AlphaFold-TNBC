---
title: "AKT1 (RAC-alpha serine/threonine-protein kinase)"
type: entity
tags:
  - kinase/human
  - structural-biology/cancer
  - oncology/tnbc
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "AKT1"
  - "RAC-alpha serine/threonine-protein kinase"
  - "PKB alpha"
  - "UniProt:P31749"
---

# AKT1 (RAC-alpha serine/threonine-protein kinase)

**AKT1** (UniProt: [P31749](https://www.uniprot.org/uniprot/P31749)) is a central human serine/threonine kinase that functions as the key nodal hub in the PI3K/Akt/mTOR signaling cascade. It regulates cell survival, proliferation, metabolic reprogramming, and migration.

Hyperactivation and dysregulation of AKT1 are hallmark drivers in Triple-Negative Breast Cancer (TNBC), where transcriptomic co-expression networks (WGCNA) identify AKT1 as a top master hub kinase.

---

## 🧬 Structural Architecture & Domains

```
N-term ─── [ PH Domain ] ─── [ Flexible Linker ] ─── [ Kinase Catalytic Domain ] ─── [ Regulatory C-tail (HM) ] ─── C-term
             (Membrane)           (Allostery)            (ATP / Substrate)               (Hydrophobic Motif)
```

1. **Pleckstrin Homology (PH) Domain**: Binds phosphatidylinositol-(3,4,5)-trisphosphate ($\text{PIP}_3$) to anchor the kinase to the plasma membrane.
2. **Inter-domain Linker**: A flexible, intrinsically disordered linker region connecting the PH domain to the catalytic domain.
3. **Kinase Catalytic Domain**: Highly conserved core containing the ATP-binding pocket, catalytic loop (HRD motif), and the [[Kinase Activation Loop & Allostery|Activation Loop]] containing Thr308.
4. **C-Terminal Regulatory Tail**: Contains the Hydrophobic Motif (HM) with Ser473, essential for maximal kinase activation when phosphorylated by mTORC2.

---

## 🔬 AlphaFold2 Structural Prediction Profile

Based on [[Source - AKT1 Kinase Modeling Project|ColabFold v1.6.1 Modeling]] of the 291-residue catalytic segment:

* **Mean Confidence**: Mean [[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]] = **65.61**, pTM-score = **0.450**.
* **Catalytic Core Resolution**: The catalytic core shows very high confidence (pLDDT $\geq 90$, 24.7% of residues; confident $70\text{–}90$, 23.7%).
* **Intrinsically Disordered Regions (IDRs)**: 30.2% of residues exhibit pLDDT $< 50$, corresponding to the regulatory flexible loop and PH linker. This high plasticity mediates allosteric conformational shifts upon membrane recruitment.
* **Secondary Structure Topology**:
  * $\alpha$-helix: 24.6%
  * $\beta$-strand: 32.9%
  * Loops & flexible turns: 39.4%
  * Left-handed / Glycine: 3.1%
* **Tertiary Contacts**: 573 active long-range $C_\alpha$ pairwise contacts ($D_{ij} \leq 8.0\text{ Å}, |i-j| \geq 6$), spanning Euclidean distances up to $76.69\text{ Å}$.

---

## 🔗 Key Cross-References
* **Modeling Methodology**: [[AlphaFold2]], [[Ramachandran Dihedral Angles]], [[Contact Map Analysis]].
* **Mechanisms**: [[Kinase Activation Loop & Allostery]].
* **Comparative Syntheses**: [[Comparative Kinase Mechanics - AKT1 vs STN7-STN8]].
* **Primary Project Source**: [[Source - AKT1 Kinase Modeling Project]].
