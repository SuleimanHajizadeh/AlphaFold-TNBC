---
title: "Local Distance Difference Test (LDDT & pLDDT)"
type: concept
tags:
  - structural-biology/metrics
  - protein-structure-prediction/validation
  - alphafold
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "LDDT"
  - "pLDDT"
  - "Local Distance Difference Test"
---

# Local Distance Difference Test (LDDT & pLDDT)

The **Local Distance Difference Test (LDDT)** (Mariani et al., 2013) is a superposition-free structural metric used to assess local model accuracy by comparing internal interatomic distance matrices rather than performing global structural alignments (such as RMSD).

**pLDDT** is the machine learning prediction of this score produced directly by deep learning systems like [[AlphaFold2]].

---

## 📐 Mathematical Formulation

Let $D$ represent all pairs of $C_\alpha$ atoms $(i, j)$ in the ground-truth structure separated by Euclidean distance $\leq R$ (inclusion radius $R = 15\text{ Å}$):

$$D = \{ (i,j) \mid d_{\text{true}}(i,j) \leq R, \ i \neq j \}$$

The distance deviation $\Delta_{ij}$ between predicted structure ($d_{\text{pred}}$) and true structure ($d_{\text{true}}$) is:

$$\Delta_{ij} = d_{\text{pred}}(i,j) - d_{\text{true}}(i,j)$$

The LDDT score calculates the fraction of preserved pairwise distances across four strict tolerance thresholds ($0.5\text{ Å}, 1.0\text{ Å}, 2.0\text{ Å}, 4.0\text{ Å}$):

$$\text{LDDT} = \frac{1}{4 |D|} \sum_{(i,j) \in D} \left[ \mathbb{I}(|\Delta_{ij}| \leq 0.5\text{Å}) + \mathbb{I}(|\Delta_{ij}| \leq 1.0\text{Å}) + \mathbb{I}(|\Delta_{ij}| \leq 2.0\text{Å}) + \mathbb{I}(|\Delta_{ij}| \leq 4.0\text{Å}) \right]$$

---

## 🎯 Confidence Bands & Interpretation

| pLDDT Interval | Color Code | Structural Interpretation |
|---|---|---|
| **$\geq 90$** | 🔵 Dark Blue | Very high confidence; accurate side-chain rotamers. |
| **$70\text{–}90$** | 🩵 Light Blue | Confident; accurate backbone topology. |
| **$50\text{–}70$** | 🟡 Yellow | Low confidence; tentative fold / flexible loop. |
| **$< 50$** | 🟠 Orange | Very low confidence; **Intrinsically Disordered Regions (IDRs)**. |

---

## 🔗 Key Cross-References
* **Modeling Engine**: [[AlphaFold2]].
* **Applications**: [[AKT1]] confidence profiling, [[STN7]] / [[STN8]] stromal vs loop validation.
* **Complementary Metrics**: [[Ramachandran Dihedral Angles]], [[Contact Map Analysis]], [[Structural Superposition & RMSD]].
