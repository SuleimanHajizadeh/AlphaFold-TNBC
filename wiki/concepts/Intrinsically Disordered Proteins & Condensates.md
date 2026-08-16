---
title: "Intrinsically Disordered Proteins & Condensates"
type: concept
tags:
  - structural-biology/disorder
  - idr/conformational-ensembles
  - phase-separation/condensates
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Intrinsically Disordered Proteins"
  - "IDP"
  - "IDR"
  - "Liquid-Liquid Phase Separation"
  - "LLPS"
---

# Intrinsically Disordered Proteins & Condensates

**Intrinsically Disordered Proteins (IDPs)** and **Intrinsically Disordered Regions (IDRs)** are functional polypeptide segments that lack a stable, unique 3D equilibrium structure under physiological conditions. Instead, they exist as dynamic, interconverting conformational ensembles.

---

## ⚡ Sequence Signatures & Biophysical Properties

```
Structured Fold (Hydrophobic Core)          Intrinsically Disordered Segment (Extended Ensemble)
  * High hydrophobic content (Leu, Ile, Val)   * Low mean hydrophobicity
  * High aromatic residues (Phe, Tyr, Trp)     * High net charge (Lys, Arg, Glu, Asp)
  * Fixed Ramachandran coordinates             * High flexibility promoters (Gly, Pro, Ser, Gln)
```

1. **AlphaFold2 Confidence Correlation**:
   * Regions with **pLDDT $< 50$** correlate strongly with true physical IDRs (demonstrated by NMR and SAXS experiments).
   * In [[AKT1]], 30.2% of residues show pLDDT $< 50$, corresponding to the regulatory flexible PH-linker and C-tail.
2. **Coupled Folding Upon Binding**:
   * Many IDRs undergo disorder-to-order transitions upon binding their physiological interaction partner (e.g., [[TP53]] transactivation domain binding MDM2).
3. **Biomolecular Condensates & Phase Separation (LLPS)**:
   * Multivalent, low-affinity interactions mediated by IDRs drive liquid-liquid phase separation (LLPS), forming membraneless organelles (nucleoli, stress granules, transcriptional hubs).

---

## 🔗 Key Cross-References
* **Confidence Metric**: [[Local Distance Difference Test (LDDT & pLDDT)]].
* **Entities with Functional IDRs**: [[AKT1]], [[TP53]], [[STN7]].
* **Kinetic Transition Modeling**: [[Markov State Models]], [[Molecular Dynamics & Enhanced Sampling]].
