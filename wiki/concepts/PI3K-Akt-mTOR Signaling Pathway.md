---
title: "PI3K-Akt-mTOR Signaling Pathway"
type: concept
tags:
  - oncology/signaling
  - kinase-cascade
  - breast-cancer/tnbc
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "PI3K-Akt Pathway"
  - "PI3K/Akt/mTOR"
  - "Akt Signaling Cascade"
---

# PI3K-Akt-mTOR Signaling Pathway

The **PI3K/Akt/mTOR signaling cascade** is one of the most critical intracellular signal transduction pathways in human physiology, orchestrating cell growth, survival, glucose metabolism, protein synthesis, motility, and angiogenesis.

Hyperactivation of this pathway occurs in over 70% of human cancers, and represents a premier molecular driver in **Triple-Negative Breast Cancer (TNBC)**.

---

## ⚡ Canonical Signaling Hierarchy

```
Receptor Tyrosine Kinase (RTK) / GPCR
                 │
                 ▼
     PI3K (Class IA: p110α/p85)
                 │  PIP2 ──► PIP3 (Lipid Second Messenger)
                 ▼         ▲
  [ PIP3 Membrane Docking ] ── PTEN (Tumor Suppressor Phosphatase)
        │             │
        ▼             ▼
      PDK1          AKT1 (via PH Domain)
        │             │
        └─────► Thr308 Phosphorylation (Activation Loop)
                      │
                      ├◄──── Ser473 Phosphorylation (mTORC2 / Hydrophobic Motif)
                      ▼
               Fully Active AKT1
                      │
      ┌───────────────┼───────────────┬───────────────┐
      ▼               ▼               ▼               ▼
GSK3β (Inhibited)  FOXO (Excluded)  TSC2 (Inhibited)  Bad (Inhibited)
(Glycogen/Survival)(Anti-apoptosis) (mTORC1 active)  (Apoptosis blocked)
```

---

## 🔬 Molecular Steps & Cross-Talk

1. **Lipid Kinase Activation**: Ligand-bound RTKs recruit Class IA Phosphoinositide 3-Kinases (PI3K), which phosphorylate $\text{PIP}_2$ into $\text{PIP}_3$.
2. **Co-localization & PH Domain Recruitment**: Both **PDK1** and **[[AKT1]]** possess Pleckstrin Homology (PH) domains that bind membrane-anchored $\text{PIP}_3$, bringing both kinases into direct spatial proximity.
3. **Dual Phosphorylation Trigger**:
   * **PDK1** phosphorylates Thr308 in the [[Kinase Activation Loop & Allostery|Activation Loop]] of AKT1.
   * **mTORC2** (mechanistic target of rapamycin complex 2) phosphorylates Ser473 in the C-terminal hydrophobic motif, conferring maximum catalytic output.
4. **PTEN Antagonism**: The tumor suppressor **PTEN** (Phosphatase and tensin homolog) acts as the primary negative regulator by dephosphorylating $\text{PIP}_3$ back to $\text{PIP}_2$. Loss of PTEN is a frequent cause of constitutive AKT1 hyperactivation in TNBC.

---

## 🔗 Key Cross-References
* **Central Kinase Hub**: [[AKT1]].
* **Structural Activation**: [[Kinase Activation Loop & Allostery]].
* **Oncogenic Mutational Context**: [[TP53]].
* **Drug Targeting & Modeling**: [[Allosteric Drug Discovery & Conformational Selection in Kinases]], [[Source - AKT1 Kinase Modeling Project]].
