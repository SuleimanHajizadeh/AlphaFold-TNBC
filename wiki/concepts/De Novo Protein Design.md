---
title: "De Novo Protein Design"
type: concept
tags:
  - protein-design/generative-ai
  - structural-biology/deep-learning
  - rfdiffusion
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Protein Design"
  - "De Novo Design"
  - "RFdiffusion"
  - "ProteinMPNN"
---

# De Novo Protein Design

**De Novo Protein Design** involves creating entirely artificial, functional proteins from scratch without relying on naturally occurring homologous templates. Modern deep learning has revolutionized the field through continuous diffusion generative models and neural inverse folding.

---

## 🏗️ Deep Learning Design Paradigm

```
Target Functional Motif / Binding Site
                 │
                 ▼
1. Backbone Generation (RFdiffusion / SE(3) Diffusion)
   Transforms random 3D Gaussian noise into compact, designable protein backbones.
                 │
                 ▼
2. Sequence Design & Inverse Folding (ProteinMPNN)
   Autoregressively predicts sequences that fold into the generated target backbone.
                 │
                 ▼
3. In Silico Filtering & Self-Consistency Validation (AlphaFold2 / ESMFold)
   Folds predicted sequences: requires RMSD < 2.0 Å, pLDDT > 80, and high contact satisfaction.
```

---

## ⚙️ Core Architectures

1. **RFdiffusion (Watson et al., 2023)**:
   * Fine-tunes RoseTTAFold structure modules into an $SE(3)$ equivariant denoising diffusion probabilistic model (DDPM).
   * Generates backbones conditioned on functional motifs, symmetric oligomers, and binder interfaces.
2. **ProteinMPNN (Dauparas et al., 2022)**:
   * Message-passing neural network trained on the Protein Data Bank (PDB) to solve the inverse folding problem: $P(\text{Sequence} \mid \text{Backbone})$.
3. **In Silico Success Metrics**:
   * **Self-Consistency RMSD (scRMSD)**: $\text{scRMSD} < 2.0\text{ Å}$ between designed backbone and [[AlphaFold2]] refolded model.
   * **sc-pLDDT**: $\geq 80\text{–}85$ across designed binding interfaces.

---

## 🔗 Key Cross-References
* **Validation & Refolding Engine**: [[AlphaFold2]], [[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]], [[Structural Superposition & RMSD]].
* **Conformational Testing**: [[Molecular Dynamics & Enhanced Sampling]], [[Molecular Docking Principles]].
* **Course Curriculum Source**: [[Source - Structural Bioinformatics Course Curriculum]].
