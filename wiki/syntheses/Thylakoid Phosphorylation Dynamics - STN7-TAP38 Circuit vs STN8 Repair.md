---
title: "Thylakoid Phosphorylation Dynamics - STN7-TAP38 Circuit vs STN8 Repair"
type: synthesis
tags:
  - synthesis/plant-biophysics
  - photosynthesis/state-transitions
  - thylakoid/phosphorylation-circuits
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Thylakoid Kinase Phosphatase Circuits"
  - "STN7 vs STN8 Functional Specialization"
---

# Thylakoid Phosphorylation Dynamics: STN7-TAP38 Circuit vs STN8 PSII Repair

This synthesis integrates the structural, evolutionary, and biophysical evidence characterizing the two distinct protein kinase/phosphatase regulatory networks of the chloroplast thylakoid membrane.

---

## ⚡ Comparative Network Architecture

```
          REDOX-DRIVEN STATE TRANSITION CIRCUIT                      HIGH-LIGHT PHOTOPROTECTION CIRCUIT
          ─────────────────────────────────────                      ──────────────────────────────────

       Excess PSII Light ──► PQ Pool Reduced (PQH2)               Excess Light Stress ──► Direct D1 Photo-oxidation
                                   │                                                            │
                                   ▼                                                            ▼
            [[Cytochrome b6f]] (Qo site activation)                                      [[STN8]] Kinase
                                   │                                                            │
                                   ▼                                                            ▼
                            [[STN7]] Kinase                                      Phosphorylates [[PSII Reaction Center D1]]
                                   │                                                            │
                                   ▼                                                            ▼
                    Phosphorylates [[LHCII]] (Thr-3)                              PSII Supercomplex Disassembly
                                   │                                                            │
                                   ▼                                                            ▼
           [[Thylakoid State Transitions|State 2 Transition (LHCII -> PSI)]]               Grana Unstacking & FtsH Cleavage
                                   │                                                            │
                                   ▼                                                            ▼
                    [[TAP38-PPH1]] Phosphatase                                   PBCP Phosphatase / D1 Re-synthesis
                                   │                                                            │
                                   ▼                                                            ▼
                          State 1 Restoration                                         Restored Active PSII
```

---

## 🔬 Evolutionary & Structural Specialization

1. **Ancient Paralog Duplication**:
   * [[STN7]] and [[STN8]] share **36.12%** sequence identity in *Arabidopsis thaliana* with **57.02%** perfect residue conservation across green algae, monocots, and dicots.
   * Superposition of their stromal catalytic cores via SVD yields an **RMSD of $4.2507\text{ Å}$** across 173 matched $C_\alpha$ positions ([[Structural Superposition & RMSD]]).
2. **Substrate Selectivity**:
   * **STN7** is evolutionary specialized for antenna mobility (phosphorylating Lhcb1/Lhcb2 in **[[LHCII]]**). Its activity is tightly gated by plastoquinol binding to the $Q_o$ site of **[[Cytochrome b6f]]**, and counteracted by the phosphatase **[[TAP38-PPH1]]**.
   * **STN8** is specialized for core protein turnover (phosphorylating Thr-2 on **[[PSII Reaction Center D1]]**, D2, and CP43). It enables unstacking of dense grana membranes to allow bulky stromal **FtsH proteases** access to damaged D1 subunits.
3. **Membrane Association Mode**:
   * Neither kinase contains a canonical hydrophobic transmembrane span (maximum Kyte-Doolittle scores are 1.49 for STN7 and 1.59 for STN8), suggesting a peripheral association mediated by amphipathic surface helices or docking to adjacent membrane complexes.

---

## 🔗 Key Cross-References
* **Kinase Entities**: [[STN7]], [[STN8]].
* **Phosphatases**: [[TAP38-PPH1]].
* **Substrates & Sensors**: [[LHCII]], [[Cytochrome b6f]], [[PSII Reaction Center D1]].
* **Core Biological Mechanisms**: [[Thylakoid State Transitions]], [[Kinase Activation Loop & Allostery]].
* **Primary Source**: [[Source - STN7-STN8 Docking Project]].
