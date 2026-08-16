---
title: "Thylakoid State Transitions"
type: concept
tags:
  - plant-biology/photosynthesis
  - thylakoid/dynamics
  - light-harvesting
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "State Transitions"
  - "Photosynthetic Acclimation"
  - "LHCII Phosphorylation"
---

# Thylakoid State Transitions

**State Transitions** represent a rapid physiological mechanism in photosynthetic organisms (plants and algae) that dynamically balances light excitation energy between **Photosystem II (PSII)** and **Photosystem I (PSI)** in response to changing light spectrum and intensity.

---

## ⚡ Molecular Mechanism

```
State 1 (PSII Favored)                            State 2 (PSI Favored)
   [ PSII - LHCII ]                                 [ PSI - P-LHCII ]
          │                                                 ▲
          │ Plastoquinone (PQ) pool reduced                 │
          ▼                                                 │
   STN7 Kinase Activation  ────────────────────────►  LHCII Phosphorylated & Detached
```

1. **State 1 $\rightarrow$ State 2 (STN7 Kinase Dependent)**:
   * Preferential PSII illumination causes reduction of the plastoquinone (PQ) pool to plastoquinol ($\text{PQH}_2$).
   * Binding of $\text{PQH}_2$ to the cytochrome $b_6f$ complex activates the thylakoid kinase [[STN7]].
   * STN7 phosphorylates the N-terminal threonine residues of Light-Harvesting Complex II (LHCII) trimers (Lhcb1/Lhcb2).
   * Phosphorylation introduces negative charges, causing LHCII to dissociate from PSII in the appressed grana thylakoids and migrate to PSI in non-appressed stroma lamellae.
2. **State 2 $\rightarrow$ State 1 (TAP38 / PPH1 Phosphatase Dependent)**:
   * When PSI is preferentially excited or in the dark, the PQ pool becomes oxidized, inactivating STN7.
   * Constitutively active phosphatase **TAP38 / PPH1** dephosphorylates LHCII, causing it to return to PSII.
3. **Contrast with PSII Repair Kinase [[STN8]]**:
   * While STN7 phosphorylates antenna complexes (LHCII), its paralog [[STN8]] phosphorylates PSII core subunits (D1, D2, CP43) to govern photo-damage repair.

---

## 🔗 Key Cross-References
* **Kinase Entities**: [[STN7]], [[STN8]].
* **Structural Pipeline**: [[Source - STN7-STN8 Docking Project]].
* **Comparative Syntheses**: [[Comparative Kinase Mechanics - AKT1 vs STN7-STN8]].
