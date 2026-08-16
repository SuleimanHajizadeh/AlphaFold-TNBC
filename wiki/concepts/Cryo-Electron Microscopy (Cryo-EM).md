---
title: "Cryo-Electron Microscopy (Cryo-EM)"
type: concept
tags:
  - structural-biology/experimental-methods
  - single-particle-analysis
  - membrane-complexes
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Cryo-EM"
  - "Single-Particle Cryo-EM"
  - "Electron Microscopy"
---

# Cryo-Electron Microscopy (Cryo-EM)

**Cryogenic Electron Microscopy (Cryo-EM)** is an experimental structural biology technique that determines the 3D structures of macromolecular complexes, membrane proteins, and multi-subunit assemblies flash-frozen in vitreous (non-crystalline) ice at liquid nitrogen temperatures (~$77\text{ K}$).

---

## 🔬 Single-Particle Analysis (SPA) Workflow

```
Purified Macromolecules ──► Vitrification (Liquid Ethane) ──► Low-Dose TEM Imaging (Direct Electron Detector)
                                                                           │
                                                                           ▼
3D Density Map (Fourier Shell Correlation) ◄── 3D Classification ◄── 2D Class Averaging & Particle Picking
                    │
                    ▼
          Atomic Model Refinement
```

1. **Vitrification**: Rapid freezing prevents water crystallization, preserving proteins in near-native hydrated conformations.
2. **Direct Electron Detectors (DED)**: Detects electrons directly with high detective quantum efficiency (DQE) and movie-mode frame alignment to correct for beam-induced specimen movement.
3. **2D Class Averaging & Contrast Transfer Function (CTF)**: Thousands of individual projection images are aligned and averaged to boost signal-to-noise ratio (SNR).
4. **Fourier Shell Correlation (FSC)**: Calculates resolution (typically gold-standard FSC = 0.143 cutoff), frequently achieving near-atomic resolutions ($1.5\text{–}3.5\text{ Å}$).

---

## ⚡ Role in Membrane Proteins & Photosystems

Unlike X-ray crystallography which requires highly ordered crystal packing, Cryo-EM excels at solving:
* Large, flexible multi-subunit assemblies (e.g. [[Cytochrome b6f]], [[LHCII]]-PSII megacomplexes).
* Membrane-bound kinase complexes and receptors in nanodiscs or detergent micelles.
* Multiple distinct conformational states from a single biochemical sample (conformational heterogeneity).

---

## 🔗 Key Cross-References
* **Experimental Counterpart**: [[X-Ray Crystallography & Electron Density]].
* **Computational Prediction**: [[AlphaFold2]], [[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]].
* **Membrane Complexes Solved**: [[Cytochrome b6f]], [[LHCII]], [[PSII Reaction Center D1]].
* **Syntheses**: [[Structural Biology Methods - Experimental Cryo-EM and X-Ray vs Deep Learning AF2]].
