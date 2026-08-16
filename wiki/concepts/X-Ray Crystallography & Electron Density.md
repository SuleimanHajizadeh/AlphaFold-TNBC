---
title: "X-Ray Crystallography & Electron Density"
type: concept
tags:
  - structural-biology/experimental-methods
  - crystallography/diffraction
  - electron-density
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "X-Ray Crystallography"
  - "Electron Density Maps"
  - "Bragg's Law"
---

# X-Ray Crystallography & Electron Density

**X-Ray Crystallography** is the foundational experimental method that determines the atomic coordinates of macromolecules by scattering X-ray beams through a periodic 3D protein crystal.

---

## 📐 Physical Principles & Fourier Transform

```
Protein Crystal ──► X-Ray Beam (λ ≈ 1 Å) ──► Diffraction Pattern (Spots / Intensities |F(hkl)|²)
                                                      │
                                                      ▼ + Phase Information (MR / SAD / MAD)
Atomic Coordinate Model ◄── Electron Density Map ρ(x,y,z) = 1/V ∑ |F| exp[iα] exp[-2πi(hx+ky+lz)]
```

1. **Bragg's Law**: Constructive interference occurs when:
   $$n\lambda = 2d \sin\theta$$
2. **The Phase Problem**: Detectors record diffraction intensities $I(hkl) \propto |F(hkl)|^2$, losing phase angles $\alpha(hkl)$. Phases are retrieved via:
   * **Molecular Replacement (MR)**: Using homologous search models (now dramatically accelerated by [[AlphaFold2]] models).
   * **Experimental Phasing**: SAD/MAD using anomalous dispersion from heavy atoms or selenium.
3. **B-Factors (Temperature Factors)**:
   * Reflect thermal isotropic vibration and static crystal disorder:
     $$B = 8\pi^2 \langle u^2 \rangle$$
   * In [[AlphaFold2]] PDB outputs, the B-factor column is repurposed to store per-residue [[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]] scores ($0\text{–}100$).

---

## 🔗 Key Cross-References
* **Alternative Experimental Method**: [[Cryo-Electron Microscopy (Cryo-EM)]].
* **Computational Prediction**: [[AlphaFold2]], [[Local Distance Difference Test (LDDT & pLDDT)]].
* **Geometric Validation**: [[Ramachandran Dihedral Angles]].
* **Syntheses**: [[Structural Biology Methods - Experimental Cryo-EM and X-Ray vs Deep Learning AF2]].
