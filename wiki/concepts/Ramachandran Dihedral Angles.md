---
title: "Ramachandran Dihedral Angles"
type: concept
tags:
  - structural-biology/geometry
  - protein-backbone
  - torsion-angles
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Ramachandran Plot"
  - "Dihedral Angles"
  - "Phi Psi Angles"
---

# Ramachandran Dihedral Angles

The **Ramachandran Plot** (Ramachandran, Ramakrishnan, & Sasisekharan, 1963) maps the distribution of the peptide backbone dihedral torsion angles **$\phi$ (phi)** and **$\psi$ (psi)**, describing the stereochemical conformations accessible to amino acid residues in polypeptide chains without steric clashes.

---

## 📐 Vector Algebra Derivation

For a peptide chain with sequential Cartesian coordinates $\mathbf{r}_1, \mathbf{r}_2, \mathbf{r}_3, \mathbf{r}_4$:

* **$\phi_i$**: Dihedral angle between planes $C_{i-1} - N_i - C_\alpha_i$ and $N_i - C_\alpha_i - C_i$.
* **$\psi_i$**: Dihedral angle between planes $N_i - C_\alpha_i - C_i$ and $C_\alpha_i - C_i - N_{i+1}$.

```
      H   H   O       H   H   O
      |   |   ||      |   |   ||
   - -N - Ca- C - - - N - Ca- C - -
        \   /           \   /
        phi_i           psi_i
```

### Computation via Cross-Products:
1. Form displacement vectors:
   $$\mathbf{v}_1 = \mathbf{r}_2 - \mathbf{r}_1, \quad \mathbf{v}_2 = \mathbf{r}_3 - \mathbf{r}_2, \quad \mathbf{v}_3 = \mathbf{r}_4 - \mathbf{r}_3$$
2. Determine plane normals:
   $$\mathbf{n}_1 = \mathbf{v}_1 \times \mathbf{v}_2, \quad \mathbf{n}_2 = \mathbf{v}_2 \times \mathbf{v}_3$$
3. Compute normalized axis vector:
   $$\mathbf{u} = \frac{\mathbf{v}_2}{|\mathbf{v}_2|}$$
4. Calculate angle $\theta \in [-\pi, \pi]$:
   $$\theta = \mathrm{atan2}\left( (\mathbf{n}_1 \times \mathbf{n}_2) \cdot \mathbf{u}, \ \mathbf{n}_1 \cdot \mathbf{n}_2 \right)$$

---

## 🗺️ Energetically Allowed Regions

* **$\alpha$-helical region**: $\phi \approx -60^\circ, \psi \approx -45^\circ$ (Right-handed $\alpha$-helices).
* **$\beta$-sheet region**: $\phi \approx -120^\circ\text{ to }-140^\circ, \psi \approx +130^\circ\text{ to }+150^\circ$ (Parallel and anti-parallel $\beta$-strands).
* **Left-handed $\alpha$-helix**: $\phi \approx +60^\circ, \psi \approx +45^\circ$ (Predominantly Glycine due to lack of a $C_\beta$ side chain).

---

## 🔗 Key Cross-References
* **Project Validation**: [[Source - AKT1 Kinase Modeling Project]] ($\alpha$-helix: 24.6%, $\beta$-sheet: 32.9%, loops: 39.4%).
* **Structure Evaluation**: [[Local Distance Difference Test (LDDT & pLDDT)|pLDDT]], [[Contact Map Analysis]].
