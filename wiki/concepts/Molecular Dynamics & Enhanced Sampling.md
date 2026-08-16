---
title: "Molecular Dynamics & Enhanced Sampling"
type: concept
tags:
  - molecular-dynamics/simulation
  - enhanced-sampling/metadynamics
  - free-energy-landscapes
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Molecular Dynamics"
  - "MD Simulations"
  - "Enhanced Sampling"
  - "Metadynamics"
  - "REMD"
---

# Molecular Dynamics & Enhanced Sampling

**Molecular Dynamics (MD)** is a computational method that simulates the time-dependent physical movements of atoms and molecules by numerically solving Newton’s classical equations of motion:

$$\mathbf{F}_i = -\nabla_i V(\mathbf{r}_1, \dots, \mathbf{r}_N) = m_i \frac{d^2 \mathbf{r}_i}{dt^2}$$

where $V(\mathbf{r})$ is the empirical potential energy function (force field).

---

## ⚡ Force Field Potential Energy Components

$$V_{\text{total}} = \sum_{\text{bonds}} k_b (r - r_0)^2 + \sum_{\text{angles}} k_\theta (\theta - \theta_0)^2 + \sum_{\text{dihedrals}} \frac{V_n}{2} [1 + \cos(n\phi - \gamma)] + \sum_{i < j} \left( 4\epsilon_{ij} \left[ \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6 \right] + \frac{q_i q_j}{4\pi\epsilon_0 \epsilon_r r_{ij}} \right)$$

1. **Bonded Interactions**: Bond stretching (harmonic), angle bending (harmonic), and dihedral angle torsions (Fourier series).
2. **Non-bonded Interactions**: Lennard-Jones 12-6 potential (van der Waals forces) and Coulomb potential (electrostatics).

---

## 🚀 Enhanced Sampling Techniques

Standard unbiased MD simulations often become trapped in local free energy minima due to high kinetic barriers ($\gg k_B T$). Enhanced sampling methods accelerate rare conformational transitions:

```
Unbiased MD (Trapped in local minimum)      Metadynamics (Fills well with repulsive Gaussian hills)
      /\                                                /\
  ___/  \___      ◄─── Kinetic Barrier ───►         ___/  \___ + ∑ Gaussians (CV)
 /          \                                      /          \
/  Min A     \  Min B                             /  Min A     \  Min B (Accessible!)
```

1. **Metadynamics**: Adds a history-dependent repulsive Gaussian bias potential along selected Collective Variables (CVs) to reconstruct the Free Energy Surface (FES):
   $$V_{\text{meta}}(\mathbf{s}, t) = \sum_{t' = \tau, 2\tau, \dots}^{t' < t} W \exp\left( -\sum_{i=1}^d \frac{(s_i - s_i(t'))^2}{2\sigma_i^2} \right)$$
2. **Replica Exchange MD (REMD)**: Runs multiple parallel non-interacting replicas at different temperatures, periodically attempting Monte Carlo exchanges to escape energetic traps.
3. **Coarse-Grained MD (Martini)**: Maps groups of ~4 heavy atoms to single interaction beads, extending simulation timescales from nanoseconds to microseconds and milliseconds.

---

## 🔗 Key Cross-References
* **Kinetic Modeling**: [[Markov State Models]].
* **Application to Kinase Plasticity**: [[Kinase Activation Loop & Allostery]], [[AKT1]], [[STN7]].
* **Course Curriculum Source**: [[Source - Structural Bioinformatics Course Curriculum]].
