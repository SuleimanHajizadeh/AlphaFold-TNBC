---
title: "Markov State Models"
type: concept
tags:
  - molecular-dynamics/kinetics
  - markov-chains
  - conformational-transitions
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Markov State Models"
  - "MSM"
  - "Kinetic Transition Networks"
---

# Markov State Models (MSMs)

**Markov State Models (MSMs)** are mathematical frameworks used to model the long-timescale conformational kinetics and thermodynamics of biomolecules from ensembles of short, distributed [[Molecular Dynamics & Enhanced Sampling|Molecular Dynamics simulations]].

---

## 📐 Mathematical Formulation

Conformational space is discretized into $k$ discrete microstates $\{S_1, S_2, \dots, S_k\}$. The evolution of the system is modeled as a memoryless (Markovian) jump process governed by the **Transition Probability Matrix** $\mathbf{T}(\tau)$:

$$T_{ij}(\tau) = P\left( x(t+\tau) \in S_j \mid x(t) \in S_i \right)$$

where $\tau$ is the lag time.

### Master Equation & Stationary Distribution
The stationary distribution $\boldsymbol{\pi}$ satisfies detailed balance (reversibility):

$$\pi_i T_{ij}(\tau) = \pi_j T_{ji}(\tau) \quad \text{and} \quad \boldsymbol{\pi} \mathbf{T}(\tau) = \boldsymbol{\pi}$$

### Implied Timescales:
The relaxation timescales $t_i$ associated with the slowest conformational transitions (eigenvalues $\lambda_i$) are given by:

$$t_i = -\frac{\tau}{\ln |\lambda_i(\tau)|}$$

---

## 🔬 MSM Construction Pipeline

1. **Dimensionality Reduction**: Time-lagged Independent Component Analysis (tICA) projects high-dimensional coordinate trajectories onto coordinates with maximum autocorrelation.
2. **Spatial Clustering**: $k$-means clustering partitions tICA space into hundreds of microstates.
3. **Transition Counting & Estimation**: Maximum Likelihood Estimation (MLE) computes the transition matrix $\mathbf{T}(\tau)$.
4. **Coarse-Graining (PCCA+)**: Perron-Cluster Cluster Analysis groups microstates into kinetically distinct macrostates (e.g. Inactive, Intermediate, Active).

---

## 🔗 Key Cross-References
* **Underlying Trajectories**: [[Molecular Dynamics & Enhanced Sampling]].
* **Biological Switches**: [[Kinase Activation Loop & Allostery]], [[Thylakoid State Transitions]].
* **Course Curriculum Source**: [[Source - Structural Bioinformatics Course Curriculum]].
