---
title: "Molecular Docking Principles"
type: concept
tags:
  - drug-discovery/docking
  - protein-ligand/interactions
  - scoring-functions
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Molecular Docking"
  - "Protein-Ligand Docking"
  - "Docking Scoring Functions"
---

# Molecular Docking Principles

**Molecular Docking** is a computational method that predicts the preferred orientation (binding pose) and binding affinity of a small molecule ligand or protein partner within the active or allosteric site of a target macromolecule.

---

## ⚙️ Core Components

```
Ligand + Receptor ──► [ Search Algorithm ] ──► [ Scoring Function ] ──► Ranked Binding Poses (ΔG)
                        (Genetic, Monte Carlo)   (Physics / Empirical)
```

1. **Conformational Search Algorithms**:
   * **Genetic Algorithms (Lamarckian GA)**: Evolves a population of ligand translations, rotations, and torsional angles.
   * **Monte Carlo / Simulated Annealing**: Randomly perturbs ligand coordinates with Metropolis acceptance criterion: $P = \min(1, \exp(-\Delta E / k_B T))$.
2. **Scoring Functions**:
   * **Physics-based**: Computes molecular mechanics terms (van der Waals, Coulomb electrostatics, desolvation).
   * **Empirical (e.g. AutoDock Vina)**:
     $$\Delta G_{\text{binding}} = w_{\text{vdw}} \sum \text{Gauss} + w_{\text{rep}} \sum \text{Repulsion} + w_{\text{hbond}} \sum \text{HBond} + w_{\text{hydro}} \sum \text{Hydrophobic} + w_{\text{rot}} N_{\text{rot}}$$
   * **Machine Learning / Knowledge-based**: Uses graph neural networks or statistical contact potentials derived from PDB crystallographic co-complexes.

---

## 🎯 Rigid vs Flexible Receptor Docking

* **Rigid Receptor**: Receptor atoms remain fixed while the ligand degrees of freedom are sampled. High throughput, but ignores induced-fit adaptations.
* **Flexible Residue Docking**: Explicitly samples side-chain rotamers (e.g., active site catalytic triad or gatekeeper residue).
* **Ensemble Docking**: Docks ligands against diverse receptor conformations extracted from [[Molecular Dynamics & Enhanced Sampling|MD trajectories]] or [[AlphaFold2]] conformational ensembles.

---

## 🔗 Key Cross-References
* **Kinase Targeting**: [[Kinase Activation Loop & Allostery]], [[AKT1]].
* **Simulation Coupling**: [[Molecular Dynamics & Enhanced Sampling]], [[Markov State Models]].
* **Course Curriculum Source**: [[Source - Structural Bioinformatics Course Curriculum]].
