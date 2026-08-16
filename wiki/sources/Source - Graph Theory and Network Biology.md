---
title: "Source - Graph Theory and Network Biology"
type: source
tags:
  - source/network-biology
  - graph-theory/algorithms
  - systems-biology/wgcna
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Source - Graph Theory"
  - "Network Biology Source"
---

# Source - Graph Theory and Network Biology

* **Location**: `Structural Biology Tutorials/graph_theory.pdf`
* **Focus**: Graph-theoretical algorithms, network topology metrics, and their direct application to biomolecular systems (transcriptomic co-expression networks and residue interaction networks).

---

## 📌 Core Mathematical Formulations

A biological network is represented as an undirected or directed graph $G = (V, E)$ with adjacency matrix $\mathbf{A}$:

1. **Node Degree ($k_i$) & Hubs**:
   $$k_i = \sum_{j \in V} A_{ij}$$
   Biological networks exhibit scale-free power-law degree distributions $P(k) \sim k^{-\gamma}$, where a small number of hub nodes (e.g. [[AKT1]] in cancer networks) maintain global system integrity.
2. **Betweenness Centrality ($C_B(v)$)**:
   $$C_B(v) = \sum_{s \neq v \neq t} \frac{\sigma_{st}(v)}{\sigma_{st}}$$
   Measures the fraction of all shortest paths passing through node $v$. In **Residue Interaction Networks (RINs)**, residues with high betweenness centrality act as critical bottlenecks for allosteric signal transmission.
3. **Clustering Coefficient ($C_i$) & Modularity**:
   $$C_i = \frac{2 e_i}{k_i (k_i - 1)}$$
   Measures local interconnectivity. WGCNA (Weighted Gene Co-expression Network Analysis) leverages soft-thresholded adjacency to cluster genes into functional modules.

---

## 🔗 Key Cross-References
* **Concepts**: [[Residue Interaction Networks & Graph Centrality]], [[Contact Map Analysis]].
* **Identified Hub Entities**: [[AKT1]].
* **Syntheses**: [[Network to Structure - From WGCNA Hubs to Allosteric Kinase Inhibitors]].
