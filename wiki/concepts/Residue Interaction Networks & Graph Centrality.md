---
title: "Residue Interaction Networks & Graph Centrality"
type: concept
tags:
  - network-biology/structure
  - graph-centrality/allostery
  - residue-interaction-networks
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Residue Interaction Networks"
  - "RIN"
  - "Protein Structure Networks"
  - "Graph Centrality in Proteins"
---

# Residue Interaction Networks & Graph Centrality

A **Residue Interaction Network (RIN)** is a graph representation of a protein structure $G = (V, E)$, where each amino acid residue is a node $v \in V$, and edges $e \in E$ represent non-covalent physical interactions (hydrogen bonds, salt bridges, hydrophobic contacts, $\pi$-stacking).

RINs apply graph theory to analyze protein stability, folding cores, and long-range allosteric communication pathways.

---

## 📐 Graph Centrality Metrics in Structural Biology

```
       (Residue A) ─── (Residue B) ─── [ Central Bottleneck Node ] ─── (Residue D) ─── (Residue E)
                                                  │
                                                  ▼
                                       High Betweenness Centrality
                                    (Allosteric Transmission Conduit)
```

1. **Degree Centrality ($k_i$)**:
   * Number of direct physical contacts made by residue $i$. Residues with very high degree form the **rigid hydrophobic core** or catalytic clusters.
2. **Betweenness Centrality ($C_B(i)$)**:
   $$C_B(i) = \sum_{s \neq i \neq t} \frac{\sigma_{st}(i)}{\sigma_{st}}$$
   * Quantifies how frequently residue $i$ lies on the shortest structural communication paths between all other residue pairs $(s, t)$.
   * High betweenness nodes with moderate degree are premier candidates for **allosteric communication switches** (e.g. DFG motif residues in [[AKT1]]).
3. **Closeness Centrality ($C_C(i)$)**:
   $$C_C(i) = \frac{1}{\sum_{j \neq i} d(i, j)}$$
   * Measures structural accessibility and compact integration within the global fold.

---

## 🔗 Key Cross-References
* **Graph Foundations**: [[Source - Graph Theory and Network Biology]].
* **Distance Geometry Matrix**: [[Contact Map Analysis]].
* **Allosteric Signaling**: [[Kinase Activation Loop & Allostery]], [[AKT1]].
* **Integrative Synthesis**: [[Network to Structure - From WGCNA Hubs to Allosteric Kinase Inhibitors]].
