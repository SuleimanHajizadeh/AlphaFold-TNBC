# 🧬 Structural Analysis of TNBC Hub Kinase AKT1

[![PDB](https://img.shields.io/badge/PDB-4EJN-blue?style=flat-square)](https://www.rcsb.org/structure/4EJN)
[![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=flat-square&logo=python&logoColor=white)]()
[![Biopython](https://img.shields.io/badge/Biopython-1.79+-green?style=flat-square)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](https://opensource.org/licenses/MIT)

## 📌 Overview

This repository extends transcriptomic TNBC analysis into the **structural biology domain**. The 3D crystal structure of **AKT1** — a central hub kinase identified through our co-expression network and PI3K-Akt signalling hyperactivation analysis — is characterized at per-residue resolution using three complementary structural analyses.

> **AKT1** (*RAC-alpha serine/threonine-protein kinase*, UniProt: [P31749](https://www.uniprot.org/uniprot/P31749)) is a master regulator of cell survival and proliferation. Its hyperactivation is a hallmark of TNBC, driving chemotherapy resistance and tumour progression.

**Structure:** [PDB 4EJN](https://www.rcsb.org/structure/4EJN) — Crystal structure of autoinhibited AKT1 in complex with a kinase inhibitor. Resolution: **2.20 Å**.

---

## 🔬 Analyses Performed

| Script | Analysis | Output |
|--------|----------|--------|
| `analysis.py` | Per-residue B-factor profile with domain annotation | `figures/AKT1_Bfactor.png` |
| `ramachandran.py` | φ/ψ backbone torsion angles → secondary structure | `figures/AKT1_Ramachandran.png` |
| `contact_map.py` | Cα pairwise distance matrix → contact topology | `figures/AKT1_ContactMap.png` |

---

## 📊 Results

### 1. B-factor Flexibility Profile

Crystallographic B-factors reveal per-residue thermal motion and structural disorder across the four functional domains of AKT1.

| Region | B-factor | Interpretation |
|--------|----------|----------------|
| Rigid (< 20 Å²) | 0.0% | — |
| Moderate (20–40 Å²) | 15.7% | Well-ordered secondary structure |
| Flexible (40–60 Å²) | 65.9% | Loop regions, interdomain linkers |
| Disordered (> 60 Å²) | 18.4% | Regulatory domain, flexible termini |

![B-factor Plot](figures/AKT1_Bfactor.png)

> **Key finding:** 84.3% of AKT1 residues exhibit high flexibility (B > 40 Å²), consistent with the known conformational dynamics of the regulatory domain that underlie allosteric activation in cancer cells.

---

### 2. Ramachandran Backbone Analysis

φ/ψ torsion angles for all 366 non-terminal residues reveal the secondary structure composition of AKT1.

| Region | Count | Percentage |
|--------|-------|------------|
| α-helix favoured | 135 | **36.9%** |
| β-strand favoured | 121 | **33.1%** |
| Left-handed (Gly) | 8 | 2.2% |
| Other allowed | 102 | 27.9% |

![Ramachandran Plot](figures/AKT1_Ramachandran.png)

> **Key finding:** AKT1 has a mixed α/β fold (36.9% helix, 33.1% strand), characteristic of bilobal kinase architecture. The high proportion of "other allowed" residues (27.9%) reflects the extensive loop regions mediating substrate binding and allosteric communication.

---

### 3. Cα Contact Map

Pairwise residue contacts (threshold: 8 Å) reveal the topological organisation of AKT1 domains and long-range tertiary interactions.

| Metric | Value |
|--------|-------|
| Residues analysed | 376 |
| Contact pairs detected | 1,017 (1.4% of all pairs) |
| Min Cα–Cα distance | 2.94 Å |
| Max Cα–Cα distance | 63.88 Å |

![Contact Map](figures/AKT1_ContactMap.png)

> **Key finding:** The block-diagonal pattern confirms three distinct structural domains (PH, Kinase, HM). Off-diagonal contacts between the PH domain (residues 1–107) and the kinase domain (153–408) are sparse, consistent with AKT1's autoinhibited conformation captured in this crystal structure.

---

## 🗂️ Repository Structure

```
AlphaFold-TNBC/
├── data/
│   ├── AKT1_P31749.fasta      ← Canonical AKT1 sequence (UniProt P31749)
│   └── AKT1_ranked_0.pdb      ← Crystal structure (RCSB PDB: 4EJN, 2.20 Å)
├── figures/
│   ├── AKT1_Bfactor.png       ← B-factor domain profile
│   ├── AKT1_Ramachandran.png  ← φ/ψ torsion angle diagram
│   └── AKT1_ContactMap.png    ← Cα pairwise contact map
├── analysis.py                ← B-factor extraction and visualization
├── ramachandran.py            ← Backbone torsion angle analysis
├── contact_map.py             ← Cα distance matrix and contact map
└── README.md
```

---

## ⚙️ Reproduce

```bash
# Install dependencies
pip install biopython matplotlib numpy

# Run all analyses
python analysis.py        # B-factor profile
python ramachandran.py    # Ramachandran plot
python contact_map.py     # Contact map
```

---

## 🔗 Connection to TNBC Transcriptomics

This structural analysis directly complements the RNA-seq findings in the [Bioinformatics-analysis](https://github.com/SuleimanHajizadeh/Bioinformatics-analysis) repository:

- **Transcriptomic level:** AKT1 identified as a hub gene upregulated in TNBC via PI3K-Akt network hyperactivation
- **Structural level:** Crystal structure reveals the autoinhibited conformation and the high conformational flexibility of the regulatory HM domain — the mechanistic basis for its oncogenic activation

---

## 🎓 Academic Context

This project demonstrates a complete **multi-scale biological analysis pipeline**:

```
RNA-seq DESeq2  →  Co-expression Network  →  3D Crystal Structure
(gene level)        (systems level)           (atomic level)
```

Covering the full chain from gene expression to protein structural dynamics — a level of analytical depth characteristic of graduate-level computational biology research.

---

**Author:** Suleiman Hajizadeh | Bioinformatician @ IMBB, Azerbaijan  
📧 suleyman.hacizade1@gmail.com
