# 🧬 Structural Characterization of TNBC Biomarkers via AlphaFold2

[![AlphaFold](https://img.shields.io/badge/Tool-AlphaFold2_|_ColabFold-darkblue?style=flat-square)](https://github.com/sokrypton/ColabFold)
[![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=flat-square&logo=python&logoColor=white)]()
[![Biopython](https://img.shields.io/badge/Biopython-Structural_Analysis-blue?style=flat-square)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](https://opensource.org/licenses/MIT)

## 📌 Overview

This repository extends the transcriptomic analysis of Triple-Negative Breast Cancer (TNBC) into the **structural biology domain**. Using **AlphaFold2 (via ColabFold)**, the 3D structure of **AKT1** — a key hub kinase identified in our TNBC co-expression network and PI3K-Akt signalling hyperactivation analysis — is predicted and characterized.

> **AKT1** (*RAC-alpha serine/threonine-protein kinase*, UniProt: P31749) plays a central oncogenic role in TNBC by promoting cell survival, proliferation, and chemotherapy resistance. Structural insight into its kinase domain informs therapeutic target identification.

---

## 🔬 Methodology

| Step | Method | Tool |
|------|--------|------|
| 1. Sequence retrieval | Canonical human AKT1 FASTA | UniProt (P31749) |
| 2. 3D structure prediction | Deep learning MSA-based structure prediction | AlphaFold2 / ColabFold |
| 3. Quality assessment | Per-residue pLDDT confidence score extraction from B-factor column | Biopython + NumPy |
| 4. Visualization (structure) | Publication-ready 3D molecular rendering | PyMOL |
| 5. Visualization (confidence) | pLDDT colour-coded bar plot | Matplotlib |

---

## 📊 Key Results

### pLDDT Confidence Profile

The per-residue **pLDDT (Predicted Local Distance Difference Test)** score from AlphaFold2 quantifies the local structural confidence:
- **≥ 90** → Very high confidence (blue)
- **70 – 90** → Confident (light blue)
- **50 – 70** → Low confidence (yellow)
- **< 50** → Very low / intrinsically disordered region (orange)

![pLDDT Plot](figures/AKT1_pLDDT.png)

---

## 🗂️ Repository Structure

```
AlphaFold-TNBC/
├── data/
│   ├── AKT1_P31749.fasta         ← Input sequence (UniProt canonical)
│   └── AKT1_ranked_0.pdb         ← AlphaFold2 top-ranked predicted structure
├── figures/
│   └── AKT1_pLDDT.png            ← Per-residue confidence plot (auto-generated)
├── results/                       ← Summary statistics output
├── analysis.py                    ← Main pLDDT analysis and plotting script
└── README.md
```

---

## ⚙️ Reproduce This Analysis

### 1. Install dependencies
```bash
conda create -n alphafold-env python=3.10 biopython matplotlib numpy
conda activate alphafold-env
```

### 2. Obtain the structure
Run **ColabFold** (Google Colab, no GPU required locally):
→ [https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)

Paste `data/AKT1_P31749.fasta` → Run → Download `ranked_0.pdb` → place in `data/AKT1_ranked_0.pdb`.

### 3. Run analysis
```bash
python analysis.py
```

---

## 🎓 Academic Context
This work demonstrates proficiency in **structural bioinformatics**, extending transcriptomic findings into 3D molecular space. It complements the RNA-seq TNBC pipeline in the [Bioinformatics-analysis](https://github.com/SuleimanHajizadeh/Bioinformatics-analysis) repository, covering the full chain from gene expression → network analysis → protein structure prediction.

---

**Author:** Suleiman Hajizadeh | Bioinformatician @ IMBB, Azerbaijan
📧 suleyman.hacizade1@gmail.com
