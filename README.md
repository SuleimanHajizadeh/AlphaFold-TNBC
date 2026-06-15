# 🧬 Computational Structural Biology Laboratory
## Advanced Pipelines for Deep Learning-Based Structure Prediction, Structural Analytics, and Automated Molecular Docking

[![AlphaFold2](https://img.shields.io/badge/Model-AlphaFold2_|_ColabFold-darkblue?style=flat-square)](https://github.com/sokrypton/ColabFold)
[![AutoDock Vina](https://img.shields.io/badge/Docking-AutoDock_Vina_v1.2.5-darkred?style=flat-square)](https://github.com/ccsb-scripps/AutoDock-Vina)
[![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=flat-square&logo=python&logoColor=white)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-green?style=flat-square)](https://opensource.org/licenses/MIT)

This repository serves as a centralized laboratory platform showcasing advanced computational structural biology and molecular modeling workflows. It integrates state-of-the-art deep learning methods for 3D macromolecular structure prediction, custom vector-algebraic structural analysis, and automated high-throughput virtual screening pipelines.

---

## 🗂️ Lab Portfolio & Project Registry

| Project Name | 📁 Directory Path | 🔬 Key Methods & Technologies | 🔗 Details |
|--------------|-------------------|-------------------------------|------------|
| **AKT1 Kinase Structural Modeling** | [`/projects/akt1-kinase-modeling`](./projects/akt1-kinase-modeling) | AlphaFold2, pLDDT & PAE matrices, Ramachandran dihedral torsion vectors, $C_\alpha$ contact maps | [README](./projects/akt1-kinase-modeling/README.md) |
| **STN7/STN8 Chloroplast Kinase Docking** | [`/projects/stn7-stn8-docking`](./projects/stn7-stn8-docking) | Receptor preparation, signal transit peptide truncation, AutoDock Vina, grid box auto-parameterization | [README](./projects/stn7-stn8-docking/README.md) |

---

## 🔬 Core Pipelines & Workflows

### 1. AKT1 Kinase Structural Modeling & Analytics
Predicts and validates the 3D structure of the Triple-Negative Breast Cancer (TNBC) oncogenic hub kinase **AKT1** using AlphaFold2. It utilizes custom Python routines to characterize prediction confidence, secondary structures, and topology.

```mermaid
graph TD
    A["FASTA Sequence (UniProt P31749)"] --> B["ColabFold v1.6.1 MMseqs2 MSA"]
    B --> C["AlphaFold2 Structure Module"]
    C --> D["Top Rank PDB Model"]
    D --> E["pLDDT Profiling (Local Distance Difference Test)"]
    D --> F["Dihedral Torsion Vectors (Ramachandran Analysis)"]
    D --> G["Cα Distance Contact Map (Euclidean Topology)"]
    D --> H["Predicted Aligned Error (PAE Modeling)"]
```

---

### 2. STN7/STN8 Chloroplast Kinase Docking Pipeline
Automates the structural filtration, receptor preparation, and high-throughput virtual screening of allosteric ligands targeting the Arabidopsis thylakoid membrane protein kinases **STN7** and **STN8**.

```mermaid
graph TD
    A["Raw Receptor Kinase PDB"] --> B["prepare_receptor.py: Water removal & heteroatom isolation"]
    B --> C["Transit peptide sequence truncation"]
    C --> D["AutoDock Vina grid box calculations"]
    D --> E["High-throughput docking (Vina 1.2.5)"]
    E --> F["Ensemble pose scoring & binding affinity mapping"]
```

---

## 🚀 Environment & Setup

Both pipelines run on a unified Python stack. Set up the Conda environment using:

```bash
conda create -n struct-bio python=3.10 -y
conda activate struct-bio
pip install numpy matplotlib biopython pandas pytest scipy
```

To run molecular docking, ensure **AutoDock Vina** and **MGLTools** (`prepare_receptor4.py`) are installed in your PATH.

---

## 🧪 Automated Testing

Unit tests for validation logic and receptor preparation are included. Run them using:

```bash
# Test AKT1 analysis logic
pytest projects/akt1-kinase-modeling/test_analysis.py

# Test STN7/STN8 receptor preparation
pytest projects/stn7-stn8-docking/scripts/docking/test_prepare_receptor.py
```

---

**Author:** Suleiman Hajizadeh  
📧 suleyman.hacizade1@gmail.com | 🔗 [GitHub Profile](https://github.com/SuleimanHajizadeh)
