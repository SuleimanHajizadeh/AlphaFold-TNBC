#!/usr/bin/env python3
"""
protein_computational_analysis.py — Computational Biophysics Analysis
======================================================================
Implements:
  1. PDB coordinate parsing (Biopython)
  2. Kyte-Doolittle hydrophobicity scoring + hydrophobic core detection
  3. Miyazawa-Jernigan (MJ) statistical contact potential
  4. B-factor vs. Cα packing density Pearson correlation
  5. Three-panel biophysical profile figure (matplotlib)

Input : sample_data/1coh.pdb  (Collagen COL1A1 chain fragment)
Output: sample_data/biophysical_profile.png

Usage:
    python protein_computational_analysis.py

Dependencies: biopython, numpy, scipy, matplotlib

Author : Suleyman Hajizadeh
Purpose: Cambridge MPhil portfolio — demonstrating mathematical implementation
         of protein biophysics from first principles
"""

import os
import sys
import warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

# âââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
# PATHS
# âââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
SAMPLE_DIR  = os.path.join(SCRIPT_DIR, "sample_data")
PDB_FILE    = os.path.join(SAMPLE_DIR, "1coh.pdb")
OUTPUT_PNG  = os.path.join(SAMPLE_DIR, "biophysical_profile.png")
os.makedirs(SAMPLE_DIR, exist_ok=True)

# âââââââââââââââââââââââââââââââââââââââââ
<truncated 20987 bytes>