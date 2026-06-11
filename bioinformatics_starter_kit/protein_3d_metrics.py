"#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ð¬ PROTEIN 3D STRUCTURE METRICS & BACKBONE VISUALIZER (PDB Parser & Euclidean Distance)
==============================================================================
Bu alÉt RCSB PDB bazasÄ±ndan protein strukturlarÄ±nÄ± avtomatik yÃ¼klÉyÉn, Bio.PDB
ilÉ atom koordinatlarÄ±nÄ± oxuyan, iki amin turÅusu arasÄ±ndakÄ± 3D Evklid mÉsafÉsini
hesablayan vÉ protein zÉncirini (backbone) 3D olaraq Matplotlib ilÉ vizuallaÅdÄ±ran
professional bioinformatika skriptidir.

MÃ¼Éllif: SÃ¼leyman HacÄ±zadÉ (Hybrid Portfolio)
Tarix: 2026-05-31
==============================================================================
"""

import os
import sys
import math
import argparse
import urllib.request

# Terminalda rÉngli Ã§Ä±xÄ±Ålar Ã¼Ã§Ã¼n ANSI kodlarÄ±
COLOR_GREEN = "\033[92m"
COLOR_RED = "\033[91m"
COLOR_YELLOW = "\033[93m"
COLOR_BLUE = "\033[94m"
COLOR_BOLD = "\033[1m"
COLOR_RESET = "\033[0m"

# KitabxanalarÄ±n yoxlanÄ±lmasÄ±
try:
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from Bio.PDB import PDBParser
except ImportError:
    print(f"{COLOR_RED}[XÆTA] LazÄ±m olan kitabxanalar tapÄ±lmadÄ±!{COLOR_RESET}")
    print(f"ZÉhmÉt olmasa asÄ±lÄ±lÄ±qlarÄ± quraÅdÄ±rÄ±n: {COLOR_BOLD}pip install -r requirements.txt{COLOR_RESET}")
    sys.exit(1)

def download_pdb(pdb_id, output_dir="sample_data"):
    """
    RCSB Protein Data Bank bazasÄ±ndan verilÉn PDB ID-yÉ uyÄun .pdb faylÄ±nÄ± yÃ¼klÉyir.
    """
    pdb_id = pdb_id.lower()
    os.makedirs(output_dir, exist_ok=True)
    file_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    
    if os.path.exists(file_path):
        print(f"   -> [OK] Protein faylÄ± lokalda mÃ¶vcuddur: {file_path}")
        return file_path
        
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"   -> {COLOR_YELLOW}[YÃKLÆNÄ°R]{COLOR_RESET} PDB faylÄ± yÃ¼klÉni
<truncated 9700 bytes>