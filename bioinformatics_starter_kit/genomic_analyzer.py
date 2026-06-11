"#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ð§¬ GENOMIC DATA ANALYZER & MUTATION DETECTOR (Hamming Distance & Translation)
==============================================================================
Bu alÉt bioloji FASTA verilÉnlÉrini analiz edÉn, transkripsiya/translasiya aparan
vÉ DNT zÉncirlÉri arasÄ±ndakÄ± mutasiyalarÄ± (SNPs) Hamming mÉsafÉsi alqoritmi ilÉ
hesablayan vÉ terminalda rÉngli ÅÉkildÉ vizuallaÅdÄ±ran professional bioinformatika skriptidir.

MÃ¼Éllif: SÃ¼leyman HacÄ±zadÉ (Hybrid Portfolio)
Tarix: 2026-05-31
==============================================================================
"""

import os
import sys
import argparse

# Terminalda rÉngli Ã§Ä±xÄ±Ålar Ã¼Ã§Ã¼n ANSI kodlarÄ±
COLOR_GREEN = "\033[92m"
COLOR_RED = "\033[91m"
COLOR_YELLOW = "\033[93m"
COLOR_BLUE = "\033[94m"
COLOR_BOLD = "\033[1m"
COLOR_RESET = "\033[0m"

# Standart Genetik Kodon CÉdvÉli (RNA -> Amino Acid)
CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',  # * = STOP kodonlar
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',  # M = START/Metionin
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def parse_fasta(file_path):
    """
    FASTA formatlÄ± faylÄ± oxuyur vÉ baÅlÄ±qlar ilÉ ardÄ±cÄ±llÄ±qlarÄ± lÃ¼ÄÉt (dict) olaraq qaytarÄ±r.
<truncated 8048 bytes>