"""
AKT1 AlphaFold2 — PyMOL Visualization Script
Author : Suleiman Hajizadeh | IMBB, Azerbaijan

Usage (Windows CMD):
    pymol -c visualize_AKT1.py

Usage (inside PyMOL GUI):
    File -> Run Script -> select this file

Generates 3 publication-ready PNG renders in the figures/ folder:
  1. AKT1_pLDDT_3D.png   — B-factor/pLDDT spectrum (blue=high, red=low)
  2. AKT1_cartoon_3D.png — Standard grey cartoon view
  3. AKT1_surface_3D.png — Molecular surface with pLDDT colouring
"""

import os
import glob

# ── Find AlphaFold rank_001 PDB ───────────────────────────────────────────
pdb_hits = glob.glob("AKT1_TNBC_42642_0/*rank_001*.pdb")
if not pdb_hits:
    # fallback to crystal structure
    PDB_FILE = "data/AKT1_ranked_0.pdb"
else:
    PDB_FILE = pdb_hits[0]

OUTPUT_DIR = "figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"[i] Using PDB: {PDB_FILE}")

# ── PyMOL API ─────────────────────────────────────────────────────────────
try:
    from pymol import cmd
except ImportError:
    print("[!] PyMOL not found. Install with: conda install -c conda-forge pymol-open-source")
    raise

# ── Setup ─────────────────────────────────────────────────────────────────
cmd.load(PDB_FILE, "AKT1")
cmd.remove("solvent")
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)

# ── RENDER 1: pLDDT spectrum (B-factor column = pLDDT in AlphaFold) ───────
cmd.show_as("cartoon", "AKT1")
cmd.spectrum("b", "red_yellow_green_cyan_blue", "AKT1", minimum=0, maximum=100)
cmd.orient("AKT1")
cmd.zoom("AKT1", 5)
cmd.ray(1200, 900)
cmd.png(f"{OUTPUT_DIR}/AKT1_pLDDT_3D.png", dpi=300)
print(f"[✓] Render 1 → {OUTPUT_DIR}/AKT1_pLDDT_3D.png")

# ── RENDER 2: Clean cartoon — grey ────────────────────────────────────────
cmd.color("grey80", "AKT1")
cmd.orient("AKT1")
cmd.zoom("AKT1", 5)
cmd.ray(1200, 900)
cmd.png(f"{OUTPUT_DIR}/AKT1_cartoon_3D.png", dpi=300)
print(f"[✓] Render 2 → {OUTPUT_DIR}/AKT1_cartoon_3D.png")

# ── RENDER 3: Surface view with pLDDT colours ─────────────────────────────
cmd.show("surface", "AKT1")
cmd.set("surface_color", "grey70")
cmd.spectrum("b", "red_yellow_green_cyan_blue", "AKT1", minimum=0, maximum=100)
cmd.set("transparency", 0.3)
cmd.orient("AKT1")
cmd.zoom("AKT1", 3)
cmd.ray(1200, 900)
cmd.png(f"{OUTPUT_DIR}/AKT1_surface_3D.png", dpi=300)
print(f"[✓] Render 3 → {OUTPUT_DIR}/AKT1_surface_3D.png")

print("\n[✓] All renders complete! Check the figures/ folder.")
cmd.quit()
