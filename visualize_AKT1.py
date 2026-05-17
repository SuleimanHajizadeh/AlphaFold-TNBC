"""
AKT1 (4EJN) — PyMOL Visualization Script
Author : Suleiman Hajizadeh | IMBB, Azerbaijan
Usage  : Run inside PyMOL: File → Run Script → select this file
         OR from terminal: pymol -c visualize_AKT1.py

Generates 3 publication-ready PNG renders:
  1. Full structure coloured by B-factor (flexibility)
  2. Domain-coloured cartoon (PH / Kinase / Regulatory)
  3. Active-site zoomed view with inhibitor
"""

import os

# ── PyMOL API ──────────────────────────────────────────────────────────────
try:
    from pymol import cmd
except ImportError:
    print("[!] This script must be run inside PyMOL.")
    print("    Launch with: pymol -c visualize_AKT1.py")
    raise

OUTPUT_DIR = "figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

PDB_FILE = "data/AKT1_ranked_0.pdb"

# ── Load structure ─────────────────────────────────────────────────────────
cmd.load(PDB_FILE, "AKT1")
cmd.remove("solvent")            # remove water molecules
cmd.remove("resn HOH")

# ── RENDER 1: B-factor spectrum ────────────────────────────────────────────
cmd.bg_color("white")
cmd.show("cartoon", "AKT1")
cmd.spectrum("b", "blue_white_red", "AKT1")   # blue=rigid, red=flexible
cmd.set("cartoon_transparency", 0.0)
cmd.orient("AKT1")
cmd.zoom("AKT1", 3)
cmd.ray(1200, 900)
cmd.png(f"{OUTPUT_DIR}/AKT1_bfactor_3D.png", dpi=300)
print("[✓] Render 1 saved → figures/AKT1_bfactor_3D.png")

# ── RENDER 2: Domain-coloured cartoon ─────────────────────────────────────
cmd.color("marine",    "AKT1 and resi 1-107")    # PH domain — blue
cmd.color("tv_orange", "AKT1 and resi 108-152")  # Linker     — orange
cmd.color("forest",    "AKT1 and resi 153-408")  # Kinase     — green
cmd.color("magenta",   "AKT1 and resi 409-480")  # HM domain  — magenta

# Show inhibitor as sticks
cmd.select("inhibitor", "AKT1 and not polymer")
cmd.show("sticks", "inhibitor")
cmd.color("yellow", "inhibitor")
cmd.set("stick_radius", 0.2)

cmd.orient("AKT1")
cmd.zoom("AKT1", 3)
cmd.ray(1200, 900)
cmd.png(f"{OUTPUT_DIR}/AKT1_domains_3D.png", dpi=300)
print("[✓] Render 2 saved → figures/AKT1_domains_3D.png")

# ── RENDER 3: Active-site zoom ─────────────────────────────────────────────
# AKT1 catalytic lysine K179 and DFG motif (D274, F275, G276)
cmd.select("active_site", "AKT1 and resi 179+274+275+276+162+293")
cmd.show("sticks", "active_site")
cmd.color("cyan", "active_site")
cmd.zoom("active_site", 5)
cmd.orient("active_site")
cmd.ray(900, 900)
cmd.png(f"{OUTPUT_DIR}/AKT1_activesite_3D.png", dpi=300)
print("[✓] Render 3 saved → figures/AKT1_activesite_3D.png")

print("\n[✓] All PyMOL renders complete. Check the figures/ directory.")
cmd.quit()
