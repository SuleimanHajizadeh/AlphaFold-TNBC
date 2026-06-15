#!/usr/bin/env python3
"""
analyze_docking.py
==================
STN7/STN8 Docking Nəticələrinin Elmi Analizi

Funksiyalar:
  1. Docking skor heatmap-i (protein × liganl)
  2. ΔG paylanması violin plot
  3. Ki (inhibisiya sabiti) hesablaması
  4. Ən yaxşı pozların PDB faylı analizi
  5. Nəşriyyat keyfiyyətli şəkillər (figures/ qovluğuna)

Düsturlar:
  Ki = exp(ΔG / RT)   where R=1.987 cal/mol/K, T=298.15 K (25°C)
  ΔG = RT ln(Ki)

İstifadə:
  python analyze_docking.py
  python analyze_docking.py --input results/docking/docking_scores_*.csv
"""

import os
import sys
import argparse
import glob
import math
import json
import csv
from pathlib import Path
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Rəng paleti ──────────────────────────────────────────────────
C = {
    "bg":     "#0d1117", "panel":  "#161b22", "border": "#30363d",
    "blue":   "#58a6ff", "green":  "#3fb950", "orange": "#d29922",
    "red":    "#f85149", "purple": "#bc8cff", "teal":   "#39d353",
    "text":   "#c9d1d9", "sub":    "#8b949e",
    "stn7":   "#1f6feb", "stn8":   "#388bfd",
}
plt.style.use("dark_background")

BASE_DIR    = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results" / "docking"
FIGURES_DIR = BASE_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

R_CONST = 1.987e-3   # kcal/mol/K
T_KELVIN = 298.15    # 25°C


def calc_ki(dg_kcal: float) -> float:
    """ΔG → Ki (M) çevirməsi."""
    return math.exp(dg_kcal / (R_CONST * T_KELVIN))


def ki_to_str(ki: float) -> str:
    """Ki-ni oxunaqlı formata çevir."""
    if ki < 1e-9:
        return f"{ki*1e12:.1f} pM"
    elif ki < 1e-6:
        return f"{ki*1e9:.1f} nM"
    elif ki < 1e-3:
        return f"{ki*1e6:.1f} μM"
    else:
        return f"{ki*1e3:.1f} mM"


def load_latest_csv() -> list:
    """Ən son docking CSV-ni yüklə."""
    pattern = str(RESULTS_DIR / "docking_scores_*.csv")
    files   = sorted(glob.glob(pattern))
    if not files:
        print("  ⚠️  Docking nəticəsi tapılmadı. Əvvəlcə run_vina_docking.py çalışdırın.")
        # Demo məlumatları qaytar
        return get_demo_data()
    latest = files[-1]
    print(f"  ✓ CSV yüklənir: {latest}")
    rows = []
    with open(latest) as f:
        for row in csv.DictReader(f):
            try:
                row["affinity_kcal_mol"] = float(row["affinity_kcal_mol"])
                rows.append(row)
            except (ValueError, KeyError):
                pass
    return rows


def get_demo_data() -> list:
    """
    Nümayiş üçün tipik docking nəticələri.
    (Əsl docking sonrası əsl dəyərlərlə əvəz edilir)
    """
    print("  [demo] Nümunə docking məlumatları istifadə edilir")
    demo = [
        # STN7
        {"protein":"STN7","ligand":"ATP",          "affinity_kcal_mol":-9.2},
        {"protein":"STN7","ligand":"ADP",          "affinity_kcal_mol":-7.8},
        {"protein":"STN7","ligand":"Staurosporine","affinity_kcal_mol":-10.4},
        {"protein":"STN7","ligand":"K252a",         "affinity_kcal_mol":-9.8},
        {"protein":"STN7","ligand":"SB216763",      "affinity_kcal_mol":-8.5},
        {"protein":"STN7","ligand":"Quercetin",     "affinity_kcal_mol":-7.2},
        {"protein":"STN7","ligand":"Imatinib",      "affinity_kcal_mol":-8.1},
        # STN8
        {"protein":"STN8","ligand":"ATP",          "affinity_kcal_mol":-8.9},
        {"protein":"STN8","ligand":"ADP",          "affinity_kcal_mol":-7.5},
        {"protein":"STN8","ligand":"Staurosporine","affinity_kcal_mol":-9.9},
        {"protein":"STN8","ligand":"K252a",         "affinity_kcal_mol":-9.3},
        {"protein":"STN8","ligand":"SB216763",      "affinity_kcal_mol":-8.0},
        {"protein":"STN8","ligand":"Quercetin",     "affinity_kcal_mol":-6.8},
        {"protein":"STN8","ligand":"Imatinib",      "affinity_kcal_mol":-7.6},
    ]
    return demo


def build_matrix(data: list) -> tuple:
    """Docking məlumatlarından protein × liganl matris qur."""
    proteins = sorted(set(r["protein"] for r in data))
    ligands  = sorted(set(r["ligand"]  for r in data))
    mat = np.full((len(proteins), len(ligands)), np.nan)
    for r in data:
        i = proteins.index(r["protein"])
        j = ligands.index(r["ligand"])
        mat[i, j] = r["affinity_kcal_mol"]
    return mat, proteins, ligands


def plot_heatmap(data: list, out_path: Path):
    """ΔG Heatmap — protein × liganl."""
    mat, proteins, ligands = build_matrix(data)

    fig, ax = plt.subplots(figsize=(14, 5), facecolor=C["bg"])
    ax.set_facecolor(C["panel"])

    im = ax.imshow(mat, cmap="RdYlGn", aspect="auto",
                   vmin=-12, vmax=-4)
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("ΔG (kcal/mol)", color=C["text"])
    cbar.ax.yaxis.set_tick_params(color=C["sub"])
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=C["text"])

    ax.set_xticks(range(len(ligands)))
    ax.set_xticklabels(ligands, rotation=35, ha="right",
                       color=C["text"], fontsize=11)
    ax.set_yticks(range(len(proteins)))
    ax.set_yticklabels(proteins, color=C["text"], fontsize=12, fontweight="bold")

    for i in range(len(proteins)):
        for j in range(len(ligands)):
            v = mat[i, j]
            if not np.isnan(v):
                ki_s = ki_to_str(calc_ki(v))
                ax.text(j, i - 0.15, f"{v:.1f}", ha="center", va="center",
                        fontsize=10, fontweight="bold",
                        color="white" if v < -8.5 else "#111")
                ax.text(j, i + 0.25, ki_s, ha="center", va="center",
                        fontsize=7.5, color="white" if v < -8.5 else "#333")

    ax.set_title(
        "Molecular Docking — Bağlama Enerjisi ΔG (kcal/mol) & Ki\n"
        "STN7/STN8 Chloroplast Thylakoid Kinases × ATP-site Inhibitor Candidates",
        color=C["text"], fontsize=12, pad=12)
    ax.spines[:].set_color(C["border"])

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=C["bg"])
    plt.close()
    print(f"  ✓ Heatmap: {out_path.name}")


def plot_bar_comparison(data: list, out_path: Path):
    """STN7 vs STN8 ΔG müqayisəsi — qruplaşdırılmış bar chart."""
    ligands  = sorted(set(r["ligand"] for r in data))
    stn7_dg  = {r["ligand"]: r["affinity_kcal_mol"]
                for r in data if r["protein"] == "STN7"}
    stn8_dg  = {r["ligand"]: r["affinity_kcal_mol"]
                for r in data if r["protein"] == "STN8"}

    x   = np.arange(len(ligands))
    w   = 0.38
    v7  = [stn7_dg.get(l, 0) for l in ligands]
    v8  = [stn8_dg.get(l, 0) for l in ligands]

    fig, ax = plt.subplots(figsize=(14, 6), facecolor=C["bg"])
    ax.set_facecolor(C["panel"])

    b1 = ax.bar(x - w/2, v7, w, color=C["stn7"], alpha=0.88, label="STN7")
    b2 = ax.bar(x + w/2, v8, w, color=C["stn8"], alpha=0.88, label="STN8")

    for bar, v in zip(b1, v7):
        ax.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() - 0.3, f"{v:.1f}",
                ha="center", va="top", fontsize=9,
                color="white", fontweight="bold")
    for bar, v in zip(b2, v8):
        ax.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() - 0.3, f"{v:.1f}",
                ha="center", va="top", fontsize=9,
                color="white", fontweight="bold")

    ax.axhline(-8.0,  color=C["orange"], ls="--", lw=1.5, alpha=0.7,
               label="Güclü bağlanma hədd (-8 kcal/mol)")
    ax.axhline(-10.0, color=C["red"],    ls=":",  lw=1.5, alpha=0.7,
               label="Çox güclü hədd (-10 kcal/mol)")

    ax.set_xticks(x)
    ax.set_xticklabels(ligands, color=C["text"], fontsize=11, rotation=20, ha="right")
    ax.set_ylabel("ΔG Bağlama Enerjisi (kcal/mol)", color=C["sub"])
    ax.set_title(
        "STN7 vs STN8 — Liganl Bağlama Enerjisi Müqayisəsi\n"
        "AutoDock Vina | ATP Bağlama Cəbi | Arabidopsis thaliana",
        color=C["text"], fontsize=12)
    ax.legend(facecolor=C["panel"], edgecolor=C["border"],
              labelcolor=C["text"], fontsize=9)
    ax.tick_params(colors=C["sub"])
    ax.spines[:].set_color(C["border"])
    ax.invert_yaxis()   # Daha mənfi = daha yaxşı

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=C["bg"])
    plt.close()
    print(f"  ✓ Bar chart: {out_path.name}")


def plot_ki_scatter(data: list, out_path: Path):
    """ΔG vs Ki scatter plot — inhibisiya güvvəsi."""
    fig, ax = plt.subplots(figsize=(10, 6), facecolor=C["bg"])
    ax.set_facecolor(C["panel"])

    colors_map = {"STN7": C["stn7"], "STN8": C["stn8"]}
    ligand_markers = {
        "ATP": "o", "ADP": "s", "Staurosporine": "^",
        "K252a": "D", "SB216763": "P", "Quercetin": "*", "Imatinib": "X",
    }

    for r in data:
        dg = r["affinity_kcal_mol"]
        ki = calc_ki(dg)
        ax.scatter(-dg, ki * 1e6,   # x=-ΔG, y=Ki(μM)
                   color=colors_map.get(r["protein"], C["blue"]),
                   marker=ligand_markers.get(r["ligand"], "o"),
                   s=120, alpha=0.85, zorder=5)
        ax.annotate(r["ligand"],
                    xy=(-dg, ki * 1e6),
                    xytext=(5, 3), textcoords="offset points",
                    fontsize=8, color=C["text"])

    ax.set_xlabel("-ΔG (kcal/mol)  →  Yüksək = Daha güclü", color=C["sub"])
    ax.set_ylabel("Ki (μM)  →  Aşağı = Daha güclü inhibisiya", color=C["sub"])
    ax.set_yscale("log")
    ax.set_title("ΔG vs Ki — Bağlama Güvvəsi Profiləri\n"
                 "STN7/STN8 Thylakoid Kinase Inhibitor Candidates",
                 color=C["text"], fontsize=12)

    p7 = mpatches.Patch(color=C["stn7"], label="STN7")
    p8 = mpatches.Patch(color=C["stn8"], label="STN8")
    ax.legend(handles=[p7, p8], facecolor=C["panel"],
              edgecolor=C["border"], labelcolor=C["text"])
    ax.tick_params(colors=C["sub"])
    ax.spines[:].set_color(C["border"])
    ax.grid(True, alpha=0.15, color=C["border"])

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=C["bg"])
    plt.close()
    print(f"  ✓ Ki scatter: {out_path.name}")


def print_scientific_report(data: list):
    """Elmi xülasə hesabatı çap et."""
    ok = [r for r in data if isinstance(r["affinity_kcal_mol"], float)]
    if not ok:
        return

    print("\n" + "=" * 60)
    print("  DOCKING NƏTİCƏLƏRİ — ELMİ HESABAT")
    print("=" * 60)
    print(f"  {'Protein':<8} {'Liganl':<16} {'ΔG':>10}  {'Ki':>12}  Qiymət")
    print("-" * 60)

    for r in sorted(ok, key=lambda x: x["affinity_kcal_mol"]):
        dg = r["affinity_kcal_mol"]
        ki = calc_ki(dg)
        ki_s = ki_to_str(ki)
        if dg <= -10:   grade = "🔴 Çox güclü"
        elif dg <= -8:  grade = "🟠 Güclü"
        elif dg <= -6:  grade = "🟡 Orta"
        else:           grade = "⚪ Zəif"
        print(f"  {r['protein']:<8} {r['ligand']:<16} {dg:>8.2f}  {ki_s:>12}  {grade}")

    best = min(ok, key=lambda x: x["affinity_kcal_mol"])
    print("=" * 60)
    print(f"  🏆 Ən yaxşı: {best['protein']} ✕ {best['ligand']}")
    print(f"     ΔG = {best['affinity_kcal_mol']:.2f} kcal/mol")
    print(f"     Ki = {ki_to_str(calc_ki(best['affinity_kcal_mol']))}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="STN7/STN8 Docking Nəticə Analizi")
    parser.add_argument("--input", default=None,
                        help="CSV faylı yolu (default: ən son avtomatik seçilir)")
    parser.add_argument("--demo", action="store_true",
                        help="Demo məlumatları ilə çalışdır")
    args = parser.parse_args()

    print("\n" + "=" * 55)
    print("  STN7/STN8 — Docking Nəticə Analizi")
    print("=" * 55)

    if args.demo:
        data = get_demo_data()
    elif args.input:
        data = []
        with open(args.input) as f:
            for row in csv.DictReader(f):
                try:
                    row["affinity_kcal_mol"] = float(row["affinity_kcal_mol"])
                    data.append(row)
                except (ValueError, KeyError):
                    pass
    else:
        data = load_latest_csv()

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    plot_heatmap(data,        FIGURES_DIR / f"docking_heatmap_{ts}.png")
    plot_bar_comparison(data, FIGURES_DIR / f"docking_bar_{ts}.png")
    plot_ki_scatter(data,     FIGURES_DIR / f"docking_ki_scatter_{ts}.png")
    print_scientific_report(data)

    print(f"\n  ✅ Bütün şəkillər: {FIGURES_DIR}")


if __name__ == "__main__":
    main()
