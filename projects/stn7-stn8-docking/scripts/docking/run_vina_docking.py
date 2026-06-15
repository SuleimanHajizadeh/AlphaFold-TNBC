#!/usr/bin/env python3
"""
run_vina_docking.py
===================
STN7/STN8 — AutoDock Vina Molecular Docking Pipeline

Bu skript bütün liganlar üçün STN7 və STN8 receptorlarına
docking hesablamalarını aparır.

Elmi metodologiya:
  - Bağlama qutusu: ATP-bağlama cəbi ətrafında 20×20×20 Å
  - Exhaustiveness: 16 (yüksək dəqiqlik, 32 CPU üçün optimallaşdırılmış)
  - Num_modes: 10 (ən yaxşı 10 poz)
  - Enerji aralığı: 3 kcal/mol

Çıxış:
  - results/docking/{protein}_{ligand}_poses.pdbqt
  - results/docking/docking_scores.csv
  - results/docking/docking_summary.txt

İstifadə:
  python run_vina_docking.py
  python run_vina_docking.py --protein STN7 --ligand ATP
  python run_vina_docking.py --dry-run
"""

import os
import sys
import argparse
import subprocess
import shutil
import csv
import json
from pathlib import Path
from datetime import datetime

# ── Yollar ──────────────────────────────────────────────────────
BASE_DIR    = Path(__file__).resolve().parent.parent
DATA_DIR    = BASE_DIR / "data" / "structures"
LIG_DIR     = DATA_DIR / "ligands"
RESULTS_DIR = BASE_DIR / "results" / "docking"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── Docking parametrləri ─────────────────────────────────────────
PROTEINS = {
    "STN7": {
        "receptor_pdbqt": DATA_DIR / "STN7_kinase_domain.pdbqt",
        "center": (-5.0, 10.0, -2.0),
        "size": (20.0, 20.0, 20.0),
        "atp_residues": [140, 141, 142, 143, 144, 145, 146, 147, 148, 167],
    },
    "STN8": {
        "receptor_pdbqt": DATA_DIR / "STN8_kinase_domain.pdbqt",
        "center": (-4.5, 11.0, -1.5),
        "size": (20.0, 20.0, 20.0),
        "atp_residues": [139, 140, 141, 142, 143, 144, 145, 146, 147, 186],
    },
}

LIGANDS_AVAILABLE = [
    "ATP", "ADP", "Staurosporine", "K252a",
    "SB216763", "Quercetin", "Imatinib",
]


def get_vina_cmd():
    """Vina icra yolunu tap."""
    if shutil.which("vina"):
        return ["vina"]
    return ["conda", "run", "-n", "bioinfo", "vina"]


def run_single_docking(
    receptor: Path,
    ligand:   Path,
    out_pose: Path,
    center: tuple,
    size: tuple,
    exhaustiveness: int = 16,
    num_modes: int = 10,
) -> dict | None:
    """
    Bir receptor–liganD cütü üçün AutoDock Vina icra et.
    Nəticəni dict olaraq qaytar: {mode, affinity, rmsd_lb, rmsd_ub}
    """
    cx, cy, cz = center
    sx, sy, sz = size
    log_path = out_pose.with_suffix(".log")

    cmd = get_vina_cmd() + [
        "--receptor", str(receptor),
        "--ligand",   str(ligand),
        "--out",      str(out_pose),
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x",   str(sx),
        "--size_y",   str(sy),
        "--size_z",   str(sz),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--energy_range", "3",
        "--cpu", "16",   # 32 CPU-nun yarısı — digər işlər üçün yer qal
    ]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600
        )
        log_path.write_text(result.stdout + result.stderr)

        # Log-dan nəticəni oxu
        return parse_vina_log(result.stdout)
    except subprocess.TimeoutExpired:
        print(f"    ✗ Vaxt aşımı (>10 dəq)")
        return None
    except Exception as e:
        print(f"    ✗ Xəta: {e}")
        return None


def parse_vina_log(log_text: str) -> dict | None:
    """Vina çıxışından ən yaxşı pozun ΔG dəyərini oxu."""
    results = []
    in_table = False
    for line in log_text.split("\n"):
        if "-----+" in line:
            in_table = True
            continue
        if in_table:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    mode     = int(parts[0])
                    affinity = float(parts[1])   # kcal/mol
                    rmsd_lb  = float(parts[2])
                    rmsd_ub  = float(parts[3])
                    results.append({
                        "mode": mode,
                        "affinity_kcal_mol": affinity,
                        "rmsd_lb": rmsd_lb,
                        "rmsd_ub": rmsd_ub,
                    })
                except (ValueError, IndexError):
                    break
    if results:
        best = results[0]
        best["all_modes"] = results
        return best
    return None


def run_all_docking(proteins, ligands, exhaustiveness=16, dry_run=False):
    """Bütün protein–liganD kombinasiyalarını işlət."""
    all_results = []
    timestamp   = datetime.now().strftime("%Y%m%d_%H%M%S")

    total = len(proteins) * len(ligands)
    done  = 0

    for pname in proteins:
        pinfo    = PROTEINS[pname]
        receptor = pinfo["receptor_pdbqt"]

        if not receptor.exists():
            print(f"\n  ⚠️  {pname} receptor tapılmadı: {receptor}")
            print(f"     Əvvəlcə çalışdırın: python prepare_receptor.py --protein {pname}")
            continue

        for lname in ligands:
            ligand_path = LIG_DIR / f"{lname}.pdbqt"
            done += 1

            print(f"\n  [{done}/{total}] {pname} ✕ {lname}")

            if not ligand_path.exists():
                print(f"    ⚠️  Liganl tapılmadı: {ligand_path}")
                print(f"       Əvvəlcə: python prepare_ligands.py --ligand {lname}")
                all_results.append({
                    "protein": pname, "ligand": lname,
                    "affinity_kcal_mol": "N/A", "rmsd_lb": "N/A",
                    "rmsd_ub": "N/A", "status": "LIGAND_MISSING",
                })
                continue

            out_pose = RESULTS_DIR / f"{pname}_{lname}_poses.pdbqt"

            if dry_run:
                print(f"    [dry-run] {pname}_{lname} — icra edilmədi")
                continue

            result = run_single_docking(
                receptor=receptor,
                ligand=ligand_path,
                out_pose=out_pose,
                center=pinfo["center"],
                size=pinfo["size"],
                exhaustiveness=exhaustiveness,
            )

            if result:
                dg = result["affinity_kcal_mol"]
                print(f"    ✅ ΔG = {dg:.2f} kcal/mol  "
                      f"(RMSD lb={result['rmsd_lb']:.2f}, ub={result['rmsd_ub']:.2f})")
                interpret = interpret_affinity(dg)
                print(f"    📊 {interpret}")
                all_results.append({
                    "protein": pname, "ligand": lname,
                    "affinity_kcal_mol": dg,
                    "rmsd_lb": result["rmsd_lb"],
                    "rmsd_ub": result["rmsd_ub"],
                    "status": "OK",
                    "interpretation": interpret,
                })
            else:
                print(f"    ✗ Docking uğursuz")
                all_results.append({
                    "protein": pname, "ligand": lname,
                    "affinity_kcal_mol": "FAILED", "rmsd_lb": "-",
                    "rmsd_ub": "-", "status": "FAILED",
                })

    if not dry_run:
        save_results(all_results, timestamp)

    return all_results


def interpret_affinity(dg: float) -> str:
    """ΔG dəyərini elmi olaraq şərh et."""
    if dg <= -10.0:
        return "Çox güclü bağlanma (nanomolar Ki gözlənilir)"
    elif dg <= -8.0:
        return "Güclü bağlanma (submikromolar Ki)"
    elif dg <= -6.0:
        return "Orta bağlanma (mikromolar Ki)"
    elif dg <= -4.0:
        return "Zəif bağlanma (millimolar Ki)"
    else:
        return "Çox zəif / bağlanma yoxdur"


def save_results(results: list, timestamp: str):
    """Nəticələri CSV və JSON kimi saxla."""
    # CSV
    csv_path = RESULTS_DIR / f"docking_scores_{timestamp}.csv"
    if results:
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\n  ✓ CSV: {csv_path}")

    # JSON
    json_path = RESULTS_DIR / f"docking_scores_{timestamp}.json"
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"  ✓ JSON: {json_path}")

    # Xülasə
    summary_path = RESULTS_DIR / "docking_summary.txt"
    lines = [
        f"STN7/STN8 AutoDock Vina Docking Xülasəsi",
        f"Tarix: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "=" * 55,
        f"{'Protein':<8} {'Liganl':<16} {'ΔG (kcal/mol)':<16} {'Status':<10}",
        "-" * 55,
    ]
    for r in results:
        dg = f"{r['affinity_kcal_mol']:.2f}" if isinstance(r['affinity_kcal_mol'], float) else str(r['affinity_kcal_mol'])
        lines.append(f"{r['protein']:<8} {r['ligand']:<16} {dg:<16} {r['status']:<10}")
    lines.append("=" * 55)
    summary_path.write_text("\n".join(lines))
    print(f"  ✓ Xülasə: {summary_path}")

    # Ən yaxşı nəticəni çap et
    ok_results = [r for r in results if isinstance(r.get("affinity_kcal_mol"), float)]
    if ok_results:
        best = min(ok_results, key=lambda x: x["affinity_kcal_mol"])
        print(f"\n  🏆 Ən güclü bağlanma: {best['protein']} ✕ {best['ligand']}")
        print(f"     ΔG = {best['affinity_kcal_mol']:.2f} kcal/mol")


def main():
    parser = argparse.ArgumentParser(
        description="STN7/STN8 AutoDock Vina Molecular Docking")
    parser.add_argument("--protein", choices=["STN7", "STN8", "both"],
                        default="both")
    parser.add_argument("--ligand", default="all",
                        choices=LIGANDS_AVAILABLE + ["all"])
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--dry-run", action="store_true",
                        help="Yalnız yoxla, docking icra etmə")
    args = parser.parse_args()

    proteins = ["STN7", "STN8"] if args.protein == "both" else [args.protein]
    ligands  = LIGANDS_AVAILABLE if args.ligand == "all" else [args.ligand]

    print("\n" + "=" * 55)
    print("  STN7/STN8 — AutoDock Vina Molecular Docking")
    print(f"  Receptorlar: {proteins}")
    print(f"  Liganlar: {ligands}")
    print(f"  Exhaustiveness: {args.exhaustiveness}")
    print(f"  CPU: 16 (32 nüvəlı serverın yarısı)")
    if args.dry_run:
        print("  *** DRY-RUN rejimi ***")
    print("=" * 55)

    run_all_docking(
        proteins=proteins,
        ligands=ligands,
        exhaustiveness=args.exhaustiveness,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
