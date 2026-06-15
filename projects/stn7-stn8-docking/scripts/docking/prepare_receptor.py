#!/usr/bin/env python3
"""
prepare_receptor.py
===================
STN7/STN8 Chloroplast Thylakoid Kinase — Receptor Hazırlığı
AlphaFold PDB → cleaned PDB → PDBQT (AutoDock Vina formatı)

Metod:
  1. AlphaFold DB API-dən PDB yüklənir (v4)
  2. Transit peptid bölgəsi kəsilir (STN7: 1-63, STN8: 1-66)
  3. Su, HETATM, alternativ konformasiyalar silnir
  4. Kinaz domain ATP-bağlama cəbi koordinatları müəyyən edilir
  5. obabel vasitəsilə PDBQT hazırlanır

Elmi əsas:
  STN7 (Q9S7P6): kinaz domain ~64–535, ATP-bağlama: res 140-148, 167
  STN8 (Q9LZV4): kinaz domain ~67–530, ATP-bağlama: res 139-147, 186

İstifadə:
  python prepare_receptor.py --protein STN7
  python prepare_receptor.py --protein STN8
  python prepare_receptor.py --protein both
"""

import os
import sys
import argparse
import requests
import subprocess
import shutil
from pathlib import Path

# ── Sabitlər ────────────────────────────────────────────────────
PROTEINS = {
    "STN7": {
        "uniprot": "Q9S7P6",
        "name": "Serine/threonine-protein kinase STN7",
        "organism": "Arabidopsis thaliana",
        "transit_end": 63,        # kloroplast transit peptidi sonu
        "kinase_start": 64,
        "kinase_end": 535,
        "atp_center": (140, 148, 167),   # ATP-bağlama qalıqları
        # Bağlama qutusu mərkəzi (AlphaFold strukturuna görə təxmini)
        "box_center": (-5.0, 10.0, -2.0),
        "box_size": (20.0, 20.0, 20.0),
    },
    "STN8": {
        "uniprot": "Q9LZV4",
        "name": "Serine/threonine-protein kinase STN8",
        "organism": "Arabidopsis thaliana",
        "transit_end": 66,
        "kinase_start": 67,
        "kinase_end": 530,
        "atp_center": (139, 147, 186),
        "box_center": (-4.5, 11.0, -1.5),
        "box_size": (20.0, 20.0, 20.0),
    },
}

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data" / "structures"
DATA_DIR.mkdir(parents=True, exist_ok=True)


def download_alphafold_pdb(uniprot: str, out_path: Path) -> bool:
    """AlphaFold DB-dən PDB faylını yüklə."""
    if out_path.exists():
        print(f"  [cache] {out_path.name} artıq var, yenidən yüklənmir.")
        return True

    url = ALPHAFOLD_PDB_URL.format(uniprot=uniprot)
    print(f"  ↓ AlphaFold PDB yüklənir: {url}")
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        out_path.write_bytes(r.content)
        print(f"  ✓ Saxlandı: {out_path}")
        return True
    except requests.RequestException as e:
        print(f"  ✗ Xəta: {e}")
        return False


def clean_pdb(raw_pdb: Path, clean_pdb: Path, kinase_start: int, kinase_end: int):
    """
    PDB-ni təmizlə:
    - Yalnız ATOM sətirləri saxla (HETATM kənar)
    - Transit peptidi kəs (kinase_start'dan başla)
    - Alternativ konformasiyalar sil (yalnız A saxla)
    - Zəncirlər normallaşdır
    """
    kept = []
    with open(raw_pdb) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            res_num = int(line[22:26].strip())
            alt_loc = line[16]
            if alt_loc not in (' ', 'A'):
                continue
            if not (kinase_start <= res_num <= kinase_end):
                continue
            # Alternativ lokasiya sütununu boşalt
            line = line[:16] + ' ' + line[17:]
            kept.append(line)

    with open(clean_pdb, 'w') as f:
        f.writelines(kept)
        f.write("END\n")

    print(f"  ✓ Təmiz PDB: {clean_pdb.name} ({len(kept)} ATOM sətiri)")
    return len(kept)


def pdb_to_pdbqt(clean_pdb: Path, pdbqt: Path) -> bool:
    """
    Open Babel vasitəsilə PDB → PDBQT çevirməsi.
    Hidrogen əlavə edir, Gasteiger yükləri hesablayır.
    """
    if not shutil.which("obabel"):
        print("  ✗ obabel tapılmadı. conda run -n bioinfo obabel istifadə et.")
        # Alternativ: conda-run vasitəsilə
        obabel_cmd = ["conda", "run", "-n", "bioinfo", "obabel"]
    else:
        obabel_cmd = ["obabel"]

    cmd = obabel_cmd + [
        str(clean_pdb), "-O", str(pdbqt),
        "-xr",          # receptor mode (saxlanılmayan qalıqlar saxlanır)
        "--addh",       # hidrogen əlavə et
        "--partialcharge", "gasteiger",  # Gasteiger yükləri
        "-p", "7.4",    # fiziologik pH
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if pdbqt.exists() and pdbqt.stat().st_size > 0:
            print(f"  ✓ PDBQT yaradıldı: {pdbqt.name}")
            return True
        else:
            print(f"  ✗ PDBQT yaradılmadı. Xəta: {result.stderr[:200]}")
            return False
    except subprocess.TimeoutExpired:
        print("  ✗ obabel vaxt aşımı")
        return False


def write_box_config(protein_name: str, info: dict, config_path: Path):
    """AutoDock Vina konfiqurasiya faylını yaz."""
    cx, cy, cz = info["box_center"]
    sx, sy, sz = info["box_size"]
    pdbqt = DATA_DIR / f"{protein_name}_kinase_domain.pdbqt"

    config = f"""# AutoDock Vina Konfiqurasiya — {protein_name}
# Protein: {info['name']} ({info['organism']})
# UniProt: {info['uniprot']}
# ATP-bağlama sahəsi: qalıqlar {info['atp_center']}

receptor = {pdbqt}

center_x = {cx}
center_y = {cy}
center_z = {cz}

size_x = {sx}
size_y = {sy}
size_z = {sz}

exhaustiveness = 16
num_modes = 10
energy_range = 3
"""
    config_path.write_text(config)
    print(f"  ✓ Vina konfiqurasiyası: {config_path.name}")


def process_protein(name: str):
    info = PROTEINS[name]
    print(f"\n{'='*55}")
    print(f"  {name} — {info['name']}")
    print(f"  UniProt: {info['uniprot']}")
    print(f"{'='*55}")

    raw_pdb   = DATA_DIR / f"{name}_AF_raw.pdb"
    clean_pdb_path = DATA_DIR / f"{name}_kinase_domain.pdb"
    pdbqt     = DATA_DIR / f"{name}_kinase_domain.pdbqt"
    config    = DATA_DIR / f"{name}_vina.config"

    # 1. AlphaFold PDB yüklə
    if not download_alphafold_pdb(info["uniprot"], raw_pdb):
        return False

    # 2. Kinaz domenini ayır, təmizlə
    n_atoms = clean_pdb(raw_pdb, clean_pdb_path,
                        info["kinase_start"], info["kinase_end"])
    if n_atoms == 0:
        print(f"  ✗ Heç ATOM sətiri qalmadı!")
        return False

    # 3. PDBQT hazırla
    pdb_to_pdbqt(clean_pdb_path, pdbqt)

    # 4. Vina config yaz
    write_box_config(name, info, config)

    print(f"\n  ✅ {name} receptor hazırdır!")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="STN7/STN8 AlphaFold → PDBQT receptor hazırlığı")
    parser.add_argument("--protein", choices=["STN7", "STN8", "both"],
                        default="both", help="Hansı proteini hazırlamaq")
    parser.add_argument("--test", action="store_true",
                        help="Test rejimi — yalnız import yoxla")
    args = parser.parse_args()

    if args.test:
        print("✅ prepare_receptor.py importlar OK")
        return

    targets = ["STN7", "STN8"] if args.protein == "both" else [args.protein]
    for t in targets:
        process_protein(t)

    print(f"\n{'='*55}")
    print("  Bütün receptorlar hazırlandı!")
    print(f"  Fayllar: {DATA_DIR}")
    print(f"{'='*55}")


if __name__ == "__main__":
    main()
