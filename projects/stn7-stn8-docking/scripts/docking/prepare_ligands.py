#!/usr/bin/env python3
"""
prepare_ligands.py
==================
STN7/STN8 Docking üçün Liganl Hazırlığı

Liganlar:
  1. ATP        — Adenozin trifosfat (təbii substrat)
  2. ADP        — Adenozin difosfat (məhsul)
  3. Staurosporine — Geniş spektrli kinaz inhibitoru (pozitiv kontrol)
  4. K252a      — Staurosporine analoqu, kinaz inhibitoru
  5. SB216763   — GSK-3β inhibitoru (Ser/Thr kinaz hədəfi)
  6. Quercetin  — Bitki polifenolu, kinaz inhibitor aktivliyi bildirilmiş
  7. Imatinib   — Tirozin kinaz inhibitoru (müqayisə üçün)

Məlumat mənbəyi: PubChem REST API (CID → SDF → 3D konformasiya)

İstifadə:
  python prepare_ligands.py
  python prepare_ligands.py --ligand ATP
"""

import os
import sys
import argparse
import requests
import subprocess
import shutil
from pathlib import Path
import time

# ── Liganlar ────────────────────────────────────────────────────
LIGANDS = {
    "ATP": {
        "pubchem_cid": 5957,
        "name": "Adenosine triphosphate",
        "description": "Kinazların təbii substratı — fosfat donor",
        "mw": 507.18,
    },
    "ADP": {
        "pubchem_cid": 644208,
        "name": "Adenosine diphosphate",
        "description": "ATP fosforilləşmə məhsulu",
        "mw": 427.20,
    },
    "Staurosporine": {
        "pubchem_cid": 44259,
        "name": "Staurosporine",
        "description": "Geniş spektrli Ser/Thr kinaz inhibitoru — pozitiv kontrol",
        "mw": 466.53,
    },
    "K252a": {
        "pubchem_cid": 3084981,
        "name": "K252a",
        "description": "Staurosporine analoqu, yüksək seçicilik",
        "mw": 467.52,
    },
    "SB216763": {
        "pubchem_cid": 176158,
        "name": "SB216763",
        "description": "GSK-3 Ser/Thr kinaz inhibitoru",
        "mw": 373.83,
    },
    "Quercetin": {
        "pubchem_cid": 5280343,
        "name": "Quercetin",
        "description": "Bitki flavonoidi — kinaz inhibitor aktivliyi",
        "mw": 302.24,
    },
    "Imatinib": {
        "pubchem_cid": 5291,
        "name": "Imatinib",
        "description": "Gleevec — kinaz inhibitor müqayisəsi üçün",
        "mw": 493.60,
    },
}

PUBCHEM_SDF_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
    "/record/SDF/?record_type=3d&response_type=save"
)

BASE_DIR = Path(__file__).resolve().parent.parent
LIG_DIR  = BASE_DIR / "data" / "structures" / "ligands"
LIG_DIR.mkdir(parents=True, exist_ok=True)


def download_sdf(name: str, cid: int) -> Path | None:
    """PubChem-dən 3D SDF yüklə."""
    sdf_path = LIG_DIR / f"{name}.sdf"
    if sdf_path.exists() and sdf_path.stat().st_size > 100:
        print(f"  [cache] {name}.sdf artıq var")
        return sdf_path

    url = PUBCHEM_SDF_URL.format(cid=cid)
    print(f"  ↓ PubChem CID {cid} → {name}.sdf")
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        sdf_path.write_bytes(r.content)
        time.sleep(0.4)   # PubChem API rate limit
        print(f"  ✓ {name}.sdf ({sdf_path.stat().st_size} bayt)")
        return sdf_path
    except requests.RequestException as e:
        print(f"  ✗ {name} yükləmə xətası: {e}")
        return None


def sdf_to_pdbqt(sdf_path: Path, pdbqt_path: Path) -> bool:
    """
    Open Babel: SDF → PDBQT
    - 3D konformasiya saxlanılır
    - Gasteiger qismli yüklər hesablanır
    - Dönüş zolaqları (rotatable bonds) aşkar edilir
    """
    if not shutil.which("obabel"):
        obabel_cmd = ["conda", "run", "-n", "bioinfo", "obabel"]
    else:
        obabel_cmd = ["obabel"]

    cmd = obabel_cmd + [
        str(sdf_path), "-O", str(pdbqt_path),
        "--gen3d",
        "--addh",
        "--partialcharge", "gasteiger",
        "-p", "7.4",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
            return True
        else:
            print(f"  ✗ {pdbqt_path.name}: {result.stderr[:150]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"  ✗ {sdf_path.name}: vaxt aşımı")
        return False


def process_ligand(name: str, info: dict) -> bool:
    print(f"\n  ── {name} ({info['name']}) ──")
    print(f"     {info['description']}")

    pdbqt_path = LIG_DIR / f"{name}.pdbqt"
    if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
        print(f"  [cache] {name}.pdbqt artıq var")
        return True

    sdf = download_sdf(name, info["pubchem_cid"])
    if sdf is None:
        return False

    ok = sdf_to_pdbqt(sdf, pdbqt_path)
    if ok:
        print(f"  ✓ PDBQT: {pdbqt_path.name}")
    return ok


def write_ligand_summary():
    """Liganlar haqqında xülasə cədvəl yaz."""
    summary_path = LIG_DIR / "ligands_summary.tsv"
    lines = ["Name\tPubChem_CID\tMW\tDescription\tPDBQT_ready\n"]
    for name, info in LIGANDS.items():
        pdbqt = LIG_DIR / f"{name}.pdbqt"
        ready = "YES" if (pdbqt.exists() and pdbqt.stat().st_size > 0) else "NO"
        lines.append(
            f"{name}\t{info['pubchem_cid']}\t{info['mw']}\t"
            f"{info['description']}\t{ready}\n"
        )
    summary_path.write_text("".join(lines))
    print(f"\n  ✓ Xülasə cədvəl: {summary_path}")


def main():
    parser = argparse.ArgumentParser(
        description="STN7/STN8 docking liganlarının hazırlanması")
    parser.add_argument("--ligand", choices=list(LIGANDS.keys()) + ["all"],
                        default="all")
    args = parser.parse_args()

    print("\n" + "=" * 55)
    print("  STN7/STN8 — Liganl Hazırlığı")
    print("  Mənbə: PubChem 3D SDF → OpenBabel PDBQT")
    print("=" * 55)

    targets = list(LIGANDS.keys()) if args.ligand == "all" else [args.ligand]
    results = {}
    for name in targets:
        results[name] = process_ligand(name, LIGANDS[name])

    write_ligand_summary()

    ok   = [n for n, v in results.items() if v]
    fail = [n for n, v in results.items() if not v]
    print(f"\n  ✅ Uğurlu: {ok}")
    if fail:
        print(f"  ✗ Uğursuz: {fail}")
    print(f"  Fayllar: {LIG_DIR}")


if __name__ == "__main__":
    main()
