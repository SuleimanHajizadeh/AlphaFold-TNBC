#!/usr/bin/env bash
# =============================================================
# STN7/STN8 Molecular Docking — Ana Pipeline Orkestratorı
# =============================================================
# Conda mühiti: bioinfo
# İstifadə:
#   ./run_docking_pipeline.sh           # tam pipeline
#   ./run_docking_pipeline.sh --demo    # demo məlumatlarla analiz
# =============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV="bioinfo"
DEMO_MODE=false

# Rənglər
GREEN='\033[0;32m'; ORANGE='\033[0;33m'; RED='\033[0;31m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'

log()  { echo -e "${GREEN}[✓]${NC} $*"; }
warn() { echo -e "${ORANGE}[⚠]${NC} $*"; }
err()  { echo -e "${RED}[✗]${NC} $*"; exit 1; }
info() { echo -e "${CYAN}[→]${NC} $*"; }

# Arqumentləri oxu
for arg in "$@"; do
  case $arg in
    --demo) DEMO_MODE=true ;;
    --help)
      echo "İstifadə: $0 [--demo] [--help]"
      echo "  --demo   Əsl docking olmadan nəticə analizini çalışdır"
      exit 0
      ;;
  esac
done

echo ""
echo -e "${BOLD}============================================================${NC}"
echo -e "${BOLD}  STN7/STN8 Chloroplast Thylakoid Kinase Docking Pipeline  ${NC}"
echo -e "${BOLD}  Conda: ${CONDA_ENV}  |  CPU: $(nproc) nüvə               ${NC}"
echo -e "${BOLD}============================================================${NC}"
echo ""

# ── Conda yoxla ──────────────────────────────────────────────────
info "Conda mühiti yoxlanır: ${CONDA_ENV}"
if ! conda env list | grep -q "^${CONDA_ENV}"; then
  err "Conda mühiti '${CONDA_ENV}' tapılmadı!"
fi
log "Conda mühiti mövcuddur"

# ── Asılılıqları yoxla ──────────────────────────────────────────
info "AutoDock Vina yoxlanır..."
if conda run -n "${CONDA_ENV}" vina --version &>/dev/null; then
  log "Vina: $(conda run -n "${CONDA_ENV}" vina --version 2>&1 | head -1)"
else
  warn "Vina tapılmadı, quraşdırılır..."
  conda install -n "${CONDA_ENV}" -c conda-forge vina -y -q
fi

info "Open Babel yoxlanır..."
if conda run -n "${CONDA_ENV}" obabel --version &>/dev/null; then
  log "Open Babel mövcuddur"
else
  warn "Open Babel tapılmadı, quraşdırılır..."
  conda install -n "${CONDA_ENV}" -c conda-forge openbabel -y -q
fi

# Demo rejimdən başqa tam pipeline
if [ "$DEMO_MODE" = false ]; then

  # ── Mərhələ 1: Receptor hazırlığı ───────────────────────────
  echo ""
  echo -e "${BOLD}── Mərhələ 1/3: Receptor Hazırlığı ──────────────────────${NC}"
  info "STN7 və STN8 AlphaFold PDB → PDBQT"
  conda run -n "${CONDA_ENV}" python \
    "${SCRIPT_DIR}/scripts/docking/prepare_receptor.py" --protein both
  log "Receptorlar hazırdır"

  # ── Mərhələ 2: Liganl hazırlığı ────────────────────────────
  echo ""
  echo -e "${BOLD}── Mərhələ 2/3: Liganl Hazırlığı ────────────────────────${NC}"
  info "PubChem SDF → PDBQT (7 liganl)"
  conda run -n "${CONDA_ENV}" python \
    "${SCRIPT_DIR}/scripts/docking/prepare_ligands.py" --ligand all
  log "Liganlar hazırdır"

  # ── Mərhələ 3: Docking ─────────────────────────────────────
  echo ""
  echo -e "${BOLD}── Mərhələ 3/3: AutoDock Vina Docking ───────────────────${NC}"
  info "14 liganl-receptor cütü üçün docking (exhaustiveness=16)"
  warn "Bu proses ~20-60 dəqiqə ala bilər (32 CPU ilə)"
  conda run -n "${CONDA_ENV}" python \
    "${SCRIPT_DIR}/scripts/docking/run_vina_docking.py" \
    --protein both --ligand all --exhaustiveness 16
  log "Docking tamamlandı"

fi

# ── Mərhələ 4: Analiz və Vizuallaşdırma ────────────────────────
echo ""
echo -e "${BOLD}── Nəticə Analizi və Vizuallaşdırma ─────────────────────${NC}"

DEMO_FLAG=""
[ "$DEMO_MODE" = true ] && DEMO_FLAG="--demo"

conda run -n "${CONDA_ENV}" python \
  "${SCRIPT_DIR}/scripts/docking/analyze_docking.py" $DEMO_FLAG
log "Analiz tamamlandı"

# ── Xülasə ─────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}============================================================${NC}"
echo -e "${BOLD}  Pipeline Tamamlandı!                                      ${NC}"
echo -e "${BOLD}============================================================${NC}"
echo ""
log "Receptor faylları: data/structures/"
log "Liganl faylları:   data/structures/ligands/"
log "Docking nəticəsi:  results/docking/"
log "Şəkillər:          figures/"
echo ""
echo -e "${CYAN}Növbəti addım — Notebook-da vizuallaşdırma:${NC}"
echo "  jupyter notebook notebooks/STN7_STN8_AutoDock_Vina_Docking.ipynb"
echo ""
