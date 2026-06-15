#!/usr/bin/env bash
# ==============================================================================
# STN7/STN8 Computational Structural Bioinformatics Research Pipeline
# Master Orchestrator Script - End-to-End Execution
# ==============================================================================

# ANSI Color Codes for beautiful terminal outputs
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Set strict execution rules
set -e
set -o pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo -e "${BLUE}======================================================================${NC}"
echo -e "${GREEN}  STN7/STN8 Computational Structural Bioinformatics Research Pipeline  ${NC}"
echo -e "${BLUE}======================================================================${NC}"
echo -e "Starting End-to-End Pipeline Execution..."
echo -e "Base working directory: ${BASE_DIR}"

# Step 1: Sequence Retrieval
echo -e "\n${YELLOW}[STAGE 1/8] Running Sequence Retrieval...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/retrieval/retrieve_sequences.py"
echo -e "${GREEN}Stage 1 completed successfully.${NC}"

# Step 2: Structure Downloader
echo -e "\n${YELLOW}[STAGE 2/8] Querying AlphaFold Database and Downloading Structures...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/retrieval/download_structures.py"
echo -e "${GREEN}Stage 2 completed successfully.${NC}"

# Step 3: Multiple Sequence Alignment & Identity Matrix & Heatmap
echo -e "\n${YELLOW}[STAGE 3/8] Running MSA, Calculating Sequence Identity & Generating Conservation Heatmap...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/alignment/run_msa.py"
conda run -n bioinfo python "${BASE_DIR}/scripts/alignment/calculate_identity.py"
conda run -n bioinfo python "${BASE_DIR}/scripts/alignment/conservation_heatmap.py"
echo -e "${GREEN}Stage 3 completed successfully.${NC}"

# Step 4: Structural Analysis
echo -e "\n${YELLOW}[STAGE 4/8] Commencing Structural Superposition, pLDDT & KD Hydrophobicity...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/structure_analysis/structure_compare.py"
conda run -n bioinfo python "${BASE_DIR}/scripts/structure_analysis/plddt_analysis.py"
conda run -n bioinfo python "${BASE_DIR}/scripts/structure_analysis/hydrophobicity.py"
echo -e "${GREEN}Stage 4 completed successfully.${NC}"

# Step 5: Statistics & Mathematics
echo -e "\n${YELLOW}[STAGE 5/8] Executing Quantitative Statistics (Shannon Entropy, t-test, Pearson)...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/statistics/run_stats.py"
echo -e "${GREEN}Stage 5 completed successfully.${NC}"

# Step 6: Evolutionary Phylogeny Tree
echo -e "\n${YELLOW}[STAGE 6/8] Reconstructing Evolutionary Phylogenetic NJ Tree...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/alignment/build_phylogeny.py"
echo -e "${GREEN}Stage 6 completed successfully.${NC}"

# Step 7: Kinase Motif Annotation, AA Composition & Secondary Structure
echo -e "\n${YELLOW}[STAGE 7/9] Annotating Kinase Motifs, AA Composition & Secondary Structure...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/structure_analysis/kinase_motif_annotation.py"
conda run -n bioinfo python "${BASE_DIR}/scripts/structure_analysis/aa_composition.py"
conda run -n bioinfo python "${BASE_DIR}/scripts/structure_analysis/secondary_structure.py"
echo -e "${GREEN}Stage 7 completed successfully.${NC}"

# Step 8: Multi-Panel Publication Summary Figure
echo -e "\n${YELLOW}[STAGE 8/9] Generating Multi-Panel Publication Summary Figure...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/summary_figure.py"
echo -e "${GREEN}Stage 8 completed successfully.${NC}"

# Step 9: Prepare ColabFold Query Inputs
echo -e "\n${YELLOW}[STAGE 9/9] Preparing ColabFold Monomer & Dimer Query FASTA Files...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/prediction/prepare_colabfold_inputs.py"
echo -e "${GREEN}Stage 9 completed successfully.${NC}"

# Step 10: AlphaFold Deep Structural Analysis
echo -e "\n${YELLOW}[STAGE 10/12] AlphaFold Deep Structural Analysis (domain map, motifs, phospho-sites)...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/structure_analysis/alphafold_deep_analysis.py"
echo -e "${GREEN}Stage 10 completed successfully.${NC}"

# Step 11: ColabFold Dimer Pipeline (mock mode — no GPU required)
echo -e "\n${YELLOW}[STAGE 11/12] ColabFold Dimer Workflow (interface contacts, PAE, composition)...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/prediction/run_colabfold_pipeline.py"
echo -e "${GREEN}Stage 11 completed successfully.${NC}"

# Step 12: Cellular Process Analysis (UniProt + STRING + curated)
echo -e "\n${YELLOW}[STAGE 12/12] Intracellular Process Analysis (GO terms, substrates, networks)...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/cellular_processes/cellular_process_analyzer.py"
echo -e "${GREEN}Stage 12 completed successfully.${NC}"

# Step 13: Integrated HTML Report
echo -e "\n${YELLOW}[FINAL] Generating Publication-Quality HTML Report...${NC}"
conda run -n bioinfo python "${BASE_DIR}/scripts/reports/generate_html_report.py"
echo -e "${GREEN}HTML Report generated.${NC}"

echo -e "\n${BLUE}======================================================================${NC}"
echo -e "${GREEN}     Pipeline End-to-End Execution Completed Flawlessly!     ${NC}"
echo -e "${BLUE}======================================================================${NC}"
echo -e ""
echo -e "  Figures    : ${BASE_DIR}/figures/"
echo -e "  Results    : ${BASE_DIR}/results/"
echo -e "  Queries    : ${BASE_DIR}/data/prediction_queries/"
echo -e "  Manuscript : ${BASE_DIR}/manuscript/manuscript_draft.md"
echo -e "  📊 HTML Report: ${BASE_DIR}/results/STN7_STN8_Complete_Analysis_Report.html"
echo -e ""
echo -e "Open the HTML report in your browser for the full publication-quality analysis."
# ==============================================================================
