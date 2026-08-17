---
title: "Wiki Log"
type: log
created: 2026-08-17
updated: 2026-08-17
---

# Wiki Log

Chronological audit ledger of all knowledge base operations (ingestion, major syntheses, and lint passes).

---

## [2026-08-17] scaffold | LLM Wiki Architecture Initialized
- **Action**: Initialized wiki file structure, index catalog, and `AGENTS.md` governance rules.
- **Layers**:
  - `raw/` configured for immutable documents (`papers/`, `articles/`, `transcripts/`, `assets/`).
  - `wiki/` configured for persistent knowledge (`entities/`, `concepts/`, `sources/`, `syntheses/`).
- **Status**: Ready for initial source ingestion.

---

## [2026-08-17] ingest | Computational Structural Biology Codebase & Projects
- **Ingested Projects**:
  - `projects/akt1-kinase-modeling/`
  - `projects/stn7-stn8-docking/`
  - `bioinformatics_starter_kit/`
- **Created Source Notes**:
  - `wiki/sources/Source - AKT1 Kinase Modeling Project.md`
  - `wiki/sources/Source - STN7-STN8 Docking Project.md`
  - `wiki/sources/Source - Bioinformatics Starter Kit.md`
- **Created Entity Notes**:
  - `wiki/entities/AKT1.md`, `wiki/entities/STN7.md`, `wiki/entities/STN8.md`, `wiki/entities/AlphaFold2.md`, `wiki/entities/TP53.md`
- **Created Concept Notes**:
  - `wiki/concepts/Local Distance Difference Test (LDDT & pLDDT).md`
  - `wiki/concepts/Kinase Activation Loop & Allostery.md`
  - `wiki/concepts/Ramachandran Dihedral Angles.md`
  - `wiki/concepts/Contact Map Analysis.md`
  - `wiki/concepts/Thylakoid State Transitions.md`
  - `wiki/concepts/Structural Superposition & RMSD.md`
- **Created Syntheses**:
  - `wiki/syntheses/Comparative Kinase Mechanics - AKT1 vs STN7-STN8.md`
- **Catalog Update**: Updated `wiki/index.md` with full categorized links and cross-references.
- **Key Takeaway**: Compiled repository codebases into an interlinked knowledge graph spanning cancer hub kinase dynamics, chloroplast thylakoid state transitions, and structural bioinformatics algorithms.

---

## [2026-08-17] ingest | Structural Bioinformatics Course Curriculum
- **Ingested Directory**: `Structural Biology Tutorials/Structural-Bioinformatics/`
- **Created Source Notes**:
  - `wiki/sources/Source - Structural Bioinformatics Course Curriculum.md`
- **Created Concept Notes**:
  - `wiki/concepts/Molecular Dynamics & Enhanced Sampling.md`
  - `wiki/concepts/Molecular Docking Principles.md`
  - `wiki/concepts/Markov State Models.md`
  - `wiki/concepts/De Novo Protein Design.md`
- **Cross-Updated**:
  - Linked new concepts to `[[AlphaFold2]]`, `[[AKT1]]`, `[[STN7]]`, `[[STN8]]`, and updated `wiki/index.md`.
- **Key Takeaway**: Integrated foundational physical and machine-learning frameworks (force fields, metadynamics, MSMs, Lamarckian docking, RFdiffusion, and ProteinMPNN) into the compiled knowledge base.

---

## [2026-08-17] synthesis | Kinase Signaling, Photosynthetic Antennas, and Allosteric Drug Discovery
- **Created Entity Notes**:
  - `wiki/entities/LHCII.md` (Light-Harvesting Complex II trimer)
  - `wiki/entities/Cytochrome b6f.md` (Integral thylakoid redox sensor)
- **Created Concept Notes**:
  - `wiki/concepts/PI3K-Akt-mTOR Signaling Pathway.md` (Signaling hierarchy & oncogenic dysregulation)
- **Created Synthesis Note**:
  - `wiki/syntheses/Allosteric Drug Discovery & Conformational Selection in Kinases.md`
- **Catalog Update**: Updated `wiki/index.md` with new entities, concepts, and syntheses.
- **Key Takeaway**: Bridged molecular signaling networks with physical conformational sampling (AF2 disorder, enhanced MD, MSMs, and ensemble docking) for allosteric kinase targeting.

---

## [2026-08-17] ingest & synthesis | Biophysical Principles, Graph Theory, and Multi-Scale Paradigms
- **Ingested Tutorials & Literature**:
  - `Structural Biology Tutorials/Principles_of_nucleic_acid_structure.pdf`
  - `Structural Biology Tutorials/graph_theory.pdf`
  - `Structural Biology Tutorials/SECOND_EDITION_Jonathan_Clayden_Nick_Greeves_Stuart_Warren_PDFDrive.pdf`
- **Created Source Notes**:
  - `wiki/sources/Source - Principles of Nucleic Acid Structure.md`
  - `wiki/sources/Source - Graph Theory and Network Biology.md`
  - `wiki/sources/Source - Chemical Principles of Enzyme Catalysis & Phosphoryl Transfer.md`
- **Created Concept Notes**:
  - `wiki/concepts/DNA-Protein Recognition & Binding Dynamics.md`
  - `wiki/concepts/Residue Interaction Networks & Graph Centrality.md`
  - `wiki/concepts/Biochemical Mechanism of Phosphoryl Transfer.md`
- **Created Synthesis Note**:
  - `wiki/syntheses/Network to Structure - From WGCNA Hubs to Allosteric Kinase Inhibitors.md`
- **Catalog Update**: Fully updated `wiki/index.md`.
- **Lint Audit**: 0 broken links, 0 orphan notes, 100% interconnected knowledge graph.
- **Key Takeaway**: Unified multi-scale biology from transcriptomic network hubs (WGCNA) to atomic kinase mechanics (AF2, RINs, MSMs, and allosteric drug design).

---

## [2026-08-17] synthesis | Experimental Structural Methods, Kinase Regulators, and Thylakoid Circuits
- **Created Entity Notes**:
  - `wiki/entities/PDK1.md` (Master AGC kinase activator of AKT1 Thr308)
  - `wiki/entities/mTORC2.md` (Hydrophobic motif kinase of AKT1 Ser473)
  - `wiki/entities/PTEN.md` (Tumor suppressor lipid phosphatase)
  - `wiki/entities/TAP38-PPH1.md` (Thylakoid phosphatase antagonizing STN7)
  - `wiki/entities/PSII Reaction Center D1.md` (Core reaction center protein phosphorylated by STN8)
- **Created Concept Notes**:
  - `wiki/concepts/Cryo-Electron Microscopy (Cryo-EM).md` (Single-particle SPA & membrane complex solving)
  - `wiki/concepts/X-Ray Crystallography & Electron Density.md` (Diffraction, phase problem & B-factors)
  - `wiki/concepts/Intrinsically Disordered Proteins & Condensates.md` (Ensembles, LLPS & AF2 disorder profiling)
- **Created Synthesis Notes**:
  - `wiki/syntheses/Structural Biology Methods - Experimental Cryo-EM and X-Ray vs Deep Learning AF2.md`
  - `wiki/syntheses/Thylakoid Phosphorylation Dynamics - STN7-TAP38 Circuit vs STN8 Repair.md`
- **Catalog Update**: Fully updated `wiki/index.md` (now containing 44 total interlinked pages).
- **Key Takeaway**: Comprehensive mapping of kinase regulatory networks, biophysical experimental standards, and plant photosynthetic photoprotection circuitry.

---

## [2026-08-17] ingest & refinement | STN7-STN8 Manuscript Findings & AutoDock Vina Pipeline
- **Ingested Manuscript & Code**:
  - `projects/stn7-stn8-docking/manuscript/manuscript_draft.md`
  - `projects/stn7-stn8-docking/scripts/docking/` (`prepare_receptor.py`, `run_vina_docking.py`, `analyze_docking.py`)
- **Updated Notes**:
  - `[[STN7]]` & `[[STN8]]`: Enriched with exact motif signatures (non-canonical `DLG` vs `DFG`, `APE` vs `PPE`, `GEGSFG`, `HRD`), domain pLDDT $t$-test values ($p < 10^{-40}$), secondary structure fractions, and ColabFold multimer heterodimer contact interface (2,190 contacts).
  - `[[Kinase Activation Loop & Allostery]]`: Added specific plant non-canonical motif variations.
  - `[[Source - STN7-STN8 Docking Project]]`: Added full manuscript conclusions and virtual screening workflow.
- **Key Takeaway**: Integrated empirical manuscript discoveries and automated virtual screening pipeline into core entity and source files.

---

## [2026-08-17] documentation | Obsidian & LLM Wiki Integration Guide
- **Created Synthesis**:
  - `wiki/syntheses/Obsidian & LLM Wiki Integration Guide.md`
- **Catalog Update**: Added guide to `wiki/index.md`.
- **Key Takeaway**: Comprehensive manual documenting the three-layer architecture, file taxonomy, graph visual navigation, and compounding workflows.
