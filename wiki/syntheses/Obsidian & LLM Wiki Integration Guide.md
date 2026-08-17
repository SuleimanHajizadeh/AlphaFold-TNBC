---
title: "Obsidian & LLM Wiki Integration Guide"
type: synthesis
tags:
  - synthesis/wiki-architecture
  - obsidian/guide
  - knowledge-management
created: 2026-08-17
updated: 2026-08-17
aliases:
  - "Obsidian LLM Wiki Guide"
  - "LLM Wiki Documentation"
  - "Wiki Manual"
---

# Obsidian & LLM Wiki Integration Guide: The Codebase Model of Knowledge

This manual documents the entire architecture, file taxonomy, operational workflows, and practical user interface techniques for navigating and compounding knowledge using **Obsidian as the IDE** and the **LLM Agent as the Compiler/Maintainer**.

---

## 🏛️ 1. The Three-Layer Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│  LAYER 1: RAW SOURCES (Immutable Ground Truth — Human-Curated)          │
│  raw/papers/   raw/articles/   raw/transcripts/   raw/assets/           │
└────────────────────────────────────┬────────────────────────────────────┘
                                     │
                                     ▼ Ingestion & Extraction (Agent)
┌─────────────────────────────────────────────────────────────────────────┐
│  LAYER 2: PERSISTENT WIKI (Compiling Knowledge Layer — Agent-Maintained)│
│  wiki/                                                                  │
│  ├── index.md             (Central Categorized Content Catalog)        │
│  ├── log.md               (Chronological Append-Only Audit Ledger)     │
│  ├── entities/            (Proteins, Genes, Complexes, AI Systems)     │
│  ├── concepts/            (Mechanisms, Biophysics, ML Architectures)   │
│  ├── sources/             (Structured Summaries of Raw Ingestions)     │
│  └── syntheses/           (Deep-Dive Cross-Cutting Compilations)       │
└────────────────────────────────────▲────────────────────────────────────┘
                                     │
                                     ▼ Governance & Conventions
┌─────────────────────────────────────────────────────────────────────────┐
│  LAYER 3: SCHEMA & RULES (AGENTS.md / LLM Instructions)                 │
│  Standardized YAML Frontmatter, [[Wikilinks]], Ingest/Query/Lint Rules  │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 📁 2. Complete File Taxonomy & Directory Map

### 2.1 Catalog & Audit Ledger
* **`wiki/index.md`**: Master table of contents. Organized into Entities, Concepts, Syntheses, and Sources. Every newly generated note is registered here.
* **`wiki/log.md`**: Chronological audit trail using parseable headings (`## [YYYY-MM-DD] <action> | <Subject>`). Records all ingests, syntheses, and lint passes.

### 2.2 Core Knowledge Pages
* **`wiki/entities/`** (12 active pages):
  * Human Oncology: [[AKT1]], [[PDK1]], [[mTORC2]], [[PTEN]], [[TP53]].
  * Photosynthesis & Thylakoid: [[STN7]], [[STN8]], [[LHCII]], [[Cytochrome b6f]], [[TAP38-PPH1]], [[PSII Reaction Center D1]].
  * Deep Learning Systems: [[AlphaFold2]].
* **`wiki/concepts/`** (14 active pages):
  * Enzymatic & Signaling Mechanisms: [[Kinase Activation Loop & Allostery]], [[PI3K-Akt-mTOR Signaling Pathway]], [[Biochemical Mechanism of Phosphoryl Transfer]], [[Thylakoid State Transitions]].
  * Biophysical & Network Principles: [[DNA-Protein Recognition & Binding Dynamics]], [[Residue Interaction Networks & Graph Centrality]], [[Intrinsically Disordered Proteins & Condensates]].
  * Experimental Methods: [[Cryo-Electron Microscopy (Cryo-EM)]], [[X-Ray Crystallography & Electron Density]].
  * Computational Structural Metrics: [[Local Distance Difference Test (LDDT & pLDDT)]], [[Ramachandran Dihedral Angles]], [[Contact Map Analysis]], [[Structural Superposition & RMSD]].
  * Machine Learning & Simulations: [[Molecular Dynamics & Enhanced Sampling]], [[Molecular Docking Principles]], [[Markov State Models]], [[De Novo Protein Design]].
* **`wiki/sources/`** (7 active pages):
  * Summaries for internal projects and external literature: [[Source - AKT1 Kinase Modeling Project]], [[Source - STN7-STN8 Docking Project]], [[Source - Bioinformatics Starter Kit]], [[Source - Structural Bioinformatics Course Curriculum]], [[Source - Principles of Nucleic Acid Structure]], [[Source - Graph Theory and Network Biology]], [[Source - Chemical Principles of Enzyme Catalysis & Phosphoryl Transfer]].
* **`wiki/syntheses/`** (6 active pages):
  * Multi-scale cross-cutting reports: [[Comparative Kinase Mechanics - AKT1 vs STN7-STN8]], [[Allosteric Drug Discovery & Conformational Selection in Kinases]], [[Network to Structure - From WGCNA Hubs to Allosteric Kinase Inhibitors]], [[Structural Biology Methods - Experimental Cryo-EM and X-Ray vs Deep Learning AF2]], [[Thylakoid Phosphorylation Dynamics - STN7-TAP38 Circuit vs STN8 Repair]], and this integration guide.

---

## 🖥️ 3. How to Use Obsidian With the LLM Wiki

Obsidian is your interactive browser, visual workbench, and thought-canvas:

```
          ┌─────────────────────────────────────────────────────────┐
          │                    OBSIDIAN WORKBENCH                   │
          │                                                         │
          │   [ Left Sidebar ]      [ Center Pane ]   [ Right Pane ]│
          │   • File Explorer       • Active Note     • Backlinks   │
          │   • Search              • Live Preview    • Outgoing    │
          │   • Tag Pane            • MathJax (LaTeX) • Local Graph │
          │                                                         │
          │                   [ Ctrl/Cmd + G ]                      │
          │                   Global Graph View                     │
          └─────────────────────────────────────────────────────────┘
```

### 3.1 Essential Obsidian Views & Features
1. **Interactive Graph View (`Cmd/Ctrl + G`)**:
   * Open Graph View to visualize the entire knowledge topology.
   * Hub notes (like [[AKT1]], [[STN7]], [[AlphaFold2]]) emerge naturally with dense link clusters.
   * Color-code nodes in Graph settings using frontmatter tags (e.g. `type: entity`, `type: concept`, `type: synthesis`).
2. **Local Graph View (Right Sidebar)**:
   * Displays the 1-hop and 2-hop neighborhood of the active note you are reading.
3. **Backlinks Pane (`Linked Mentions`)**:
   * When reading any entity (e.g. [[STN8]]), the Backlinks pane automatically reveals every concept, source, and synthesis referencing it.
4. **Obsidian Dataview Queries**:
   * Every note includes YAML frontmatter (`type`, `tags`, `created`, `updated`). You can write dynamic Dataview blocks directly in notes:
   ```markdown
   ```dataview
   TABLE updated as "Last Updated", tags as "Tags"
   FROM "wiki/entities"
   SORT updated desc
   ```
   ```

---

## ⚡ 4. The Human-Agent Operational Workflows

### 4.1 The `Ingest` Workflow
```
Drop PDF / Note into raw/ ──► Ask Agent to Ingest ──► Agent Updates Sources, Entities, Concepts, Index & Log
```
* **Your role**: Add new raw sources (`raw/papers/`, web clippings, or discussion transcripts).
* **Agent role**: Reads the document, extracts key quantitative parameters, creates `wiki/sources/<source>.md`, updates relevant entity and concept pages with cross-references, updates `wiki/index.md`, and logs the operation.

### 4.2 The `Query` Workflow
```
Ask Cross-Cutting Question ──► Agent Traverses Wiki ──► Synthesizes Answer with [[Wikilinks]] ──► Persists in syntheses/
```
* Instead of losing valuable research answers in transient chat logs, the agent compiles deep answers into persistent notes in `wiki/syntheses/`. Knowledge accumulates permanently.

### 4.3 The `Lint` Workflow
```
Request Wiki Lint ──► Agent Scans All Links & Indexes ──► Resolves Orphans/Gaps ──► Commits to Git
```
* The agent executes automated health-checks to ensure 0 broken links, 0 orphaned notes, and consistency across all definitions.

---

## 🔗 5. Summary
The LLM Wiki eliminates the classic maintenance decay of personal wikis. By delegating file management, cross-referencing, mathematical documentation, and logging to the LLM agent, you maintain an ever-growing, publication-ready second brain in computational structural biology.
