# LLM Wiki Schema & Agent Governance

This document defines the architecture, conventions, and operational workflows for maintaining the LLM Wiki in this repository.

---

## 1. Directory Structure

```
.
├── raw/                      # Immutable source documents (Human-curated)
│   ├── papers/               # PDF, markdown, or text papers & preprints
│   ├── articles/             # Clipped web pages, articles, notes
│   ├── transcripts/          # Meeting or lecture notes, discussions
│   └── assets/               # Local images and figures referenced by sources
│
├── wiki/                     # Persistent, compiled knowledge base (Agent-maintained)
│   ├── index.md              # Catalog of all wiki notes categorized by type
│   ├── log.md                # Chronological append-only record of all actions
│   ├── entities/             # Specific entities (e.g. [[AKT1]], [[STN7]], [[AlphaFold3]])
│   ├── concepts/             # Broad concepts (e.g. [[Kinase Activation Loop]], [[Molecular Docking]])
│   ├── sources/              # Structured notes for each ingested raw document
│   └── syntheses/            # Deep-dive analyses, comparisons, and answered queries
│
└── AGENTS.md                 # This specification & operational instructions
```

---

## 2. Note Conventions & Formatting

### 2.1 Frontmatter
Every wiki note must begin with YAML frontmatter for compatibility with Obsidian Dataview:

```yaml
---
title: "Page Title"
type: entity | concept | source | synthesis
tags:
  - domain/subdomain
created: YYYY-MM-DD
updated: YYYY-MM-DD
aliases:
  - "Alias 1"
  - "Alias 2"
---
```

### 2.2 Cross-Linking & Synthesis
* Use standard Obsidian wikilinks: `[[Target Note Name]]` or `[[Target Note Name|Display Text]]`.
* When writing or updating a note, proactively link to existing concepts and entities.
* Do not link generic or common English words; only link meaningful domain terms, entities, and methods.
* Highlight contradictions, unresolved questions, or contrasting evidence when new sources challenge existing notes.

---

## 3. Core Operational Workflows

### 3.1 `Ingest` Workflow
When the user asks to ingest a new source (or drops a file into `raw/`):
1. **Read & Extract**: Parse the source material, extracting key claims, methodology, quantitative results, and citations.
2. **Create Source Note**: Create `wiki/sources/<source-identifier>.md` with frontmatter, executive summary, key findings, and references.
3. **Cross-Update Entity & Concept Pages**:
   - Check existing notes in `wiki/entities/` and `wiki/concepts/`.
   - Update existing notes with new insights, updating the `updated` timestamp.
   - Create new entity/concept notes if a core concept is missing.
4. **Update Catalog**: Add the new pages to `wiki/index.md` under their appropriate categories.
5. **Log Entry**: Append an entry to `wiki/log.md` with format:
   ```markdown
   ## [YYYY-MM-DD] ingest | <Source Title>
   - Ingested: `raw/...`
   - Created: `wiki/sources/...`
   - Updated: `[[Entity 1]]`, `[[Concept 2]]`
   - Key Takeaway: One sentence summary.
   ```

### 3.2 `Query` Workflow
When the user asks a question against the knowledge base:
1. **Catalog Scan**: Consult `wiki/index.md` to identify relevant candidate notes.
2. **Deep Synthesis**: Read the specific wiki notes (and raw sources if primary evidence is needed).
3. **Formulate Answer**: Provide a coherent synthesis with inline `[[wikilinks]]` to the supporting notes.
4. **Persist Valuable Syntheses**: If the answer is an in-depth comparison, architectural overview, or novel insight, offer or automatically create a persistent note in `wiki/syntheses/<topic>.md` and link it in `wiki/index.md` and `wiki/log.md`.

### 3.3 `Lint` Workflow
When asked to health-check or lint the wiki:
1. Scan all markdown files in `wiki/`.
2. Check for:
   - **Broken Links**: Wikilinks pointing to non-existent notes.
   - **Orphan Pages**: Notes with 0 incoming backlinks.
   - **Missing Key Pages**: Terms linked frequently that lack their own dedicated page.
   - **Stale / Conflicting Claims**: Divergent statements across different summaries.
3. Provide a structured remediation report and propose batch fixes.

---

## 4. Log Format
Every mutating operation (ingest, major synthesis, lint pass) MUST append to `wiki/log.md` using the parseable heading pattern:
`## [YYYY-MM-DD] <action> | <Subject>`
