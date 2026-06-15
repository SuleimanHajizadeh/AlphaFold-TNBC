#!/usr/bin/env python3
"""
STN7/STN8 Intracellular Process Analyzer
Fetches UniProt + STRING data, builds comprehensive cellular process catalogue.
"""
import os, sys, json, time
import urllib.request, urllib.error
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE_DIR    = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

PROTEINS = {
    "STN7": {"uniprot": "Q9S713", "tax_id": "3702", "string_id": "3702.AT1G68830.1"},
    "STN8": {"uniprot": "Q9LZV4", "tax_id": "3702", "string_id": "3702.AT5G01920.1"},
}

# ── Curated Literature Data (authoritative fallback) ──────────────────────
CURATED = {
    "STN7": {
        "full_name": "Serine/threonine-protein kinase STN7, chloroplastic",
        "gene_name": "STN7 / AT1G68830",
        "organism":  "Arabidopsis thaliana",
        "subcellular_locations": [
            "Chloroplast thylakoid membrane (stromal face)",
            "Thylakoid lumen (transit peptide cleaved)",
        ],
        "molecular_functions": [
            "Serine/threonine kinase activity (EC 2.7.11.1)",
            "ATP binding (P-loop GxGxxG motif)",
            "Protein substrate recognition (LHCII proteins)",
            "Redox-sensing via disulfide bond (Cys residues in lumen domain)",
        ],
        "biological_processes": [
            "State transition: State 1 → State 2 (light-harvesting regulation)",
            "Phosphorylation of LHCII proteins (LHCB1, LHCB2, LHCB3)",
            "Regulation of photosynthetic electron flow (PSI/PSII balancing)",
            "Chloroplast retrograde signaling (kinase-mediated)",
            "Plant acclimation to fluctuating light conditions",
            "Regulation of cyclic electron flow around PSI",
        ],
        "substrates": [
            {"name": "LHCB1.1", "uniprot": "P04778", "site": "Thr-5",  "effect": "LHCII migration from PSII to PSI"},
            {"name": "LHCB1.2", "uniprot": "P0CJ48", "site": "Thr-5",  "effect": "State 1→2 transition"},
            {"name": "LHCB2.1", "uniprot": "P04779", "site": "Thr-5",  "effect": "Grana unstacking"},
            {"name": "LHCB3",   "uniprot": "Q9SYW8", "site": "Thr-5",  "effect": "Antenna remodeling"},
            {"name": "STN7",    "uniprot": "Q9S713", "site": "Thr-309", "effect": "Autophosphorylation (activation loop)"},
        ],
        "regulatory_inputs": [
            "Activated by: reduced plastoquinone pool (high PSII light)",
            "Inhibited by: oxidized PQ pool (high PSI light)",
            "Inhibited by: thioredoxin (Trx-f) in stroma — reduces disulfide in lumen domain",
            "Inhibited by: TAP38/PPH1 phosphatase (dephosphorylates LHCII)",
            "Regulated by: PSI-H and PSI-L subunits (kinase activation scaffold)",
        ],
        "interaction_partners": [
            {"name": "TAP38/PPH1",  "role": "Antagonistic phosphatase — reverses STN7 action"},
            {"name": "PSI-H (PSAH)","role": "Required for STN7 kinase activation"},
            {"name": "PSI-L (PSAL)","role": "Required for LHCII–PSI docking"},
            {"name": "LHCB1",       "role": "Primary substrate — LHCII major component"},
            {"name": "LHCB2",       "role": "Primary substrate — LHCII major component"},
            {"name": "Cytochrome b6f","role": "Activation signal source (PQ oxidation sensor)"},
            {"name": "STN8",        "role": "Functional paralog — shared regulatory network"},
        ],
        "go_terms": {
            "molecular_function": ["GO:0004674 (protein serine/threonine kinase)","GO:0005524 (ATP binding)","GO:0046983 (protein dimerization)"],
            "biological_process": ["GO:0009651 (state transition)","GO:0018105 (peptidyl-serine phosphorylation)","GO:0009768 (photosynthesis, light harvesting)","GO:0010218 (response to far red light)","GO:0048437 (floral organ development)"],
            "cellular_component": ["GO:0009535 (chloroplast thylakoid membrane)","GO:0009534 (chloroplast thylakoid)","GO:0009570 (chloroplast stroma)"],
        },
        "cellular_pathway": [
            "1. High PSII illumination → PQ pool reduction",
            "2. Cyt b6f complex senses reduced PQ → activates STN7",
            "3. STN7 kinase domain phosphorylates LHCII Thr-5",
            "4. Phospho-LHCII detaches from PSII (grana)",
            "5. Phospho-LHCII migrates to stroma lamellae",
            "6. Phospho-LHCII docks at PSI via PSI-H/PSI-L",
            "7. PSI antenna size increases → balanced excitation",
            "8. Reversed by TAP38/PPH1 phosphatase under State 1 conditions",
        ],
    },
    "STN8": {
        "full_name": "Serine/threonine-protein kinase STN8, chloroplastic",
        "gene_name": "STN8 / AT5G01920",
        "organism":  "Arabidopsis thaliana",
        "subcellular_locations": [
            "Chloroplast thylakoid membrane (stromal face)",
            "Appressed thylakoid membrane (grana)",
        ],
        "molecular_functions": [
            "Serine/threonine kinase activity (EC 2.7.11.1)",
            "ATP binding (P-loop GxGxxG motif)",
            "PSII core protein substrate recognition",
            "Regulation of PSII repair cycle via D1 phosphorylation",
        ],
        "biological_processes": [
            "Phosphorylation of PSII core proteins (D1, D2, CP43, PsbH)",
            "Regulation of PSII repair cycle under high-light stress",
            "Photoinhibition recovery (D1 turnover facilitation)",
            "Grana margin dynamics (PSII monomer release)",
            "Chloroplast quality control (damaged PSII targeting)",
            "High-light acclimation response",
        ],
        "substrates": [
            {"name": "D1 (PsbA)",  "uniprot": "P83755", "site": "Thr-2 (N-term)", "effect": "D1 targeting to FtsH protease degradation"},
            {"name": "D2 (PsbD)",  "uniprot": "P56762", "site": "Ser-2 (N-term)", "effect": "PSII core destabilization"},
            {"name": "CP43 (PsbC)","uniprot": "P06172", "site": "Thr (N-term)",   "effect": "PSII monomer formation"},
            {"name": "PsbH",       "uniprot": "P06173", "site": "Thr-2",           "effect": "PSII repair — FtsH accessibility"},
            {"name": "STN8",       "uniprot": "Q9LZV4", "site": "Thr-289",         "effect": "Autophosphorylation (activation loop)"},
        ],
        "regulatory_inputs": [
            "Activated by: high-light stress (photoinhibitory conditions)",
            "Activated by: ROS production at PSII (singlet oxygen signal)",
            "Regulated by: PBCP phosphatase (reverses PSII phosphorylation)",
            "Independent of: PQ pool redox state (unlike STN7)",
            "Partially redundant with STN7 for LHCII phosphorylation under extreme conditions",
        ],
        "interaction_partners": [
            {"name": "PBCP",          "role": "Antagonistic phosphatase — reverses PSII core phosphorylation"},
            {"name": "FtsH protease", "role": "Degrades phospho-D1 → enables PSII repair"},
            {"name": "D1 (PsbA)",     "role": "Primary substrate — PSII reaction center"},
            {"name": "D2 (PsbD)",     "role": "Primary substrate — PSII reaction center"},
            {"name": "PsbH",          "role": "PSII core substrate — repair regulation"},
            {"name": "CP43 (PsbC)",   "role": "PSII inner antenna substrate"},
            {"name": "STN7",          "role": "Functional paralog — partial substrate overlap"},
        ],
        "go_terms": {
            "molecular_function": ["GO:0004674 (protein serine/threonine kinase)","GO:0005524 (ATP binding)","GO:0042803 (protein homodimerization)"],
            "biological_process": ["GO:0010098 (photosystem II repair)","GO:0018105 (peptidyl-serine phosphorylation)","GO:0006979 (response to oxidative stress)","GO:0009644 (response to high light intensity)","GO:0030154 (cell differentiation)"],
            "cellular_component": ["GO:0009535 (chloroplast thylakoid membrane)","GO:0009534 (chloroplast thylakoid)","GO:0009654 (photosystem II oxygen evolving complex)"],
        },
        "cellular_pathway": [
            "1. High-light stress → excess photons at PSII → D1 photodamage",
            "2. ROS/singlet oxygen signal activates STN8 kinase",
            "3. STN8 phosphorylates D1 Thr-2, D2 Ser-2, CP43, PsbH",
            "4. Phospho-PSII core destabilizes grana stacking",
            "5. Damaged PSII migrates to stroma lamellae / grana margins",
            "6. FtsH protease (Zn-metalloprotease) degrades phospho-D1",
            "7. New D1 is inserted → PSII reassembly",
            "8. PBCP phosphatase reverses remaining phosphorylation",
            "9. Repaired PSII re-enters grana stacks",
        ],
    },
}

def fetch_uniprot(uniprot_id, timeout=10):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read().decode())
    except Exception as e:
        print(f"    UniProt API unavailable ({e}), using curated data")
        return None

def fetch_string_partners(string_id, tax_id, timeout=10, score=700):
    url = (f"https://string-db.org/api/json/interaction_partners"
           f"?identifiers={string_id}&species={tax_id}&required_score={score}&limit=20")
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read().decode())
    except Exception as e:
        print(f"    STRING API unavailable ({e}), using curated data")
        return None

def extract_uniprot_data(data):
    """Extract GO terms and PTMs from UniProt JSON."""
    if not data:
        return {}, []
    go = {"molecular_function": [], "biological_process": [], "cellular_component": []}
    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "GO":
            go_id   = ref.get("id", "")
            props   = {p["key"]: p["value"] for p in ref.get("properties", [])}
            go_term = props.get("GoTerm", "")
            aspect  = props.get("GoAspect", "")
            if "F" in aspect:   go["molecular_function"].append(f"{go_id} ({go_term})")
            elif "P" in aspect: go["biological_process"].append(f"{go_id} ({go_term})")
            elif "C" in aspect: go["cellular_component"].append(f"{go_id} ({go_term})")

    ptms = []
    for feat in data.get("features", []):
        if feat.get("type") in ("Modified residue", "Natural variant"):
            desc = feat.get("description", "")
            pos  = feat.get("location", {}).get("start", {}).get("value", "?")
            ptms.append({"position": pos, "description": desc})
    return go, ptms

def plot_pathway(name, steps, out_path):
    n = len(steps)
    fig, ax = plt.subplots(figsize=(10, max(6, n * 0.9)), facecolor="#0f172a")
    ax.set_facecolor("#0f172a")
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.5, n + 0.5)
    ax.axis("off")

    colors = ["#1e40af","#1d4ed8","#2563eb","#3b82f6","#60a5fa",
              "#93c5fd","#bfdbfe","#dbeafe","#eff6ff"]
    for i, step in enumerate(reversed(steps)):
        y = i
        col = colors[i % len(colors)]
        rect = mpatches.FancyBboxPatch((0.3, y + 0.05), 9.4, 0.75,
            boxstyle="round,pad=0.08", facecolor=col, edgecolor="#334155",
            linewidth=1.2, zorder=2)
        ax.add_patch(rect)
        ax.text(5, y + 0.43, step, color="white", fontsize=8.5,
                ha="center", va="center", fontweight="bold", wrap=True, zorder=3)
        if i < n - 1:
            ax.annotate("", xy=(5, y + 0.85), xytext=(5, y + 0.95),
                arrowprops=dict(arrowstyle="->", color="#60a5fa", lw=1.5))

    ax.set_title(f"{name} — Intracellular Signaling Pathway", color="white",
                 fontsize=13, fontweight="bold", pad=14)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="#0f172a")
    plt.close()
    print(f"  Saved pathway figure → {os.path.basename(out_path)}")

def plot_interaction_network(name, partners, out_path):
    n = len(partners)
    if n == 0:
        return
    fig, ax = plt.subplots(figsize=(10, 10), facecolor="#0f172a")
    ax.set_facecolor("#0f172a")
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.axis("off")

    ax.add_patch(plt.Circle((0, 0), 0.28, color="#2563eb", zorder=5))
    ax.text(0, 0, name, color="white", fontsize=13, fontweight="bold",
            ha="center", va="center", zorder=6)

    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    role_colors = {
        "Antagonistic phosphatase": "#ef4444",
        "Antagonistic": "#ef4444",
        "Primary substrate": "#22c55e",
        "substrate": "#22c55e",
        "scaffold": "#f59e0b",
        "Required": "#f59e0b",
        "protease": "#a78bfa",
        "Functional paralog": "#94a3b8",
        "source": "#fb923c",
    }

    for i, partner in enumerate(partners):
        angle = angles[i]
        x = 1.2 * np.cos(angle)
        y = 1.2 * np.sin(angle)

        role = partner.get("role", "")
        col  = "#64748b"
        for k, v in role_colors.items():
            if k.lower() in role.lower():
                col = v; break

        ax.plot([0, x * 0.85], [0, y * 0.85], color=col, linewidth=1.4,
                alpha=0.7, zorder=2)
        ax.add_patch(plt.Circle((x, y), 0.18, color=col, zorder=4, alpha=0.9))
        ax.text(x, y, partner["name"], color="white", fontsize=6.5,
                ha="center", va="center", fontweight="bold", zorder=5)
        lx = 1.42 * np.cos(angle)
        ly = 1.42 * np.sin(angle)
        short_role = role[:35] + "…" if len(role) > 35 else role
        ax.text(lx, ly, short_role, color="#94a3b8", fontsize=5.5,
                ha="center", va="center", style="italic")

    ax.set_title(f"{name} Protein Interaction Network", color="white",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="#0f172a")
    plt.close()
    print(f"  Saved network figure → {os.path.basename(out_path)}")

def build_report_text(results):
    lines = []
    lines.append("=" * 70)
    lines.append("  STN7/STN8 INTRACELLULAR PROCESS COMPREHENSIVE REPORT")
    lines.append("  Arabidopsis thaliana Chloroplast Thylakoid Kinases")
    lines.append("=" * 70)

    for name, data in results.items():
        cur = data["curated"]
        lines.append(f"\n{'─'*70}")
        lines.append(f"  {name} — {cur['full_name']}")
        lines.append(f"  Gene: {cur['gene_name']}  |  UniProt: {data['uniprot_id']}")
        lines.append(f"{'─'*70}")

        lines.append("\n  SUBCELLULAR LOCALIZATION:")
        for loc in cur["subcellular_locations"]:
            lines.append(f"    • {loc}")

        lines.append("\n  MOLECULAR FUNCTIONS:")
        for mf in cur["molecular_functions"]:
            lines.append(f"    • {mf}")

        lines.append("\n  BIOLOGICAL PROCESSES:")
        for bp in cur["biological_processes"]:
            lines.append(f"    • {bp}")

        lines.append(f"\n  KNOWN SUBSTRATES ({len(cur['substrates'])}):")
        for s in cur["substrates"]:
            lines.append(f"    • {s['name']} ({s['uniprot']})  site: {s['site']}")
            lines.append(f"      Effect: {s['effect']}")

        lines.append("\n  REGULATORY INPUTS:")
        for ri in cur["regulatory_inputs"]:
            lines.append(f"    • {ri}")

        lines.append(f"\n  INTERACTION PARTNERS ({len(cur['interaction_partners'])}):")
        for ip in cur["interaction_partners"]:
            lines.append(f"    • {ip['name']:<20} → {ip['role']}")

        lines.append("\n  CELLULAR SIGNALING PATHWAY:")
        for step in cur["cellular_pathway"]:
            lines.append(f"    {step}")

        lines.append("\n  GO TERMS:")
        for aspect, terms in cur["go_terms"].items():
            lines.append(f"    {aspect.replace('_',' ').title()}:")
            for t in terms:
                lines.append(f"      - {t}")

    lines.append(f"\n{'='*70}")
    lines.append("  STRUCTURAL-FUNCTION CORRELATION SUMMARY")
    lines.append(f"{'='*70}")
    lines.append("""
  Domain               | STN7 Function              | STN8 Function
  ─────────────────────┼────────────────────────────┼────────────────────────────
  Transit Peptide      | Chloroplast import signal   | Chloroplast import signal
  Transmembrane Anchor | Thylakoid membrane anchor   | Thylakoid membrane anchor
  Kinase Domain        | LHCII phosphorylation       | PSII core phosphorylation
  Activation Loop      | Autophosphorylation (Thr309)| Autophosphorylation (Thr289)
  C-terminal Tail      | Substrate specificity       | Substrate specificity
""")
    return "\n".join(lines)

def main():
    print("\n" + "="*60)
    print("  STN7/STN8 Intracellular Process Analyzer")
    print("="*60)

    results = {}
    for name, meta in PROTEINS.items():
        print(f"\n── {name} ({meta['uniprot']}) ─────────────────────────────")

        # Try live UniProt fetch
        print("  Fetching UniProt data...")
        uniprot_data = fetch_uniprot(meta["uniprot"])
        go_live, ptms_live = extract_uniprot_data(uniprot_data)

        # Try STRING
        print("  Fetching STRING interaction data...")
        string_data = fetch_string_partners(meta["string_id"], meta["tax_id"])

        # Build result merging live + curated
        cur = CURATED[name]
        go_final = go_live if any(go_live.values()) else cur["go_terms"]

        string_partners = []
        if string_data:
            for p in string_data[:10]:
                string_partners.append({
                    "name":  p.get("preferredName_B", p.get("stringId_B", "?")),
                    "score": p.get("score", 0),
                    "role":  "STRING interaction partner",
                })

        results[name] = {
            "uniprot_id":      meta["uniprot"],
            "curated":         cur,
            "go_terms":        go_final,
            "live_ptms":       ptms_live[:30] if ptms_live else [],
            "string_partners": string_partners,
            "n_substrates":    len(cur["substrates"]),
            "n_partners":      len(cur["interaction_partners"]),
            "n_pathway_steps": len(cur["cellular_pathway"]),
        }

        print(f"  GO terms: {sum(len(v) for v in go_final.values())} total")
        print(f"  Substrates: {len(cur['substrates'])}")
        print(f"  Interaction partners: {len(cur['interaction_partners'])}")

    # ── Save JSON ────────────────────────────────────────────────
    json_path = os.path.join(RESULTS_DIR, "cellular_process_report.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  Saved JSON → cellular_process_report.json")

    # ── Save text report ─────────────────────────────────────────
    txt_path = os.path.join(RESULTS_DIR, "cellular_process_report.txt")
    with open(txt_path, "w") as f:
        f.write(build_report_text(results))
    print(f"  Saved TXT  → cellular_process_report.txt")

    # ── Figures ──────────────────────────────────────────────────
    for name in ["STN7", "STN8"]:
        cur = CURATED[name]
        plot_pathway(name, cur["cellular_pathway"],
                     os.path.join(FIGURES_DIR, f"cellular_pathway_{name}.png"))
        plot_interaction_network(name, cur["interaction_partners"],
                                 os.path.join(FIGURES_DIR, f"interaction_network_{name}.png"))

    print("\n✓ Cellular process analysis complete.")

if __name__ == "__main__":
    main()
