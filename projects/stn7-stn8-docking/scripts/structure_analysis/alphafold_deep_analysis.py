#!/usr/bin/env python3
"""
STN7/STN8 AlphaFold Deep Structural Characterization
Domain segmentation, pLDDT stats, active site mapping, phospho-site context.
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import PDB
from Bio.PDB import PPBuilder

BASE_DIR    = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
STRUCTURES  = {
    "STN7": os.path.join(BASE_DIR, "data/structures/Q9S713_AlphaFold.pdb"),
    "STN8": os.path.join(BASE_DIR, "data/structures/Q9LZV4_AlphaFold.pdb"),
}
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

DOMAIN_ANNOTATIONS = {
    "STN7": {
        "transit_peptide":      (1,   59),
        "transmembrane_anchor": (60,  100),
        "kinase_domain":        (101, 430),
        "activation_loop":      (281, 310),
        "c_terminal_tail":      (431, 565),
    },
    "STN8": {
        "transit_peptide":      (1,   59),
        "transmembrane_anchor": (60,  95),
        "kinase_domain":        (96,  390),
        "activation_loop":      (261, 290),
        "c_terminal_tail":      (391, 517),
    },
}

CATALYTIC_MOTIFS = {
    "STN7": {
        "P_loop_GxGxxG":  (130, 136),
        "HxD_catalytic":  (228, 230),
        "DFG_motif":      (284, 286),
        "APE_motif":      (309, 311),
    },
    "STN8": {
        "P_loop_GxGxxG":  (112, 118),
        "HxD_catalytic":  (208, 210),
        "DFG_motif":      (264, 266),
        "APE_motif":      (289, 291),
    },
}

D3_TO_1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
    'GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
    'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'
}
PHOSPHO_AA = {'SER', 'THR', 'TYR'}
DOMAIN_COLORS = {
    "transit_peptide":      "#94a3b8",
    "transmembrane_anchor": "#f97316",
    "kinase_domain":        "#3b82f6",
    "activation_loop":      "#ef4444",
    "c_terminal_tail":      "#8b5cf6",
}

def plddt_tier(s):
    if s >= 90: return "Very High (>90)"
    if s >= 70: return "High (70-90)"
    if s >= 50: return "Low (50-70)"
    return "Very Low (<50)"

def load_residues(name, pdb_path):
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure(name, pdb_path)
    chain  = list(struct[0].get_chains())[0]
    residues = []
    for res in chain:
        if not PDB.Polypeptide.is_aa(res):
            continue
        bfs = [a.get_bfactor() for a in res.get_atoms()]
        residues.append({
            "res_num":    res.get_id()[1],
            "res_name":   res.get_resname().upper(),
            "one_letter": D3_TO_1.get(res.get_resname().upper(), 'X'),
            "plddt":      float(np.mean(bfs)) if bfs else 0.0,
            "bf_std":     float(np.std(bfs))  if len(bfs) > 1 else 0.0,
        })
    return residues

def assign_domain(res_num, domains):
    for d in ["activation_loop","transit_peptide","transmembrane_anchor","kinase_domain","c_terminal_tail"]:
        if d in domains:
            s, e = domains[d]
            if s <= res_num <= e:
                return d
    return "unassigned"

def domain_stats(residues, domains):
    buckets = {d: [] for d in list(domains.keys()) + ["unassigned"]}
    for r in residues:
        buckets[assign_domain(r["res_num"], domains)].append(r["plddt"])
    out = {}
    for d, scores in buckets.items():
        if not scores: continue
        a = np.array(scores)
        out[d] = {
            "n_residues": len(scores),
            "mean_plddt": round(float(a.mean()), 2),
            "std_plddt":  round(float(a.std()),  2),
            "confidence_tier": plddt_tier(a.mean()),
        }
    return out

def phospho_sites(residues, domains):
    return [{
        "res_num": r["res_num"], "res_name": r["res_name"],
        "one_letter": r["one_letter"], "plddt": r["plddt"],
        "domain": assign_domain(r["res_num"], domains),
        "confidence_tier": plddt_tier(r["plddt"]),
        "surface_exposure": "exposed" if r["bf_std"] > 8 else "buried",
    } for r in residues if r["res_name"] in PHOSPHO_AA]

def motif_report(residues, motifs, domains):
    res_map = {r["res_num"]: r for r in residues}
    report  = {}
    for name, (s, e) in motifs.items():
        window = [res_map[n] for n in range(s, e+1) if n in res_map]
        if not window:
            report[name] = {"residues": f"{s}-{e}", "sequence": "", "mean_plddt": None}
            continue
        seq  = "".join(r["one_letter"] for r in window)
        mean = round(float(np.mean([r["plddt"] for r in window])), 2)
        report[name] = {
            "residues": f"{s}-{e}", "sequence": seq,
            "mean_plddt": mean, "domain": assign_domain(s, domains),
            "confidence_tier": plddt_tier(mean),
        }
    return report

def ss_counts(pdb_path, domains):
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("tmp", pdb_path)
    ppb    = PPBuilder()
    counts = {d: {"helix": 0, "strand": 0} for d in list(domains.keys()) + ["unassigned"]}
    for pp in ppb.build_peptides(struct[0]):
        for i, (phi, psi) in enumerate(pp.get_phi_psi_list()):
            if phi is None or psi is None: continue
            rn  = pp[i].get_id()[1]
            dom = assign_domain(rn, domains)
            pd, ps = np.degrees(phi), np.degrees(psi)
            if -80 <= pd <= -40 and -60 <= ps <= -20:
                counts[dom]["helix"]  += 1
            elif -160 <= pd <= -60 and 90 <= ps <= 170:
                counts[dom]["strand"] += 1
    return counts

def plot_domain_map(name, residues, domains, motifs, out_path):
    res_nums   = [r["res_num"] for r in residues]
    plddt_vals = [r["plddt"]   for r in residues]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 7),
        gridspec_kw={"height_ratios": [1, 0.18]}, facecolor="#0f172a")
    fig.subplots_adjust(hspace=0.06)

    ax1.set_facecolor("#1e293b")
    ax1.fill_between(res_nums, plddt_vals, alpha=0.22, color="#3b82f6")
    ax1.plot(res_nums, plddt_vals, color="#60a5fa", linewidth=1.2, zorder=3)

    for y, label, col in [(90,"Very High","#22c55e"),(70,"High","#facc15"),(50,"Low","#f97316")]:
        ax1.axhline(y, color=col, linewidth=0.8, linestyle="--", alpha=0.55)
        ax1.text(res_nums[-1]+2, y+0.5, label, color=col, fontsize=7)

    for mname, (ms, me) in motifs.items():
        mid = (ms + me) / 2
        ax1.axvspan(ms, me, color="#f59e0b", alpha=0.30, zorder=2)
        ax1.text(mid, 96, mname.replace("_"," ").replace("motif","").strip(),
                 color="#fbbf24", fontsize=6, ha="center", fontweight="bold", rotation=40)

    ax1.set_xlim(res_nums[0], res_nums[-1])
    ax1.set_ylim(0, 106)
    ax1.set_ylabel("pLDDT Confidence", color="white", fontsize=10)
    ax1.set_xlabel("Residue Position", color="white", fontsize=10)
    ax1.tick_params(colors="white", labelsize=8)
    for sp in ax1.spines.values(): sp.set_edgecolor("#334155")
    ax1.set_title(f"{name} — AlphaFold Deep Structural Analysis",
                  color="white", fontsize=13, fontweight="bold", pad=10)

    ax2.set_facecolor("#0f172a")
    for r in residues:
        dom = assign_domain(r["res_num"], domains)
        ax2.barh(0, 1, left=r["res_num"], height=0.6,
                 color=DOMAIN_COLORS.get(dom, "#e2e8f0"), align="center")
    ax2.set_xlim(res_nums[0], res_nums[-1])
    ax2.axis("off")

    handles = [mpatches.Patch(color=col, label=dom.replace("_"," ").title())
               for dom, col in DOMAIN_COLORS.items()]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=8, labelcolor="white", bbox_to_anchor=(0.5, -0.02))

    plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="#0f172a")
    plt.close()
    print(f"  Saved → {os.path.basename(out_path)}")

def analyze(name):
    pdb_path = STRUCTURES[name]
    print(f"\n{'='*58}\n  Analyzing {name}: {os.path.basename(pdb_path)}\n{'='*58}")
    if not os.path.exists(pdb_path):
        print(f"  ERROR: {pdb_path} not found"); return None

    domains  = DOMAIN_ANNOTATIONS[name]
    motifs   = CATALYTIC_MOTIFS[name]
    residues = load_residues(name, pdb_path)
    print(f"  Loaded {len(residues)} residues")

    ds   = domain_stats(residues, domains)
    ps   = phospho_sites(residues, domains)
    mr   = motif_report(residues, motifs, domains)
    ssc  = ss_counts(pdb_path, domains)

    all_p = [r["plddt"] for r in residues]
    gs = {
        "total_residues":    len(residues),
        "global_mean_plddt": round(float(np.mean(all_p)), 2),
        "global_std_plddt":  round(float(np.std(all_p)),  2),
        "very_high_pct": round(sum(1 for p in all_p if p>=90)/len(all_p)*100, 1),
        "high_pct":      round(sum(1 for p in all_p if 70<=p<90)/len(all_p)*100, 1),
        "low_pct":       round(sum(1 for p in all_p if 50<=p<70)/len(all_p)*100, 1),
        "very_low_pct":  round(sum(1 for p in all_p if p<50)/len(all_p)*100, 1),
    }

    result = {
        "protein": name,
        "uniprot_id": "Q9S713" if name=="STN7" else "Q9LZV4",
        "pdb_source": os.path.basename(pdb_path),
        "global_stats": gs,
        "domain_boundaries": {d: list(v) for d, v in domains.items()},
        "domain_plddt_stats": ds,
        "secondary_structure_per_domain": ssc,
        "catalytic_motifs": mr,
        "phospho_sites_count": len(ps),
        "top_phospho_sites": sorted(ps, key=lambda x: -x["plddt"])[:20],
        "all_phospho_sites": ps,
    }

    json_path = os.path.join(RESULTS_DIR, f"alphafold_deep_analysis_{name}.json")
    with open(json_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Saved JSON → {os.path.basename(json_path)}")

    fig_path = os.path.join(FIGURES_DIR, f"alphafold_domain_map_{name}.png")
    plot_domain_map(name, residues, domains, motifs, fig_path)

    # Print summary
    print(f"\n  Global pLDDT: {gs['global_mean_plddt']:.1f} ± {gs['global_std_plddt']:.1f}")
    print(f"  Confidence dist: VH={gs['very_high_pct']}%  H={gs['high_pct']}%  L={gs['low_pct']}%  VL={gs['very_low_pct']}%")
    print(f"  Phosphorylatable residues: {len(ps)}")
    for dom, stats in ds.items():
        print(f"    {dom:<26}: {stats['mean_plddt']:.1f} ± {stats['std_plddt']:.1f}  [{stats['confidence_tier']}]")
    for mname, data in mr.items():
        if data["mean_plddt"]:
            print(f"    {mname:<24}: pLDDT={data['mean_plddt']:.1f}  seq={data['sequence']}")

    return result

def main():
    print("\n" + "="*58)
    print("  STN7/STN8 AlphaFold Deep Structural Characterization")
    print("="*58)
    for name in ["STN7", "STN8"]:
        analyze(name)
    print("\n✓ Deep analysis complete.")

if __name__ == "__main__":
    main()
