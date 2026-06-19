"""
Step 09: GATE G4 — Formal Proviral Gene Enrichment (SOP Phase 5)
Fisher exact test: Are shared DEGs significantly enriched for proviral host factors?
Also tests antiviral enrichment, cross-flavivirus proviral, and directional analysis.
Checkpoint-based: safe to restart.
"""

import json, time, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import fisher_exact

warnings.filterwarnings("ignore")

BASE_DIR  = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
DEG_DIR   = BASE_DIR / "01_processed_data" / "deg_tables"
LIT_DIR   = BASE_DIR / "02_literature_resources"
RES_DIR   = BASE_DIR / "03_results" / "phase5_host_factors"
FIG_MAIN  = BASE_DIR / "04_figures" / "main"
CKPT_DIR  = BASE_DIR / "checkpoints"
LOG_FILE  = BASE_DIR / "logs" / "step09_gate_g4.log"

for d in [RES_DIR, FIG_MAIN, CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE = CKPT_DIR / "step09_checkpoint.json"

def load_ckpt():
    return json.load(open(CKPT_FILE)) if CKPT_FILE.exists() else {}

def save_ckpt(d):
    json.dump(d, open(CKPT_FILE, "w"), indent=2)

def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")


def run_fisher(shared_genes, gene_set, background_genes):
    """2×2 Fisher exact test for enrichment of gene_set in shared_genes vs background."""
    shared  = set(shared_genes)
    gset    = set(gene_set)
    bg      = set(background_genes)

    a = len(shared & gset)               # shared DEGs in set
    b = len(shared - gset)               # shared DEGs NOT in set
    c = len((bg - shared) & gset)        # background (non-shared) in set
    d = len((bg - shared) - gset)        # background NOT in set

    if a == 0:
        return {"a": a, "b": b, "c": c, "d": d, "odds_ratio": 0.0, "p_value": 1.0,
                "fold_enrichment": 0.0, "in_shared": [], "expected": 0.0}

    oddsratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
    expected = len(shared) * len(gset) / max(len(bg), 1)
    fold_enrich = a / expected if expected > 0 else float("inf")

    return {
        "a": a, "b": b, "c": c, "d": d,
        "odds_ratio": round(float(oddsratio), 3),
        "p_value": float(pvalue),
        "fold_enrichment": round(fold_enrich, 2),
        "expected": round(expected, 2),
        "in_shared": sorted(shared & gset),
    }


def main():
    ckpt = load_ckpt()

    log("=" * 60)
    log("Step 09: GATE G4 — Proviral Enrichment Fisher Test")
    log("=" * 60)

    # ─── Load shared DEGs (gene symbols) ──────────────────────────────────────
    log("Loading shared DEGs ...")
    shared_all = pd.read_csv(BASE_DIR / "03_results" / "phase3_shared_degs" / "shared_DEGs_annotated.csv")
    shared_up  = shared_all[shared_all["log2FC_DENV"] > 0]
    shared_dn  = shared_all[shared_all["log2FC_DENV"] < 0]
    log(f"  Shared all: {len(shared_all)} genes")
    log(f"  Shared up:  {len(shared_up)} genes → {shared_up['symbol'].tolist()}")
    log(f"  Shared down:{len(shared_dn)} genes")

    # ─── Load background (all genes tested in DESeq2) ─────────────────────────
    log("Loading background gene set ...")
    denv_all = pd.read_csv(DEG_DIR / "DEGs_DENV_vs_Control_annotated_full.csv") \
               if (DEG_DIR / "DEGs_DENV_vs_Control_annotated_full.csv").exists() \
               else pd.read_csv(DEG_DIR / "DEGs_DENV_vs_Control_annotated.csv")
    bg_genes = denv_all["symbol"].dropna().unique().tolist()
    log(f"  Background: {len(bg_genes)} genes")

    # ─── Load host factor sets ────────────────────────────────────────────────
    log("Loading host factor sets ...")
    proviral_all   = [g.strip() for g in open(LIT_DIR / "host_factors_proviral.txt").readlines() if g.strip()]
    antiviral_all  = [g.strip() for g in open(LIT_DIR / "host_factors_antiviral.txt").readlines() if g.strip()]
    cross_proviral = [g.strip() for g in open(LIT_DIR / "host_factors_cross_proviral.txt").readlines() if g.strip()]
    cross_antiviral= [g.strip() for g in open(LIT_DIR / "host_factors_cross_antiviral.txt").readlines() if g.strip()]
    curated        = pd.read_csv(LIT_DIR / "host_factors_curated.csv")
    log(f"  Proviral (all): {len(proviral_all)} genes")
    log(f"  Antiviral (all): {len(antiviral_all)} genes")
    log(f"  Cross-flavivirus proviral: {len(cross_proviral)} genes")
    log(f"  Cross-flavivirus antiviral: {len(cross_antiviral)} genes")

    # ─── Fisher tests ──────────────────────────────────────────────────────────
    log("\nRunning Fisher exact tests (GATE G4) ...")

    tests = {
        "shared_up_vs_proviral_all":      (shared_up["symbol"],  proviral_all),
        "shared_up_vs_antiviral_all":     (shared_up["symbol"],  antiviral_all),
        "shared_up_vs_cross_proviral":    (shared_up["symbol"],  cross_proviral),
        "shared_up_vs_cross_antiviral":   (shared_up["symbol"],  cross_antiviral),
        "shared_all_vs_proviral_all":     (shared_all["symbol"], proviral_all),
        "shared_all_vs_antiviral_all":    (shared_all["symbol"], antiviral_all),
    }

    results = {}
    for test_name, (genes, gset) in tests.items():
        res = run_fisher(genes.tolist(), gset, bg_genes)
        results[test_name] = res
        sig = "*** SIGNIFICANT ***" if res["p_value"] < 0.05 else "ns"
        log(f"\n  [{test_name}]")
        log(f"    Overlap (a): {res['a']} | Not in set (b): {res['b']}")
        log(f"    Background in set (c): {res['c']} | Neither (d): {res['d']}")
        log(f"    Expected overlap: {res['expected']}")
        log(f"    Odds ratio: {res['odds_ratio']}  Fold enrichment: {res['fold_enrichment']}x")
        log(f"    Fisher p = {res['p_value']:.4e}  {sig}")
        if res["in_shared"]:
            log(f"    Genes in intersection: {res['in_shared']}")

    # ─── Annotate shared DEGs with proviral/antiviral labels ──────────────────
    log("\nAnnotating shared DEGs with proviral/antiviral classifications ...")
    shared_ann = shared_all.copy()
    shared_ann["in_proviral_all"]      = shared_ann["symbol"].isin(proviral_all)
    shared_ann["in_antiviral_all"]     = shared_ann["symbol"].isin(antiviral_all)
    shared_ann["in_cross_proviral"]    = shared_ann["symbol"].isin(cross_proviral)
    shared_ann["in_cross_antiviral"]   = shared_ann["symbol"].isin(cross_antiviral)
    shared_ann["host_factor_class"] = "Novel"
    shared_ann.loc[shared_ann["in_antiviral_all"],  "host_factor_class"] = "Antiviral"
    shared_ann.loc[shared_ann["in_proviral_all"],   "host_factor_class"] = "Proviral"

    shared_ann.to_csv(RES_DIR / "shared_DEGs_hostfactor_annotated.csv", index=False)
    log(f"  Saved: shared_DEGs_hostfactor_annotated.csv")
    log(f"  Novel: {(shared_ann['host_factor_class']=='Novel').sum()}")
    log(f"  Proviral: {(shared_ann['host_factor_class']=='Proviral').sum()}")
    log(f"  Antiviral: {(shared_ann['host_factor_class']=='Antiviral').sum()}")

    # ─── Build GATE G4 summary table ──────────────────────────────────────────
    gate_rows = []
    for test_name, res in results.items():
        gate_rows.append({
            "test": test_name,
            "n_shared": res["a"] + res["b"],
            "overlap": res["a"],
            "expected": res["expected"],
            "fold_enrichment": res["fold_enrichment"],
            "odds_ratio": res["odds_ratio"],
            "fisher_p": res["p_value"],
            "significant": res["p_value"] < 0.05,
            "genes_in_overlap": "; ".join(res["in_shared"]),
        })
    gate_df = pd.DataFrame(gate_rows)
    gate_df.to_csv(RES_DIR / "gate_g4_results.csv", index=False)
    log(f"\nGATE G4 results saved → {RES_DIR}/gate_g4_results.csv")

    # ─── Figure ───────────────────────────────────────────────────────────────
    log("\nGenerating GATE G4 figure ...")

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    # Panel A — Fisher test summary (forest-style)
    ax = axes[0]
    test_labels = [
        "Shared up vs\nProviral (all)",
        "Shared up vs\nAntiviral (all)",
        "Shared up vs\nCross-PROV",
        "Shared up vs\nCross-ANTI",
    ]
    test_keys = ["shared_up_vs_proviral_all", "shared_up_vs_antiviral_all",
                 "shared_up_vs_cross_proviral", "shared_up_vs_cross_antiviral"]
    fe_vals = [results[k]["fold_enrichment"] for k in test_keys]
    p_vals  = [results[k]["p_value"] for k in test_keys]
    colors  = ["#D32F2F" if p < 0.05 else "#9E9E9E" for p in p_vals]

    y = range(len(test_labels))
    bars = ax.barh(list(y), fe_vals, color=colors, edgecolor="black", alpha=0.85)
    ax.axvline(1.0, color="black", linestyle="--", linewidth=1.2, label="Expected (1×)")
    for i, (fe, p) in enumerate(zip(fe_vals, p_vals)):
        label = f"{fe:.1f}× (p={p:.3f})" if fe > 0 else "0× (p=1.0)"
        ax.text(fe + 0.05, i, label, va="center", fontsize=9)
    ax.set_yticks(list(y))
    ax.set_yticklabels(test_labels, fontsize=10)
    ax.set_xlabel("Fold Enrichment", fontsize=11)
    ax.set_title("A  GATE G4: Proviral/Antiviral Enrichment\nin Shared DEGs (Fisher Exact Test)",
                 fontsize=11, fontweight="bold")
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color="#D32F2F", label="p < 0.05 (significant)"),
        Patch(color="#9E9E9E", label="p ≥ 0.05 (not significant)"),
    ], fontsize=9, loc="lower right")

    # Panel B — Shared DEG annotation pie/bar chart
    ax2 = axes[1]
    class_counts = shared_ann["host_factor_class"].value_counts()
    colors_pie = {"Novel": "#607D8B", "Proviral": "#D32F2F", "Antiviral": "#1565C0"}
    pie_colors = [colors_pie.get(c, "#999") for c in class_counts.index]
    wedges, texts, autotexts = ax2.pie(
        class_counts.values, labels=class_counts.index,
        colors=pie_colors, autopct="%1.0f%%", startangle=140,
        textprops={"fontsize": 11},
        wedgeprops={"edgecolor": "white", "linewidth": 2}
    )
    for at in autotexts:
        at.set_fontsize(10)
        at.set_fontweight("bold")
    ax2.set_title(f"B  Shared DEGs (n={len(shared_all)}):\nHost Factor Classification",
                  fontsize=11, fontweight="bold")

    plt.suptitle("GATE G4 — Host Factor Enrichment Analysis", fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig_path = FIG_MAIN / "Figure_GATE_G4_Proviral_Enrichment.png"
    plt.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close()
    log(f"Figure saved → {fig_path}")

    # ─── GATE G4 Decision ─────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("GATE G4 DECISION")
    log("=" * 60)

    g4_prov = results["shared_up_vs_proviral_all"]["p_value"]
    g4_anti = results["shared_up_vs_antiviral_all"]["p_value"]
    g4_cross = results["shared_up_vs_cross_proviral"]["p_value"]

    log(f"  Proviral enrichment p = {g4_prov:.4e}  (threshold: p < 0.05)")
    log(f"  Antiviral enrichment p = {g4_anti:.4e}  (threshold: p < 0.05)")
    log(f"  Cross-proviral enrichment p = {g4_cross:.4e}  (threshold: p < 0.05)")

    any_sig = any(p < 0.05 for p in [g4_prov, g4_anti, g4_cross])
    if g4_prov < 0.05:
        log("  GATE G4 STATUS: PASSED — Shared DEGs enriched for proviral genes")
        log("  → Can proceed to Phase 6 (miRNA integration)")
        gate_status = "PASSED"
    elif g4_anti < 0.05:
        log("  GATE G4 STATUS: PARTIAL PASS — Enriched for antiviral (not proviral)")
        log("  → Convergence reflects immune defense, not proviral exploitation")
        gate_status = "PARTIAL"
    else:
        log("  GATE G4 STATUS: NOT PASSED — No significant host factor enrichment")
        log("  → The 15 shared genes are novel candidates outside known factor lists")
        log("  → Per SOP: proceed but note G4 status; miRNA integration still warranted")
        gate_status = "NOT PASSED"

    ckpt["gate_g4_done"] = True
    ckpt["gate_g4_status"] = gate_status
    ckpt["gate_g4_proviral_p"] = g4_prov
    ckpt["gate_g4_antiviral_p"] = g4_anti
    save_ckpt(ckpt)
    log("\nNext: run step10_gate_g5_mirna_integration.py")


if __name__ == "__main__":
    main()
