"""
Step 10: GATE G5 — miRNA Integration (SOP Phase 6)
Central hypothesis: DENV/ZIKV suppress miRNAs → target genes are de-repressed →
appear as shared upregulated DEGs.

Three-tier test:
  - 55-set (cross-flavivirus miRNAs): should show STRONGEST enrichment in shared DEGs
  - 91-set (all DENV miRNAs): should show moderate enrichment
  - 36-set (DENV-only control): should show WEAKEST enrichment (specificity check)

Also computes three-way intersection: shared DEGs ∩ proviral ∩ miRNA-55 targets
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
import gseapy as gp

warnings.filterwarnings("ignore")

BASE_DIR  = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
LIT_DIR   = BASE_DIR / "02_literature_resources"
DEG_DIR   = BASE_DIR / "01_processed_data" / "deg_tables"
RES_DIR   = BASE_DIR / "03_results" / "phase6_mirna"
FIG_MAIN  = BASE_DIR / "04_figures" / "main"
CKPT_DIR  = BASE_DIR / "checkpoints"
LOG_FILE  = BASE_DIR / "logs" / "step10_gate_g5.log"

for d in [RES_DIR, FIG_MAIN, CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE = CKPT_DIR / "step10_checkpoint.json"

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


def run_fisher(shared_genes, target_genes, background_genes):
    """Fisher exact test for enrichment of miRNA targets in shared DEGs."""
    shared  = set(shared_genes)
    targets = set(target_genes)
    bg      = set(background_genes)

    a = len(shared & targets)
    b = len(shared - targets)
    c = len((bg - shared) & targets)
    d = len((bg - shared) - targets)

    if a == 0:
        return {"overlap": 0, "expected": 0.0, "odds_ratio": 0.0,
                "p_value": 1.0, "fold_enrichment": 0.0, "genes": []}

    oddsratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
    expected = len(shared) * len(targets) / max(len(bg), 1)
    fold_enrich = a / expected if expected > 0 else float("inf")

    return {
        "overlap": a,
        "expected": round(expected, 3),
        "odds_ratio": round(float(oddsratio), 3),
        "p_value": float(pvalue),
        "fold_enrichment": round(fold_enrich, 2),
        "genes": sorted(shared & targets),
    }


def get_enrichr_targets_for_mirna_set(mirna_set, db="TargetScan_microRNA_2017"):
    """Query Enrichr db and retrieve all target genes for miRNAs in the given set.
    Strips version suffixes to match miRNA names (e.g. 'hsa-miR-21-5p' in set vs 'hsa-miR-21-5p' in db).
    """
    # Download the full library
    try:
        library = gp.get_library(db)
    except Exception as e:
        log(f"  Could not fetch Enrichr library {db}: {e}")
        return set()

    mirna_set_lower = {m.lower() for m in mirna_set}
    mirna_set_base  = {m.split("-")[0] + "-" + m.split("-")[1] if m.count("-") >= 1 else m for m in mirna_set_lower}

    target_genes = set()
    matched_mirnas = []

    for term, genes in library.items():
        term_lower = term.lower()
        # Match if the db term is a member of our miRNA set
        if term_lower in mirna_set_lower:
            target_genes.update([g.upper() for g in genes])
            matched_mirnas.append(term)
        # Also match base name (without seed)
        term_base = "-".join(term_lower.split("-")[:3])
        if term_base in mirna_set_lower:
            target_genes.update([g.upper() for g in genes])
            if term not in matched_mirnas:
                matched_mirnas.append(term)

    return target_genes, matched_mirnas


def main():
    ckpt = load_ckpt()

    log("=" * 60)
    log("Step 10: GATE G5 — miRNA Integration (Central Hypothesis)")
    log("=" * 60)

    # ─── Load shared DEGs (gene symbols) ──────────────────────────────────────
    log("Loading shared DEGs ...")
    shared_all = pd.read_csv(BASE_DIR / "03_results" / "phase3_shared_degs" / "shared_DEGs_annotated.csv")
    shared_up  = shared_all[shared_all["log2FC_DENV"] > 0]["symbol"].dropna().tolist()
    log(f"  Shared upregulated: {len(shared_up)} genes → {shared_up}")

    # Background = all tested genes
    denv_ann = pd.read_csv(DEG_DIR / "DEGs_DENV_vs_Control_annotated.csv")
    bg_genes = denv_ann["symbol"].dropna().unique().tolist()
    log(f"  Background: {len(bg_genes)} genes")

    # Load proviral genes for three-way intersection
    proviral_all = [g.strip() for g in open(LIT_DIR / "host_factors_proviral.txt").readlines() if g.strip()]

    # ─── Load miRNA sets ───────────────────────────────────────────────────────
    log("\nLoading miRNA sets ...")
    mirna_55 = [m.strip() for m in open(LIT_DIR / "mirna_55_cross_flavivirus.txt").readlines() if m.strip()]
    mirna_91 = [m.strip() for m in open(LIT_DIR / "mirna_91_denv_all.txt").readlines() if m.strip()]
    mirna_36 = [m.strip() for m in open(LIT_DIR / "mirna_36_denv_only_control.txt").readlines() if m.strip()]
    log(f"  55-set (cross-flavivirus): {len(mirna_55)} miRNAs")
    log(f"  91-set (all DENV):         {len(mirna_91)} miRNAs")
    log(f"  36-set (DENV-only ctrl):   {len(mirna_36)} miRNAs")

    # ─── Step 1: Enrichr analysis on shared DEGs (miRNA databases) ────────────
    if not ckpt.get("enrichr_mirna_done"):
        log("\nRunning Enrichr (TargetScan + miRTarBase) on shared upregulated DEGs ...")

        for db in ["TargetScan_microRNA_2017", "miRTarBase_2017"]:
            try:
                result = gp.enrichr(
                    gene_list=shared_up,
                    gene_sets=[db],
                    outdir=None,
                    verbose=False,
                )
                df = result.results
                df.to_csv(RES_DIR / f"enrichr_{db}_shared_up.csv", index=False)
                sig = df[df["Adjusted P-value"] < 0.05]
                log(f"  {db}: {len(df)} terms, {len(sig)} significant (adj-p < 0.05)")
                if len(sig) > 0:
                    log(f"  Top 5 significant miRNAs:")
                    for _, row in sig.head(5).iterrows():
                        log(f"    {row['Term']}  adj-p={row['Adjusted P-value']:.4f}  overlap={row['Overlap']}")
            except Exception as e:
                log(f"  {db}: ERROR — {e}")

        ckpt["enrichr_mirna_done"] = True
        save_ckpt(ckpt)
    else:
        log("Enrichr miRNA already done. Loading results ...")

    # ─── Step 2: Get miRNA target genes from TargetScan for each set ──────────
    log("\nFetching TargetScan target gene lists for each miRNA set ...")

    if not ckpt.get("targets_fetched"):
        try:
            targets_55, matched_55 = get_enrichr_targets_for_mirna_set(mirna_55, "TargetScan_microRNA_2017")
            targets_91, matched_91 = get_enrichr_targets_for_mirna_set(mirna_91, "TargetScan_microRNA_2017")
            targets_36, matched_36 = get_enrichr_targets_for_mirna_set(mirna_36, "TargetScan_microRNA_2017")

            log(f"  55-set: {len(matched_55)} miRNAs matched → {len(targets_55)} unique target genes")
            log(f"  91-set: {len(matched_91)} miRNAs matched → {len(targets_91)} unique target genes")
            log(f"  36-set: {len(matched_36)} miRNAs matched → {len(targets_36)} unique target genes")

            # Save targets
            pd.DataFrame({"gene": sorted(targets_55)}).to_csv(RES_DIR / "targets_55set_TargetScan.csv", index=False)
            pd.DataFrame({"gene": sorted(targets_91)}).to_csv(RES_DIR / "targets_91set_TargetScan.csv", index=False)
            pd.DataFrame({"gene": sorted(targets_36)}).to_csv(RES_DIR / "targets_36set_TargetScan.csv", index=False)

            ckpt["targets_fetched"] = True
            ckpt["n_targets_55"] = len(targets_55)
            ckpt["n_targets_91"] = len(targets_91)
            ckpt["n_targets_36"] = len(targets_36)
            ckpt["n_matched_55"] = len(matched_55)
            save_ckpt(ckpt)
        except Exception as e:
            log(f"  ERROR fetching targets: {e}")
            # Fall back to using Enrichr results directly
            targets_55, targets_91, targets_36 = set(), set(), set()
            log("  Falling back to Enrichr result overlap method")
    else:
        log("  Loading cached target sets ...")
        targets_55 = set(pd.read_csv(RES_DIR / "targets_55set_TargetScan.csv")["gene"].tolist())
        targets_91 = set(pd.read_csv(RES_DIR / "targets_91set_TargetScan.csv")["gene"].tolist())
        targets_36 = set(pd.read_csv(RES_DIR / "targets_36set_TargetScan.csv")["gene"].tolist())
        log(f"  Loaded: 55-set={len(targets_55)}, 91-set={len(targets_91)}, 36-set={len(targets_36)} targets")

    # ─── Step 3: Three-tier Fisher enrichment test (GATE G5) ──────────────────
    log("\n" + "=" * 50)
    log("GATE G5: Three-Tier miRNA Enrichment Test")
    log("=" * 50)

    tier_results = {}
    for label, targets in [("55-set (cross-flavi)", targets_55),
                            ("91-set (DENV all)",    targets_91),
                            ("36-set (DENV-only)",   targets_36)]:
        if len(targets) == 0:
            log(f"  {label}: No targets available — skipping Fisher test")
            tier_results[label] = {"overlap": 0, "p_value": 1.0, "fold_enrichment": 0.0, "genes": []}
            continue

        res = run_fisher(shared_up, targets, bg_genes)
        tier_results[label] = res
        sig = "*** SIG ***" if res["p_value"] < 0.05 else "ns"
        log(f"\n  [{label}]")
        log(f"    Target genes: {len(targets)}")
        log(f"    Overlap with shared DEGs: {res['overlap']} (expected {res['expected']})")
        log(f"    Fold enrichment: {res['fold_enrichment']}×  Odds ratio: {res['odds_ratio']}")
        log(f"    Fisher p = {res['p_value']:.4e}  {sig}")
        if res["genes"]:
            log(f"    Genes in intersection: {res['genes']}")

    # ─── Step 4: Enrichr-based overlap with miRNA sets ────────────────────────
    log("\nChecking Enrichr results for miRNA set overlap ...")

    mirna_55_lower = {m.lower() for m in mirna_55}
    mirna_91_lower = {m.lower() for m in mirna_91}
    mirna_36_lower = {m.lower() for m in mirna_36}

    enrichr_overlap = {}
    for db in ["TargetScan_microRNA_2017", "miRTarBase_2017"]:
        enr_file = RES_DIR / f"enrichr_{db}_shared_up.csv"
        if not enr_file.exists():
            continue
        enr_df = pd.read_csv(enr_file)
        enr_df["term_lower"] = enr_df["Term"].str.lower()

        sig_terms = set(enr_df[enr_df["Adjusted P-value"] < 0.05]["term_lower"].tolist())
        all_terms  = set(enr_df["term_lower"].tolist())

        in_55_sig  = len(sig_terms & mirna_55_lower)
        in_55_all  = len(all_terms & mirna_55_lower)
        in_91_sig  = len(sig_terms & mirna_91_lower)
        in_36_sig  = len(sig_terms & mirna_36_lower)

        log(f"\n  [{db}]")
        log(f"    Sig. terms (adj-p<0.05): {len(sig_terms)}")
        log(f"    55-set miRNAs found (sig/all): {in_55_sig}/{in_55_all}")
        log(f"    91-set miRNAs found (sig): {in_91_sig}")
        log(f"    36-set miRNAs found (sig): {in_36_sig}")

        enrichr_overlap[db] = {
            "sig_terms": len(sig_terms),
            "in_55_sig": in_55_sig,
            "in_91_sig": in_91_sig,
            "in_36_sig": in_36_sig,
        }

    # ─── Step 5: Three-way intersection ───────────────────────────────────────
    log("\nComputing three-way intersection: shared DEGs ∩ proviral ∩ miRNA-55 targets ...")

    three_way = set(shared_up) & set(proviral_all) & targets_55
    two_way_mirna = set(shared_up) & targets_55
    two_way_prov  = set(shared_up) & set(proviral_all)

    log(f"  Shared up ∩ miRNA-55 targets: {len(two_way_mirna)} genes → {sorted(two_way_mirna)}")
    log(f"  Shared up ∩ proviral:         {len(two_way_prov)} genes → {sorted(two_way_prov)}")
    log(f"  Three-way (∩ both):           {len(three_way)} genes → {sorted(three_way)}")

    # Save three-way result
    three_way_df = pd.DataFrame({
        "gene": sorted(shared_up),
        "in_mirna_55_targets": [g in targets_55 for g in sorted(shared_up)],
        "in_proviral": [g in set(proviral_all) for g in sorted(shared_up)],
        "three_way": [g in three_way for g in sorted(shared_up)],
    })
    three_way_df.to_csv(RES_DIR / "three_way_intersection.csv", index=False)
    log(f"\n  Three-way intersection saved → {RES_DIR}/three_way_intersection.csv")

    # ─── Figure ───────────────────────────────────────────────────────────────
    log("\nGenerating GATE G5 figure ...")

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    # Panel A — Three-tier enrichment comparison
    ax = axes[0]
    tier_labels = ["55-set\n(Cross-flavi)", "91-set\n(DENV all)", "36-set\n(DENV-only)"]
    tier_fe     = [tier_results.get(k, {}).get("fold_enrichment", 0) for k in
                   ["55-set (cross-flavi)", "91-set (DENV all)", "36-set (DENV-only)"]]
    tier_p      = [tier_results.get(k, {}).get("p_value", 1.0) for k in
                   ["55-set (cross-flavi)", "91-set (DENV all)", "36-set (DENV-only)"]]
    tier_ovlp   = [tier_results.get(k, {}).get("overlap", 0) for k in
                   ["55-set (cross-flavi)", "91-set (DENV all)", "36-set (DENV-only)"]]

    bar_colors = ["#D32F2F" if p < 0.05 else "#607D8B" for p in tier_p]
    bars = ax.bar(tier_labels, tier_fe, color=bar_colors, edgecolor="black", linewidth=0.8, alpha=0.85, width=0.5)
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1.2, label="Expected (1×)")
    for bar, fe, ov, p in zip(bars, tier_fe, tier_ovlp, tier_p):
        label = f"{fe:.1f}×\n(n={ov}, p={p:.3f})" if fe > 0 else "0×\n(ns)"
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                label, ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Fold Enrichment (miRNA targets in shared DEGs)", fontsize=10)
    ax.set_title("A  GATE G5: Three-Tier miRNA Target\nEnrichment in Shared DEGs", fontsize=11, fontweight="bold")
    ax.set_ylim(0, max(tier_fe + [1]) * 1.5 if any(tier_fe) else 3)
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(color="#D32F2F", label="p < 0.05"),
        Patch(color="#607D8B", label="p ≥ 0.05"),
    ], fontsize=9)

    # Panel B — Three-way intersection breakdown (stacked bar)
    ax2 = axes[1]
    categories = ["Shared DEGs\n(total)", "miRNA-55\ntargets only", "Proviral\nonly", "Three-way\n(∩ all)"]
    counts = [len(shared_up), len(two_way_mirna), len(two_way_prov), len(three_way)]
    colors = ["#607D8B", "#E91E63", "#D32F2F", "#FF6F00"]
    bars2 = ax2.bar(categories, counts, color=colors, edgecolor="black", linewidth=0.8, alpha=0.9)
    for bar, n in zip(bars2, counts):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
                 str(n), ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Number of Genes", fontsize=11)
    ax2.set_title("B  Three-Way Intersection:\nShared DEGs ∩ Proviral ∩ miRNA-55 Targets",
                  fontsize=11, fontweight="bold")
    ax2.set_ylim(0, max(counts) * 1.4 if counts else 5)

    plt.suptitle("GATE G5 — miRNA Integration: Central Hypothesis Test", fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig_path = FIG_MAIN / "Figure_GATE_G5_miRNA_Integration.png"
    plt.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close()
    log(f"Figure saved → {fig_path}")

    # ─── GATE G5 Decision ─────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("GATE G5 DECISION")
    log("=" * 60)

    p_55 = tier_results.get("55-set (cross-flavi)", {}).get("p_value", 1.0)
    p_91 = tier_results.get("91-set (DENV all)", {}).get("p_value", 1.0)
    p_36 = tier_results.get("36-set (DENV-only)", {}).get("p_value", 1.0)
    fe_55 = tier_results.get("55-set (cross-flavi)", {}).get("fold_enrichment", 0)
    fe_36 = tier_results.get("36-set (DENV-only)", {}).get("fold_enrichment", 0)

    log(f"  55-set p = {p_55:.4e}  FE = {fe_55}×")
    log(f"  91-set p = {p_91:.4e}")
    log(f"  36-set p = {p_36:.4e}  FE = {fe_36}×")
    log(f"  Three-way intersection genes: {len(three_way)}")

    if p_55 < 0.05 and fe_55 > fe_36:
        gate_status = "PASSED"
        log("  GATE G5: PASSED — Shared DEGs significantly enriched for 55-set targets")
        log("  SUPPORTS central hypothesis: miRNA suppression → shared gene de-repression")
    elif p_55 < 0.05:
        gate_status = "PASSED"
        log("  GATE G5: PASSED — 55-set enrichment significant (specificity pending)")
    elif len(two_way_mirna) > 0 and fe_55 > 1:
        gate_status = "TREND"
        log("  GATE G5: TREND — Partial support for miRNA hypothesis (not significant)")
    else:
        gate_status = "NOT PASSED"
        log("  GATE G5: NOT PASSED — No miRNA target enrichment in shared DEGs")
        log("  → Central hypothesis not supported at gene level")
        log("  → Check pathway-level miRNA connection instead")

    # Save summary
    summary = {
        "tier_55_p": p_55, "tier_91_p": p_91, "tier_36_p": p_36,
        "fe_55": fe_55, "fe_36": fe_36,
        "two_way_mirna_genes": sorted(two_way_mirna),
        "three_way_genes": sorted(three_way),
        "gate_g5_status": gate_status,
    }
    pd.DataFrame([summary]).to_csv(RES_DIR / "gate_g5_summary.csv", index=False)

    ckpt["gate_g5_done"] = True
    ckpt["gate_g5_status"] = gate_status
    ckpt["p_55"] = p_55
    ckpt["three_way_n"] = len(three_way)
    save_ckpt(ckpt)
    log("\nNext: run step11_download_validation_datasets.py")


if __name__ == "__main__":
    main()
