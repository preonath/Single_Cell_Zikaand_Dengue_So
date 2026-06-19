"""
Step 03b: Fair comparison — DENV MOI=1 only vs ZIKV MOI=1
Rationale: original DENV used MOI=1+10 combined; ZIKV only MOI=1.
This creates an unfair comparison. This script redoes DENV DEG with MOI=1 only,
then compares results side-by-side with the original.
Checkpoint-based.
"""

import json, time, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import sparse
from scipy.stats import pearsonr, fisher_exact
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")

BASE_DIR = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
ANN_DIR  = BASE_DIR / "01_processed_data" / "anndata_objects"
PB_DIR   = BASE_DIR / "01_processed_data" / "pseudobulk"
DEG_DIR  = BASE_DIR / "01_processed_data" / "deg_tables"
RES_DIR  = BASE_DIR / "03_results" / "phase3_shared_degs"
FIG_MAIN = BASE_DIR / "04_figures" / "main"
FIG_SUPP = BASE_DIR / "04_figures" / "supplementary"
CKPT_DIR = BASE_DIR / "checkpoints"
LOG_FILE = BASE_DIR / "logs" / "step03b_fair.log"

for d in [PB_DIR, DEG_DIR, RES_DIR, FIG_MAIN, FIG_SUPP, CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE = CKPT_DIR / "step03b_checkpoint.json"
FC_THRESHOLD = 1.0
P_THRESHOLD  = 0.05

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


# ─── Build pseudobulk for MOI=1 only ─────────────────────────────────────────
def make_pseudobulk_moi1(adata_raw: sc.AnnData,
                          adata_filt: sc.AnnData) -> tuple:
    log("Building pseudobulk — MOI=1 only (fair comparison) ...")

    keep = adata_filt.obs_names
    adata_raw = adata_raw[adata_raw.obs_names.isin(keep)].copy()

    # MOI=1 or MOI=0 (control) only
    mask_moi = (
        (adata_raw.obs["moi"].astype(int) == 1) |
        (adata_raw.obs["moi"].astype(int) == 0)
    )
    adata_moi1 = adata_raw[mask_moi].copy()

    log(f"  Cells after MOI filter: {adata_moi1.n_obs}")
    log(f"  Condition breakdown: {adata_moi1.obs['condition'].value_counts().to_dict()}")

    pb_counts  = {}
    sample_rows = []
    for cond in ["Control", "DENV", "ZIKV"]:
        for tp in ["4h", "12h", "24h", "48h"]:
            mask = (
                (adata_moi1.obs["condition"] == cond) &
                (adata_moi1.obs["timepoint"] == tp)
            )
            n = mask.sum()
            if n == 0:
                continue
            X_sub = adata_moi1[mask].X
            agg = np.array(X_sub.sum(axis=0)).flatten() if sparse.issparse(X_sub) \
                  else X_sub.sum(axis=0).flatten()
            name = f"{cond}_{tp}"
            pb_counts[name] = agg.astype(int)
            sample_rows.append({
                "sample": name, "condition": cond,
                "timepoint": tp, "n_cells": int(n)
            })
            log(f"  {name}: {n} cells, total counts = {int(agg.sum()):,}")

    count_df = pd.DataFrame(pb_counts, index=adata_moi1.var_names)
    info_df  = pd.DataFrame(sample_rows).set_index("sample")

    count_df.to_csv(PB_DIR / "pseudobulk_counts_moi1.csv")
    info_df.to_csv(PB_DIR / "sample_info_moi1.csv")
    log(f"\n{info_df.to_string()}\n")
    return count_df, info_df


# ─── PyDESeq2 ─────────────────────────────────────────────────────────────────
def run_deseq2(count_df, info_df, contrast_condition, label):
    log(f"\nRunning DESeq2: {label} ...")
    keep = info_df[info_df["condition"].isin(["Control", contrast_condition])].index
    counts_sub = count_df[keep].T
    info_sub   = info_df.loc[keep].copy()

    counts_sub = counts_sub.loc[:, (counts_sub > 0).any(axis=0)]
    counts_sub = counts_sub.loc[:, counts_sub.sum(axis=0) >= 10]
    log(f"  Samples: {len(keep)}, Genes: {counts_sub.shape[1]}")

    dds = DeseqDataSet(
        counts=counts_sub, metadata=info_sub[["condition","timepoint"]],
        design_factors=["timepoint","condition"],
        refit_cooks=True, quiet=True
    )
    dds.deseq2()
    stat_res = DeseqStats(
        dds, contrast=["condition", contrast_condition, "Control"], quiet=True
    )
    stat_res.summary()
    res = stat_res.results_df.copy().reset_index().rename(columns={"index":"gene_id"})
    n_sig = ((res["padj"] < P_THRESHOLD) &
             (res["log2FoldChange"].abs() >= FC_THRESHOLD)).sum()
    log(f"  Significant DEGs: {n_sig}")
    return res


# ─── Classify ─────────────────────────────────────────────────────────────────
def classify(res, label):
    sig = res[(res["padj"] < P_THRESHOLD) &
              (res["log2FoldChange"].abs() >= FC_THRESHOLD) &
              res["padj"].notna()]
    up   = sig[sig["log2FoldChange"] > 0]["gene_id"].tolist()
    down = sig[sig["log2FoldChange"] < 0]["gene_id"].tolist()
    log(f"  {label} — Up: {len(up)}, Down: {len(down)}, Total: {len(up)+len(down)}")
    return {"up": up, "down": down, "all": up + down}


# ─── Shared DEGs ──────────────────────────────────────────────────────────────
def find_shared(denv_degs, zikv_degs, all_genes, label):
    shared_up   = list(set(denv_degs["up"])   & set(zikv_degs["up"]))
    shared_down = list(set(denv_degs["down"]) & set(zikv_degs["down"]))
    shared_all  = shared_up + shared_down
    discordant  = list(
        (set(denv_degs["up"]) & set(zikv_degs["down"])) |
        (set(denv_degs["down"]) & set(zikv_degs["up"]))
    )

    N  = len(set(all_genes))
    nd = len(denv_degs["up"]); nz = len(zikv_degs["up"]); ns = len(shared_up)
    expected = nd * nz / N if N > 0 else 0
    fold_enr = ns / expected if expected > 0 else float("inf")
    contingency = [[ns, nd-ns],[nz-ns, N-nd-nz+ns]]
    _, pval = fisher_exact(contingency, alternative="greater")

    log(f"\n  [{label}] Shared up: {len(shared_up)}, Shared down: {len(shared_down)}")
    log(f"  [{label}] Discordant: {len(discordant)}")
    log(f"  [{label}] Fisher p: {pval:.2e}, Fold enrichment: {fold_enr:.2f}x")
    return shared_up, shared_down, shared_all, pval, fold_enr


# ─── Comparison plot: Original vs Fair ────────────────────────────────────────
def make_comparison_figure(orig: dict, fair: dict):
    log("\nGenerating comparison figure: Original vs Fair (MOI=1) ...")

    labels   = ["DENV DEGs", "ZIKV DEGs", "Shared\n(up)", "Shared\n(down)",
                 "Pearson r", "Fold\nenrichment"]
    orig_vals = [orig["n_denv"], orig["n_zikv"], orig["shared_up"],
                 orig["shared_down"], orig["pearson_r"], min(orig["fold_enr"], 30)]
    fair_vals = [fair["n_denv"], fair["n_zikv"], fair["shared_up"],
                 fair["shared_down"], fair["pearson_r"], min(fair["fold_enr"], 30)]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Bar chart comparison
    ax = axes[0]
    x = np.arange(4)
    w = 0.35
    ax.bar(x - w/2, [orig_vals[i] for i in range(4)], w,
           label="Original (DENV MOI=1+10)", color="#E41A1C", alpha=0.8)
    ax.bar(x + w/2, [fair_vals[i] for i in range(4)], w,
           label="Fair (DENV MOI=1 only)", color="#377EB8", alpha=0.8)
    for bars in ax.containers:
        ax.bar_label(bars, padding=2, fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(["DENV\nDEGs","ZIKV\nDEGs","Shared\nUp","Shared\nDown"],
                       fontsize=10)
    ax.set_ylabel("Number of genes"); ax.legend(fontsize=9)
    ax.set_title("DEG counts: Original vs Fair comparison", fontweight="bold")

    # Pearson r + Fold enrichment
    ax = axes[1]
    metrics = ["Pearson r", "Fold enrichment (up)"]
    o_vals  = [orig["pearson_r"], min(orig["fold_enr"], 30)]
    f_vals  = [fair["pearson_r"], min(fair["fold_enr"], 30)]
    x2 = np.arange(2)
    ax.bar(x2 - w/2, o_vals, w, label="Original", color="#E41A1C", alpha=0.8)
    ax.bar(x2 + w/2, f_vals, w, label="Fair MOI=1", color="#377EB8", alpha=0.8)
    for bars in ax.containers:
        ax.bar_label(bars, fmt="%.2f", padding=2, fontsize=9)
    ax.set_xticks(x2); ax.set_xticklabels(metrics, fontsize=10)
    ax.set_title("Convergence metrics", fontweight="bold")
    ax.legend(fontsize=9)

    fig.suptitle("Effect of MOI correction on DENV-ZIKV convergence analysis",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    fig.savefig(FIG_SUPP / "MOI_correction_comparison.pdf", bbox_inches="tight")
    fig.savefig(FIG_SUPP / "MOI_correction_comparison.png", bbox_inches="tight", dpi=150)
    plt.close()
    log("  Comparison figure saved")


# ─── FC Correlation ───────────────────────────────────────────────────────────
def fc_correlation(denv_res, zikv_res, shared_up, label):
    merged = pd.merge(
        denv_res[["gene_id","log2FoldChange"]].rename(columns={"log2FoldChange":"FC_DENV"}),
        zikv_res[["gene_id","log2FoldChange"]].rename(columns={"log2FoldChange":"FC_ZIKV"}),
        on="gene_id"
    ).dropna()
    r, p = pearsonr(merged["FC_DENV"], merged["FC_ZIKV"])
    log(f"  [{label}] Pearson r = {r:.4f}, p = {p:.2e}")

    fig, ax = plt.subplots(figsize=(7, 6))
    is_shared = merged["gene_id"].isin(shared_up)
    ax.scatter(merged.loc[~is_shared,"FC_DENV"], merged.loc[~is_shared,"FC_ZIKV"],
               c="lightgrey", s=4, alpha=0.3, rasterized=True, label="Other")
    ax.scatter(merged.loc[is_shared,"FC_DENV"], merged.loc[is_shared,"FC_ZIKV"],
               c="#D73027", s=40, alpha=0.9, zorder=5,
               edgecolors="black", linewidths=0.5,
               label=f"Shared up (n={is_shared.sum()})")
    lim = max(abs(merged["FC_DENV"].quantile([0.01,0.99])).max(),
              abs(merged["FC_ZIKV"].quantile([0.01,0.99])).max()) + 0.5
    ax.axline((0,0), slope=1, color="grey", linestyle="--", linewidth=0.8)
    z = np.polyfit(merged["FC_DENV"], merged["FC_ZIKV"], 1)
    xr = np.linspace(-lim, lim, 200)
    ax.plot(xr, np.poly1d(z)(xr), color="black", linewidth=1.2)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("log₂FC (DENV vs Control)"); ax.set_ylabel("log₂FC (ZIKV vs Control)")
    ax.set_title(f"FC Correlation — {label}\nPearson r = {r:.3f}, p = {p:.1e}",
                 fontweight="bold")
    ax.legend(fontsize=9); ax.set_aspect("equal")
    fname = f"FC_correlation_{label.replace(' ','_')}"
    fig.tight_layout()
    fig.savefig(FIG_MAIN / f"{fname}.pdf", bbox_inches="tight")
    fig.savefig(FIG_MAIN / f"{fname}.png", bbox_inches="tight", dpi=150)
    plt.close()
    return r


# ─── Updated shared gene heatmap ──────────────────────────────────────────────
def make_heatmap(shared_ids, denv_res, zikv_res, annot_df, label):
    log(f"Generating heatmap for {label} ...")
    denv_sub = denv_res[denv_res["gene_id"].isin(shared_ids)][
        ["gene_id","log2FoldChange"]].rename(columns={"log2FoldChange":"FC_DENV"})
    zikv_sub = zikv_res[zikv_res["gene_id"].isin(shared_ids)][
        ["gene_id","log2FoldChange"]].rename(columns={"log2FoldChange":"FC_ZIKV"})
    merged = denv_sub.merge(zikv_sub, on="gene_id").merge(
        annot_df[["gene_id","symbol"]], on="gene_id", how="left"
    )
    merged["label"] = merged["symbol"].fillna(merged["gene_id"])
    merged = merged.set_index("label")[["FC_DENV","FC_ZIKV"]].sort_values("FC_DENV", ascending=False)
    merged.columns = ["DENV vs Control","ZIKV vs Control"]

    fig, ax = plt.subplots(figsize=(5, max(4, len(merged)*0.38)))
    sns.heatmap(merged, ax=ax, cmap="RdBu_r", center=0,
                annot=True, fmt=".2f", annot_kws={"size":8},
                linewidths=0.5, linecolor="white",
                cbar_kws={"label":"log₂FC"})
    ax.set_title(f"Shared upregulated genes\n({label})", fontsize=11, fontweight="bold")
    ax.set_ylabel("")
    fig.tight_layout()
    fname = f"SharedGenes_Heatmap_{label.replace(' ','_')}"
    fig.savefig(FIG_MAIN / f"{fname}.pdf", bbox_inches="tight")
    fig.savefig(FIG_MAIN / f"{fname}.png", bbox_inches="tight", dpi=150)
    plt.close()
    log(f"  Heatmap saved → {fname}")


# ─── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    ckpt = load_ckpt()
    log("=" * 60)
    log("Step 03b: Fair DEG comparison — DENV MOI=1 vs ZIKV MOI=1")
    log("=" * 60)

    adata_filt = sc.read_h5ad(ANN_DIR / "adata_processed.h5ad")
    adata_raw  = sc.read_h5ad(ANN_DIR / "adata_raw.h5ad")
    annot_df   = pd.read_csv(BASE_DIR / "02_literature_resources" / "ensembl_annotation.csv")

    # ── A: Pseudobulk MOI=1 ───────────────────────────────────────────────────
    pb_path = PB_DIR / "pseudobulk_counts_moi1.csv"
    if ckpt.get("pb_moi1_done") and pb_path.exists():
        log("✓ MOI=1 pseudobulk already built")
        count_df = pd.read_csv(pb_path, index_col=0)
        info_df  = pd.read_csv(PB_DIR / "sample_info_moi1.csv", index_col=0)
    else:
        count_df, info_df = make_pseudobulk_moi1(adata_raw, adata_filt)
        ckpt["pb_moi1_done"] = True; save_ckpt(ckpt)

    # ── B: DESeq2 DENV MOI=1 ──────────────────────────────────────────────────
    denv_moi1_path = DEG_DIR / "DEGs_DENV_MOI1_vs_Control.csv"
    if ckpt.get("denv_moi1_done") and denv_moi1_path.exists():
        log("✓ DENV MOI=1 DEGs already computed")
        denv_moi1 = pd.read_csv(denv_moi1_path)
    else:
        denv_moi1 = run_deseq2(count_df, info_df, "DENV", "DENV MOI=1 vs Control")
        denv_moi1.to_csv(denv_moi1_path, index=False)
        ckpt["denv_moi1_done"] = True; save_ckpt(ckpt)

    # ── C: Load original ZIKV DEGs (already fair — MOI=1 only) ───────────────
    zikv_res = pd.read_csv(DEG_DIR / "DEGs_ZIKV_vs_Control.csv")

    # ── D: Classify ───────────────────────────────────────────────────────────
    log("\nClassifying DEGs ...")
    denv_moi1_degs = classify(denv_moi1, "DENV MOI=1")
    zikv_degs      = classify(zikv_res,  "ZIKV MOI=1")

    # ── E: Shared DEGs — fair comparison ─────────────────────────────────────
    all_genes = list(set(denv_moi1["gene_id"]) & set(zikv_res["gene_id"]))
    shared_up_fair, shared_down_fair, shared_all_fair, pval_fair, fold_fair = \
        find_shared(denv_moi1_degs, zikv_degs, all_genes, "FAIR MOI=1")

    # Save shared gene lists
    pd.DataFrame({"gene_id": shared_up_fair}).to_csv(
        RES_DIR / "shared_DEGs_up_fair_moi1.csv", index=False)
    pd.DataFrame({"gene_id": shared_down_fair}).to_csv(
        RES_DIR / "shared_DEGs_down_fair_moi1.csv", index=False)

    # ── F: FC correlation — fair ──────────────────────────────────────────────
    if not ckpt.get("fc_fair_done"):
        r_fair = fc_correlation(denv_moi1, zikv_res, shared_up_fair, "Fair MOI=1")
        ckpt["pearson_r_fair"] = float(r_fair); ckpt["fc_fair_done"] = True
        save_ckpt(ckpt)
    else:
        r_fair = ckpt.get("pearson_r_fair", 0)
        log(f"✓ FC correlation done — Pearson r (fair) = {r_fair:.4f}")

    # ── G: Load original results for comparison ───────────────────────────────
    denv_orig   = pd.read_csv(DEG_DIR / "DEGs_DENV_vs_Control.csv")
    orig_ckpt   = json.load(open(CKPT_DIR / "step03_checkpoint.json"))
    denv_o_degs = classify(denv_orig, "DENV Original")
    shared_up_orig = pd.read_csv(RES_DIR / "shared_DEGs_shared_up.csv")["gene_id"].tolist()
    r_orig = orig_ckpt.get("pearson_r", 0)
    fold_orig = orig_ckpt.get("fold_enrichment", 0)

    orig = {"n_denv": len(denv_o_degs["all"]), "n_zikv": len(zikv_degs["all"]),
            "shared_up": len(shared_up_orig), "shared_down": 0,
            "pearson_r": r_orig, "fold_enr": fold_orig}
    fair = {"n_denv": len(denv_moi1_degs["all"]), "n_zikv": len(zikv_degs["all"]),
            "shared_up": len(shared_up_fair), "shared_down": len(shared_down_fair),
            "pearson_r": r_fair, "fold_enr": fold_fair}

    # ── H: Comparison figure ──────────────────────────────────────────────────
    if not ckpt.get("comparison_fig_done"):
        make_comparison_figure(orig, fair)
        ckpt["comparison_fig_done"] = True; save_ckpt(ckpt)

    # ── I: Updated heatmap ────────────────────────────────────────────────────
    if not ckpt.get("heatmap_fair_done"):
        if shared_up_fair:
            make_heatmap(shared_up_fair, denv_moi1, zikv_res, annot_df, "Fair MOI=1")
        if shared_down_fair:
            make_heatmap(shared_down_fair, denv_moi1, zikv_res, annot_df, "Fair MOI=1 Down")
        ckpt["heatmap_fair_done"] = True; save_ckpt(ckpt)

    # ── J: Annotate fair shared genes ─────────────────────────────────────────
    if shared_up_fair:
        shared_ann = denv_moi1[denv_moi1["gene_id"].isin(shared_up_fair)][
            ["gene_id","log2FoldChange","padj"]].rename(
            columns={"log2FoldChange":"FC_DENV","padj":"padj_DENV"})
        shared_ann = shared_ann.merge(
            zikv_res[["gene_id","log2FoldChange","padj"]].rename(
                columns={"log2FoldChange":"FC_ZIKV","padj":"padj_ZIKV"}), on="gene_id")
        shared_ann = shared_ann.merge(annot_df[["gene_id","symbol","name"]], on="gene_id", how="left")
        shared_ann = shared_ann.sort_values("FC_DENV", ascending=False)
        shared_ann.to_csv(RES_DIR / "shared_DEGs_annotated_fair_moi1.csv", index=False)

        log("\n" + "=" * 60)
        log("FAIR COMPARISON — Shared upregulated genes:")
        log("=" * 60)
        for _, row in shared_ann.iterrows():
            log(f"  {str(row.get('symbol','?')):<14} "
                f"FC_DENV={row['FC_DENV']:+.2f}  "
                f"FC_ZIKV={row['FC_ZIKV']:+.2f}  "
                f"{str(row.get('name',''))[:55]}")

    # ── Summary ────────────────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("STEP 03b COMPLETE — Original vs Fair Comparison:")
    log(f"{'Metric':<25} {'Original (MOI=1+10)':>20} {'Fair (MOI=1)':>15}")
    log("-" * 62)
    log(f"{'DENV DEGs':<25} {orig['n_denv']:>20} {fair['n_denv']:>15}")
    log(f"{'ZIKV DEGs':<25} {orig['n_zikv']:>20} {fair['n_zikv']:>15}")
    log(f"{'Shared Up':<25} {orig['shared_up']:>20} {fair['shared_up']:>15}")
    log(f"{'Shared Down':<25} {orig['shared_down']:>20} {fair['shared_down']:>15}")
    log(f"{'Pearson r':<25} {orig['pearson_r']:>20.4f} {fair['pearson_r']:>15.4f}")
    log(f"{'Fold enrichment':<25} {orig['fold_enr']:>20.2f} {fair['fold_enr']:>15.2f}")
    log("=" * 62)

    if fair["pearson_r"] > orig["pearson_r"]:
        log("✓ MOI correction IMPROVES convergence signal — use fair results for paper")
    else:
        log("→ MOI correction did not change signal — original results are valid")
    log("\nNext: run step05_pathway_enrichment.py")


if __name__ == "__main__":
    main()
