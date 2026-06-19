"""
Step 03: Pseudobulk aggregation → DESeq2-style DEG analysis → Shared DEGs
DENV vs Control | ZIKV vs Control | Overlap + Fold-change correlation
Checkpoint-based: safe to restart.
"""

import json, time, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from scipy import sparse
from scipy.stats import pearsonr, spearmanr, fisher_exact
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")

# ─── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR  = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
ANN_DIR   = BASE_DIR / "01_processed_data" / "anndata_objects"
PB_DIR    = BASE_DIR / "01_processed_data" / "pseudobulk"
DEG_DIR   = BASE_DIR / "01_processed_data" / "deg_tables"
RES_DIR   = BASE_DIR / "03_results" / "phase3_shared_degs"
FIG_MAIN  = BASE_DIR / "04_figures" / "main"
FIG_SUPP  = BASE_DIR / "04_figures" / "supplementary"
CKPT_DIR  = BASE_DIR / "checkpoints"
LOG_FILE  = BASE_DIR / "logs" / "step03_deg.log"

for d in [PB_DIR, DEG_DIR, RES_DIR, FIG_MAIN, FIG_SUPP, CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE = CKPT_DIR / "step03_checkpoint.json"

# ─── Thresholds ───────────────────────────────────────────────────────────────
FC_THRESHOLD  = 1.0    # |log2FC| >= 1
P_THRESHOLD   = 0.05   # padj < 0.05

# ─── Helpers ──────────────────────────────────────────────────────────────────
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


# ─── Step A: Pseudobulk aggregation ──────────────────────────────────────────
def make_pseudobulk(adata: sc.AnnData) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregate raw counts per condition × timepoint.
    Returns: (count_matrix df, sample_info df)
    """
    log("Building pseudobulk samples (condition × timepoint) ...")

    # Use raw counts (stored in layers or X before normalization)
    # We saved raw in adata_raw.h5ad; reload it
    adata_raw = sc.read_h5ad(ANN_DIR / "adata_raw.h5ad")

    # Keep only cells that passed QC (match obs index)
    keep_cells = adata.obs_names
    adata_raw  = adata_raw[adata_raw.obs_names.isin(keep_cells)].copy()

    conditions = ["Control", "DENV", "ZIKV"]
    timepoints = ["4h", "12h", "24h", "48h"]

    pb_counts  = {}
    sample_info_rows = []

    for cond in conditions:
        for tp in timepoints:
            mask = (
                (adata_raw.obs["condition"] == cond) &
                (adata_raw.obs["timepoint"] == tp)
            )
            n_cells = mask.sum()
            if n_cells == 0:
                continue

            sample_name = f"{cond}_{tp}"
            # Aggregate: sum raw counts across cells
            X_sub = adata_raw[mask].X
            if sparse.issparse(X_sub):
                agg = np.array(X_sub.sum(axis=0)).flatten()
            else:
                agg = X_sub.sum(axis=0).flatten()

            pb_counts[sample_name] = agg.astype(int)
            sample_info_rows.append({
                "sample":    sample_name,
                "condition": cond,
                "timepoint": tp,
                "n_cells":   int(n_cells)
            })
            log(f"  {sample_name}: {n_cells} cells aggregated, "
                f"total counts = {int(agg.sum()):,}")

    count_df = pd.DataFrame(pb_counts, index=adata_raw.var_names)
    info_df  = pd.DataFrame(sample_info_rows).set_index("sample")

    count_df.to_csv(PB_DIR / "pseudobulk_counts.csv")
    info_df.to_csv(PB_DIR / "sample_info.csv")
    log(f"Pseudobulk matrix: {count_df.shape[1]} samples × {count_df.shape[0]} genes")
    log(f"\n{info_df.to_string()}\n")
    return count_df, info_df


# ─── Step B: PyDESeq2 — combined model (condition + timepoint) ────────────────
def run_deseq2_combined(count_df: pd.DataFrame,
                        info_df: pd.DataFrame,
                        contrast_condition: str) -> pd.DataFrame:
    """
    Run DESeq2 combined model: ~timepoint + condition
    contrast_condition: "DENV" or "ZIKV"
    """
    log(f"\nRunning DESeq2: {contrast_condition} vs Control (combined model) ...")

    # Filter to Control + contrast condition
    keep_samples = info_df[
        info_df["condition"].isin(["Control", contrast_condition])
    ].index.tolist()

    counts_sub = count_df[keep_samples].T  # samples × genes
    info_sub   = info_df.loc[keep_samples].copy()

    # Remove genes with zero counts in all samples
    counts_sub = counts_sub.loc[:, (counts_sub > 0).any(axis=0)]
    # Remove genes with very low total counts
    counts_sub = counts_sub.loc[:, counts_sub.sum(axis=0) >= 10]

    log(f"  Samples: {len(keep_samples)}, Genes after filtering: {counts_sub.shape[1]}")

    # Convert timepoint to ordered factor
    info_sub["timepoint_num"] = info_sub["timepoint"].str.replace("h","").astype(int)

    # PyDESeq2
    dds = DeseqDataSet(
        counts   = counts_sub,
        metadata = info_sub[["condition", "timepoint"]],
        design_factors = ["timepoint", "condition"],
        refit_cooks = True,
        quiet = True
    )
    dds.deseq2()

    stat_res = DeseqStats(
        dds,
        contrast = ["condition", contrast_condition, "Control"],
        quiet = True
    )
    stat_res.summary()

    res_df = stat_res.results_df.copy()
    res_df.index.name = "gene_id"
    res_df = res_df.reset_index()
    res_df["comparison"] = f"{contrast_condition}_vs_Control"

    n_sig = ((res_df["padj"] < P_THRESHOLD) &
             (res_df["log2FoldChange"].abs() >= FC_THRESHOLD)).sum()
    log(f"  Significant DEGs (|log2FC|>={FC_THRESHOLD}, padj<{P_THRESHOLD}): {n_sig}")

    return res_df


# ─── Step C: Classify DEGs ────────────────────────────────────────────────────
def classify_degs(res_df: pd.DataFrame, label: str) -> dict:
    sig = res_df[
        (res_df["padj"] < P_THRESHOLD) &
        (res_df["log2FoldChange"].abs() >= FC_THRESHOLD) &
        res_df["padj"].notna()
    ]
    up   = sig[sig["log2FoldChange"] > 0]["gene_id"].tolist()
    down = sig[sig["log2FoldChange"] < 0]["gene_id"].tolist()
    log(f"  {label} — Up: {len(up)}, Down: {len(down)}, Total: {len(up)+len(down)}")
    return {"up": up, "down": down, "all": up + down}


# ─── Step D: Volcano plots ────────────────────────────────────────────────────
def make_volcano(res_df: pd.DataFrame, title: str,
                 color_up: str, color_down: str,
                 out_path: Path):
    df = res_df.dropna(subset=["padj", "log2FoldChange"]).copy()
    df["-log10padj"] = -np.log10(df["padj"].clip(lower=1e-300))

    sig_up   = (df["padj"] < P_THRESHOLD) & (df["log2FoldChange"] >= FC_THRESHOLD)
    sig_down = (df["padj"] < P_THRESHOLD) & (df["log2FoldChange"] <= -FC_THRESHOLD)
    ns       = ~(sig_up | sig_down)

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(df.loc[ns,  "log2FoldChange"], df.loc[ns,  "-log10padj"],
               c="lightgrey", s=5, alpha=0.5, rasterized=True, label="NS")
    ax.scatter(df.loc[sig_up,  "log2FoldChange"], df.loc[sig_up,  "-log10padj"],
               c=color_up,   s=8, alpha=0.8, rasterized=True,
               label=f"Up ({sig_up.sum()})")
    ax.scatter(df.loc[sig_down,"log2FoldChange"], df.loc[sig_down,"-log10padj"],
               c=color_down, s=8, alpha=0.8, rasterized=True,
               label=f"Down ({sig_down.sum()})")

    ax.axvline(x= FC_THRESHOLD, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(x=-FC_THRESHOLD, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(y=-np.log10(P_THRESHOLD), color="black", linestyle=":", linewidth=0.8)

    ax.set_xlabel("log₂ Fold Change", fontsize=12)
    ax.set_ylabel("-log₁₀(padj)",     fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.legend(framealpha=0.8, fontsize=10)
    ax.set_xlim(-10, 10)

    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close()
    log(f"  Volcano saved → {out_path.name}")


# ─── Step E: Shared DEG analysis ─────────────────────────────────────────────
def find_shared_degs(denv_degs: dict, zikv_degs: dict,
                     all_genes: list) -> dict:
    shared_up   = list(set(denv_degs["up"])   & set(zikv_degs["up"]))
    shared_down = list(set(denv_degs["down"]) & set(zikv_degs["down"]))
    shared_all  = shared_up + shared_down
    discordant  = list(
        (set(denv_degs["up"]) & set(zikv_degs["down"])) |
        (set(denv_degs["down"]) & set(zikv_degs["up"]))
    )
    denv_only = list(set(denv_degs["all"]) - set(zikv_degs["all"]))
    zikv_only = list(set(zikv_degs["all"]) - set(denv_degs["all"]))

    log(f"\nShared DEG Summary:")
    log(f"  Shared up     : {len(shared_up)}")
    log(f"  Shared down   : {len(shared_down)}")
    log(f"  Shared total  : {len(shared_all)}")
    log(f"  Discordant    : {len(discordant)}")
    log(f"  DENV only     : {len(denv_only)}")
    log(f"  ZIKV only     : {len(zikv_only)}")

    # Fisher's exact test (upregulated)
    N  = len(set(all_genes))
    nd = len(denv_degs["up"])
    nz = len(zikv_degs["up"])
    ns = len(shared_up)
    expected = nd * nz / N
    contingency = [
        [ns,          nd - ns],
        [nz - ns,     N - nd - nz + ns]
    ]
    odds, pval = fisher_exact(contingency, alternative="greater")
    fold_enr   = ns / expected if expected > 0 else float("inf")
    log(f"\n  Fisher test (upregulated):")
    log(f"    Expected by chance : {expected:.1f}")
    log(f"    Fold enrichment    : {fold_enr:.2f}x")
    log(f"    Odds ratio         : {odds:.2f}")
    log(f"    p-value            : {pval:.2e}")

    summary = pd.DataFrame([{
        "Category":         cat,
        "Count":            cnt
    } for cat, cnt in [
        ("DENV_up",    len(denv_degs["up"])),
        ("DENV_down",  len(denv_degs["down"])),
        ("ZIKV_up",    len(zikv_degs["up"])),
        ("ZIKV_down",  len(zikv_degs["down"])),
        ("Shared_up",  len(shared_up)),
        ("Shared_down",len(shared_down)),
        ("Shared_all", len(shared_all)),
        ("Discordant", len(discordant)),
        ("DENV_only",  len(denv_only)),
        ("ZIKV_only",  len(zikv_only)),
        ("Fisher_p_up",pval),
        ("FoldEnr_up", round(fold_enr, 3)),
        ("OddsRatio_up", round(odds, 3)),
    ]])
    summary.to_csv(RES_DIR / "deg_summary.csv", index=False)

    return {
        "shared_up":   shared_up,
        "shared_down": shared_down,
        "shared_all":  shared_all,
        "discordant":  discordant,
        "denv_only":   denv_only,
        "zikv_only":   zikv_only,
        "fisher_p":    pval,
        "fold_enr":    fold_enr
    }


# ─── Step F: Fold-change correlation plot ─────────────────────────────────────
def make_fc_correlation(denv_res: pd.DataFrame,
                        zikv_res: pd.DataFrame,
                        shared: dict):
    log("\nGenerating fold-change correlation plot ...")

    merged = pd.merge(
        denv_res[["gene_id","log2FoldChange","padj"]].rename(
            columns={"log2FoldChange":"FC_DENV","padj":"padj_DENV"}),
        zikv_res[["gene_id","log2FoldChange","padj"]].rename(
            columns={"log2FoldChange":"FC_ZIKV","padj":"padj_ZIKV"}),
        on="gene_id"
    ).dropna(subset=["FC_DENV","FC_ZIKV"])

    # Category labels
    merged["category"] = "Other"
    merged.loc[merged["gene_id"].isin(shared["shared_up"]),   "category"] = "Shared_up"
    merged.loc[merged["gene_id"].isin(shared["shared_down"]), "category"] = "Shared_down"

    pearson_r, pearson_p = pearsonr(merged["FC_DENV"], merged["FC_ZIKV"])
    spearman_r, spearman_p = spearmanr(merged["FC_DENV"], merged["FC_ZIKV"])
    log(f"  Pearson r  = {pearson_r:.4f}, p = {pearson_p:.2e}")
    log(f"  Spearman r = {spearman_r:.4f}, p = {spearman_p:.2e}")

    merged.to_csv(RES_DIR / "merged_foldchanges.csv", index=False)
    pd.DataFrame([{
        "Pearson_r":pearson_r, "Pearson_p":pearson_p,
        "Spearman_r":spearman_r,"Spearman_p":spearman_p
    }]).to_csv(RES_DIR / "correlation_results.csv", index=False)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 7))

    cat_style = {
        "Other":      {"c":"lightgrey", "s":4,  "alpha":0.3, "zorder":1},
        "Shared_up":  {"c":"#D73027",   "s":20, "alpha":0.8, "zorder":3},
        "Shared_down":{"c":"#4575B4",   "s":20, "alpha":0.8, "zorder":3},
    }
    handles = []
    for cat, style in cat_style.items():
        sub = merged[merged["category"] == cat]
        ax.scatter(sub["FC_DENV"], sub["FC_ZIKV"],
                   label=f"{cat} (n={len(sub)})", **style, rasterized=True)

    # Reference lines
    lim = max(abs(merged["FC_DENV"].quantile(0.99)),
              abs(merged["FC_ZIKV"].quantile(0.99))) + 0.5
    ax.axline((0,0), slope=1, color="grey", linestyle="--",
              linewidth=0.8, label="y=x (perfect conservation)")
    ax.axhline(0, color="black", linewidth=0.4)
    ax.axvline(0, color="black", linewidth=0.4)

    # Regression line
    z = np.polyfit(merged["FC_DENV"], merged["FC_ZIKV"], 1)
    xr = np.linspace(-lim, lim, 200)
    ax.plot(xr, np.poly1d(z)(xr), color="black", linewidth=1.2, label="Regression")

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log₂FC (DENV vs Control)", fontsize=12)
    ax.set_ylabel("log₂FC (ZIKV vs Control)", fontsize=12)
    ax.set_title(
        f"Genome-wide fold-change correlation\n"
        f"Pearson r = {pearson_r:.3f}, p = {pearson_p:.1e}",
        fontsize=12, fontweight="bold"
    )
    ax.legend(fontsize=9, framealpha=0.8)
    ax.set_aspect("equal")

    fig.tight_layout()
    fig.savefig(FIG_MAIN / "Figure_FoldChange_Correlation.pdf", bbox_inches="tight")
    fig.savefig(FIG_MAIN / "Figure_FoldChange_Correlation.png", bbox_inches="tight", dpi=150)
    plt.close()
    log("  Fold-change correlation plot saved → FIG_MAIN")

    return merged, pearson_r


# ─── Step G: Venn-style overlap bar chart ─────────────────────────────────────
def make_overlap_figure(denv_degs: dict, zikv_degs: dict, shared: dict):
    log("Generating DEG overlap figure ...")

    categories = ["DENV only", "Shared\n(concordant)", "ZIKV only", "Discordant"]
    up_vals    = [
        len(set(denv_degs["up"]) - set(zikv_degs["up"])),
        len(shared["shared_up"]),
        len(set(zikv_degs["up"]) - set(denv_degs["up"])),
        len(set(denv_degs["up"]) & set(zikv_degs["down"]))
    ]
    down_vals  = [
        len(set(denv_degs["down"]) - set(zikv_degs["down"])),
        len(shared["shared_down"]),
        len(set(zikv_degs["down"]) - set(denv_degs["down"])),
        len(set(denv_degs["down"]) & set(zikv_degs["up"]))
    ]

    x     = np.arange(len(categories))
    width = 0.35
    fig, ax = plt.subplots(figsize=(10, 6))

    bars1 = ax.bar(x - width/2, up_vals,   width, label="Upregulated",
                   color=["#E41A1C","#D73027","#377EB8","#FF7F00"], alpha=0.8)
    bars2 = ax.bar(x + width/2, down_vals, width, label="Downregulated",
                   color=["#E41A1C","#4575B4","#377EB8","#FF7F00"], alpha=0.5,
                   hatch="///")

    for bar in list(bars1) + list(bars2):
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width()/2, h + 2,
                    str(int(h)), ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x); ax.set_xticklabels(categories, fontsize=11)
    ax.set_ylabel("Number of genes", fontsize=12)
    ax.set_title("DEG overlap: DENV vs ZIKV (vs Control)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_MAIN / "Figure_DEG_overlap.pdf", bbox_inches="tight")
    fig.savefig(FIG_MAIN / "Figure_DEG_overlap.png", bbox_inches="tight", dpi=150)
    plt.close()
    log("  DEG overlap figure saved → FIG_MAIN")


# ─── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    ckpt = load_ckpt()
    log("=" * 60)
    log("Step 03: Pseudobulk DEG → Shared DEGs → Correlation")
    log("=" * 60)

    # Load processed AnnData
    log("Loading processed AnnData ...")
    adata = sc.read_h5ad(ANN_DIR / "adata_processed.h5ad")
    log(f"  Shape: {adata.shape}")

    # ── A: Pseudobulk ─────────────────────────────────────────────────────────
    pb_path   = PB_DIR / "pseudobulk_counts.csv"
    info_path = PB_DIR / "sample_info.csv"
    if ckpt.get("pseudobulk_done") and pb_path.exists():
        log("✓ Pseudobulk already built — loading ...")
        count_df = pd.read_csv(pb_path, index_col=0)
        info_df  = pd.read_csv(info_path, index_col=0)
        log(f"  {count_df.shape[1]} samples × {count_df.shape[0]} genes")
    else:
        count_df, info_df = make_pseudobulk(adata)
        ckpt["pseudobulk_done"] = True
        save_ckpt(ckpt)

    # ── B: DESeq2 — DENV vs Control ───────────────────────────────────────────
    denv_path = DEG_DIR / "DEGs_DENV_vs_Control.csv"
    if ckpt.get("denv_deg_done") and denv_path.exists():
        log("✓ DENV DEGs already computed — loading ...")
        denv_res = pd.read_csv(denv_path)
    else:
        denv_res = run_deseq2_combined(count_df, info_df, "DENV")
        denv_res.to_csv(denv_path, index=False)
        ckpt["denv_deg_done"] = True
        save_ckpt(ckpt)

    # ── C: DESeq2 — ZIKV vs Control ───────────────────────────────────────────
    zikv_path = DEG_DIR / "DEGs_ZIKV_vs_Control.csv"
    if ckpt.get("zikv_deg_done") and zikv_path.exists():
        log("✓ ZIKV DEGs already computed — loading ...")
        zikv_res = pd.read_csv(zikv_path)
    else:
        zikv_res = run_deseq2_combined(count_df, info_df, "ZIKV")
        zikv_res.to_csv(zikv_path, index=False)
        ckpt["zikv_deg_done"] = True
        save_ckpt(ckpt)

    # ── D: Volcano plots ───────────────────────────────────────────────────────
    if not ckpt.get("volcano_done"):
        make_volcano(denv_res, "DENV vs Control",
                     "#E41A1C", "#4575B4",
                     FIG_SUPP / "volcano_DENV")
        make_volcano(zikv_res, "ZIKV vs Control",
                     "#E41A1C", "#4575B4",
                     FIG_SUPP / "volcano_ZIKV")
        ckpt["volcano_done"] = True
        save_ckpt(ckpt)
    else:
        log("✓ Volcano plots already saved")

    # ── E: Classify DEGs ──────────────────────────────────────────────────────
    log("\nClassifying DEGs ...")
    all_genes  = list(set(denv_res["gene_id"].tolist()) &
                      set(zikv_res["gene_id"].tolist()))
    denv_degs  = classify_degs(denv_res, "DENV")
    zikv_degs  = classify_degs(zikv_res, "ZIKV")

    # ── F: Shared DEG analysis ─────────────────────────────────────────────────
    if not ckpt.get("shared_done"):
        shared = find_shared_degs(denv_degs, zikv_degs, all_genes)

        # Save gene lists
        for key, genes in shared.items():
            if isinstance(genes, list):
                pd.DataFrame({"gene_id": genes}).to_csv(
                    RES_DIR / f"shared_DEGs_{key}.csv", index=False)

        ckpt["shared_done"]     = True
        ckpt["n_shared_up"]     = len(shared["shared_up"])
        ckpt["n_shared_down"]   = len(shared["shared_down"])
        ckpt["fisher_p"]        = float(shared["fisher_p"])
        ckpt["fold_enrichment"] = float(shared["fold_enr"])
        save_ckpt(ckpt)
    else:
        log("✓ Shared DEG analysis already done — loading ...")
        shared = {
            "shared_up":  pd.read_csv(RES_DIR/"shared_DEGs_shared_up.csv")["gene_id"].tolist(),
            "shared_down":pd.read_csv(RES_DIR/"shared_DEGs_shared_down.csv")["gene_id"].tolist(),
        }
        shared["shared_all"] = shared["shared_up"] + shared["shared_down"]
        denv_only = list(set(denv_degs["all"]) - set(shared["shared_all"]))
        zikv_only = list(set(zikv_degs["all"]) - set(shared["shared_all"]))
        shared["denv_only"]   = denv_only
        shared["zikv_only"]   = zikv_only
        shared["discordant"]  = list(
            (set(denv_degs["up"]) & set(zikv_degs["down"])) |
            (set(denv_degs["down"]) & set(zikv_degs["up"]))
        )

    # ── G: Fold-change correlation ─────────────────────────────────────────────
    if not ckpt.get("fc_corr_done"):
        merged, pearson_r = make_fc_correlation(denv_res, zikv_res, shared)
        ckpt["fc_corr_done"] = True
        ckpt["pearson_r"]    = float(pearson_r)
        save_ckpt(ckpt)
        log(f"\n  GATE G3 — Pearson r = {pearson_r:.3f}")
        if pearson_r > 0.5:
            log("  ✓ STRONG SIGNAL — proceed with full confidence")
        elif pearson_r > 0.3:
            log("  ⚠ MODERATE SIGNAL — proceed, note as moderate convergence")
        else:
            log("  ✗ WEAK/NO SIGNAL — investigate before proceeding")
    else:
        log(f"✓ FC correlation done — Pearson r = {ckpt.get('pearson_r','?')}")

    # ── H: Overlap figure ─────────────────────────────────────────────────────
    if not ckpt.get("overlap_fig_done"):
        make_overlap_figure(denv_degs, zikv_degs, shared)
        ckpt["overlap_fig_done"] = True
        save_ckpt(ckpt)
    else:
        log("✓ Overlap figure already saved")

    # ── Summary ────────────────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("STEP 03 COMPLETE — Summary:")
    log(f"  DENV DEGs (up/down): {len(denv_degs['up'])} / {len(denv_degs['down'])}")
    log(f"  ZIKV DEGs (up/down): {len(zikv_degs['up'])} / {len(zikv_degs['down'])}")
    log(f"  Shared up           : {len(shared['shared_up'])}")
    log(f"  Shared down         : {len(shared['shared_down'])}")
    log(f"  Shared total        : {len(shared['shared_all'])}")
    log(f"  Pearson r           : {ckpt.get('pearson_r','?')}")
    log(f"  Fisher p (up)       : {ckpt.get('fisher_p','?')}")
    log(f"  Fold enrichment     : {ckpt.get('fold_enrichment','?')}")
    log(f"\n  DEG tables  → {DEG_DIR}")
    log(f"  Shared lists→ {RES_DIR}")
    log(f"  Figures     → {FIG_MAIN}")
    log("\nNext: run step04_gene_annotation.py")
    log("=" * 60)


if __name__ == "__main__":
    main()
