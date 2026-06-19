"""
Step 15: GSE118305 DEA — ZIKV Macrophages (SOP Phase 8, Step 8.1)
Compares ZIKV-infected (4G2pos) vs Mock at 24h using log2-FPKM + Welch t-test.
Adds GSE118305 to the GATE G6 replication analysis.
"""

import json, time, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

BASE_DIR = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
PROC_DIR = BASE_DIR / "01_processed_data" / "deg_tables"
RES_DIR  = BASE_DIR / "03_results" / "phase8_validation"
FIG_MAIN = BASE_DIR / "04_figures" / "main"
LOG_FILE = BASE_DIR / "logs" / "step15_gse118305.log"

for d in [RES_DIR, FIG_MAIN, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

def main():
    log("=" * 60)
    log("Step 15: GSE118305 DEA — ZIKV Macrophages (Phase 8.1)")
    log("=" * 60)

    # ─── Load FPKM matrix ────────────────────────────────────────────────────
    log("Loading FPKM matrix ...")
    expr = pd.read_csv(PROC_DIR / "GSE118305_fpkm_genes.csv", index_col=0)
    log(f"  Shape: {expr.shape}")

    # ─── Identify RNA-seq sample columns ─────────────────────────────────────
    # Column naming: "HMDM-RNAseq-{condition}-{time}-{4G2pos/neg}-{donor} FPKM"
    all_cols = expr.columns.tolist()
    rna_cols = [c for c in all_cols if "RNAseq" in c]
    log(f"  RNA-seq columns: {len(rna_cols)}")

    mock_cols   = [c for c in rna_cols if "mock" in c.lower() and "pooled" not in c.lower()]
    zikv_cols   = [c for c in rna_cols if "4G2pos" in c and "pooled" not in c and "24h" in c]
    bystand_cols= [c for c in rna_cols if "4G2neg" in c and "pooled" not in c and "24h" in c]

    log(f"  Mock 24h:         {len(mock_cols)} samples → {[c.split()[0] for c in mock_cols]}")
    log(f"  ZIKV-infected 24h:{len(zikv_cols)} samples → {[c.split()[0] for c in zikv_cols]}")
    log(f"  ZIKV-bystander 24h:{len(bystand_cols)} samples")

    # ─── Log2-transform (FPKM → log2(FPKM+1)) ───────────────────────────────
    log2_expr = np.log2(expr + 1)

    # Filter: expressed in at least 3 samples (log2FPKM > 0.5)
    keep = ((log2_expr[mock_cols] > 0.5).sum(axis=1) >= 2) | \
           ((log2_expr[zikv_cols] > 0.5).sum(axis=1) >= 2)
    log2_expr = log2_expr[keep]
    log(f"  Genes after expression filter: {log2_expr.shape[0]}")

    # ─── Welch t-test: ZIKV-infected vs Mock at 24h ──────────────────────────
    log("Running Welch t-test: ZIKV-infected vs Mock (24h) ...")
    results = []
    for gene in log2_expr.index:
        mock_vals = log2_expr.loc[gene, mock_cols].dropna().values
        zikv_vals = log2_expr.loc[gene, zikv_cols].dropna().values
        if len(mock_vals) < 2 or len(zikv_vals) < 2:
            continue
        lfc = np.mean(zikv_vals) - np.mean(mock_vals)
        stat, pval = ttest_ind(zikv_vals, mock_vals, equal_var=False)
        results.append({"symbol": gene, "log2FoldChange": lfc, "pvalue": pval,
                        "mean_zikv": np.mean(zikv_vals), "mean_mock": np.mean(mock_vals)})

    res_df = pd.DataFrame(results).dropna(subset=["pvalue"])

    # BH FDR correction
    _, padj, _, _ = multipletests(res_df["pvalue"], method="fdr_bh")
    res_df["padj"] = padj
    res_df = res_df.sort_values("pvalue")

    degs = res_df[(res_df["padj"] < 0.05) & (res_df["log2FoldChange"].abs() >= 1)]
    degs_up = degs[degs["log2FoldChange"] > 0]
    degs_dn = degs[degs["log2FoldChange"] < 0]

    log(f"  Total genes tested: {len(res_df)}")
    log(f"  Significant DEGs (padj<0.05, |lFC|≥1): {len(degs)}")
    log(f"  Upregulated: {len(degs_up)}  Downregulated: {len(degs_dn)}")
    log(f"  Top upregulated: {degs_up.head(10)['symbol'].tolist()}")

    res_df.to_csv(PROC_DIR / "DEGs_ZIKV_macrophages_GSE118305.csv", index=False)
    log(f"  Saved: DEGs_ZIKV_macrophages_GSE118305.csv")

    # ─── GATE G6 update: check discovery shared genes ────────────────────────
    shared_ann = pd.read_csv(BASE_DIR / "03_results" / "phase3_shared_degs" / "shared_DEGs_annotated.csv")
    shared_up_genes = set(shared_ann[shared_ann["log2FC_DENV"] > 0]["symbol"].str.upper().tolist())
    val_up_genes    = set(degs_up["symbol"].str.upper().tolist())
    bg_genes        = set(res_df["symbol"].str.upper().tolist())

    overlap = shared_up_genes & val_up_genes
    log(f"\n  Discovery genes replicated in GSE118305 macrophages: {len(overlap)}")
    log(f"  Replicated: {sorted(overlap)}")

    from scipy.stats import fisher_exact
    a = len(overlap)
    b = len(shared_up_genes - val_up_genes)
    c = len((bg_genes - shared_up_genes) & val_up_genes)
    d = len((bg_genes - shared_up_genes) - val_up_genes)
    if a > 0:
        or_v, pval = fisher_exact([[a,b],[c,d]], alternative="greater")
        fe = (a / len(shared_up_genes)) / (len(val_up_genes) / max(len(bg_genes),1))
    else:
        pval, fe = 1.0, 0.0

    rep_rate = round(100 * a / max(len(shared_up_genes), 1), 1)
    log(f"  Replication rate: {rep_rate}%  Fisher p = {pval:.4e}  FE = {fe:.2f}×")

    # ─── Volcano plot ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 6))
    res_plot = res_df.copy()
    res_plot["-log10p"] = -np.log10(res_plot["pvalue"].clip(lower=1e-30))
    colors = np.where((res_plot["padj"] < 0.05) & (res_plot["log2FoldChange"] > 1), "#D32F2F",
             np.where((res_plot["padj"] < 0.05) & (res_plot["log2FoldChange"] < -1), "#1565C0", "#BDBDBD"))
    ax.scatter(res_plot["log2FoldChange"], res_plot["-log10p"], c=colors, s=4, alpha=0.5, rasterized=True)
    ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axvline(1, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axvline(-1, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)

    # Highlight shared genes
    for gene in shared_up_genes:
        row = res_plot[res_plot["symbol"] == gene]
        if len(row) > 0:
            ax.scatter(row["log2FoldChange"].values[0], row["-log10p"].values[0],
                       color="#FF6F00", s=80, zorder=5, edgecolors="black", linewidths=0.5)
            ax.annotate(gene, (row["log2FoldChange"].values[0], row["-log10p"].values[0]),
                        fontsize=7, xytext=(4, 2), textcoords="offset points")

    ax.set_xlabel("log2 Fold Change (ZIKV-infected vs Mock)", fontsize=11)
    ax.set_ylabel("-log10(p-value)", fontsize=11)
    ax.set_title(f"ZIKV Macrophages GSE118305\nVolcano Plot (24h, {len(degs)} DEGs, {rep_rate}% shared genes replicated)",
                 fontsize=11, fontweight="bold")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#D32F2F", label=f"Up ({len(degs_up)})"),
                        Patch(color="#1565C0", label=f"Down ({len(degs_dn)})"),
                        Patch(color="#FF6F00", label=f"Discovery shared ({len(overlap)} replicated)")],
              fontsize=9)
    plt.tight_layout()
    plt.savefig(FIG_MAIN / "Figure_GSE118305_ZIKV_Macrophages_Volcano.png", dpi=200, bbox_inches="tight")
    plt.close()
    log(f"  Volcano saved → Figure_GSE118305_ZIKV_Macrophages_Volcano.png")

    log(f"\nSTEP 15 COMPLETE")

if __name__ == "__main__":
    main()
