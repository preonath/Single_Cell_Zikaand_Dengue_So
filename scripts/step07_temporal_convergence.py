"""
Step 07: Temporal Convergence Trajectory (SOP Phase 4, Step 4.4)
Computes DENV vs ZIKV log2FC Pearson correlation at each timepoint
(4h, 12h, 24h, 48h) using normalized pseudobulk counts.
One pseudobulk sample per condition×timepoint → simple log2FC ratio method.
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
from scipy.stats import pearsonr, spearmanr

warnings.filterwarnings("ignore")

# ─── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
PB_DIR   = BASE_DIR / "01_processed_data" / "pseudobulk"
DEG_DIR  = BASE_DIR / "01_processed_data" / "deg_tables"
RES_DIR  = BASE_DIR / "03_results" / "phase4_temporal"
FIG_MAIN = BASE_DIR / "04_figures" / "main"
CKPT_DIR = BASE_DIR / "checkpoints"
LOG_FILE = BASE_DIR / "logs" / "step07_temporal.log"

for d in [RES_DIR, FIG_MAIN, CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE  = CKPT_DIR / "step07_checkpoint.json"
TIMEPOINTS = ["4h", "12h", "24h", "48h"]
PSEUDOCOUNT = 1.0

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


def main():
    ckpt = load_ckpt()

    log("=" * 60)
    log("Step 07: Temporal Convergence Trajectory (Phase 4, Step 4.4)")
    log("=" * 60)

    # ─── Load pseudobulk counts ───────────────────────────────────────────────
    log("Loading pseudobulk counts matrix ...")
    counts = pd.read_csv(PB_DIR / "pseudobulk_counts.csv", index_col=0)
    log(f"  Shape: {counts.shape[0]} genes × {counts.shape[1]} samples")
    log(f"  Columns: {list(counts.columns)}")

    # ─── Normalize to CPM (counts per million) ────────────────────────────────
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    log(f"  Normalized to CPM")

    # Filter low-expression genes (CPM > 1 in at least 3 samples)
    keep = (cpm > 1).sum(axis=1) >= 3
    cpm = cpm[keep]
    log(f"  Genes after CPM>1 filter: {cpm.shape[0]}")

    # Log2 transform
    lcpm = np.log2(cpm + PSEUDOCOUNT)

    # ─── Load overall DEG results for context ─────────────────────────────────
    denv_degs = pd.read_csv(DEG_DIR / "DEGs_DENV_vs_Control.csv", index_col="gene_id")
    zikv_degs = pd.read_csv(DEG_DIR / "DEGs_ZIKV_vs_Control.csv", index_col="gene_id")

    shared_genes_file = BASE_DIR / "03_results" / "phase3_shared_degs" / "shared_DEGs_all.csv"
    if shared_genes_file.exists():
        shared_df = pd.read_csv(shared_genes_file)
        shared_genes = set(shared_df["gene_id"].tolist()) if "gene_id" in shared_df.columns else set()
    else:
        shared_genes = set()
    log(f"  Shared DEGs loaded: {len(shared_genes)} genes")

    # ─── Compute log2FC at each timepoint ─────────────────────────────────────
    log("\nComputing log2FC (DENV/Control) and (ZIKV/Control) per timepoint ...")

    tp_results = []

    for tp in TIMEPOINTS:
        ctrl_col = f"Control_{tp}"
        denv_col = f"DENV_{tp}"
        zikv_col = f"ZIKV_{tp}"

        if not all(c in lcpm.columns for c in [ctrl_col, denv_col, zikv_col]):
            log(f"  {tp}: Missing columns — skipping")
            continue

        lfc_denv = lcpm[denv_col] - lcpm[ctrl_col]
        lfc_zikv = lcpm[zikv_col] - lcpm[ctrl_col]

        # Align and keep finite values
        fc_df = pd.DataFrame({"lfc_DENV": lfc_denv, "lfc_ZIKV": lfc_zikv}).dropna()
        fc_df = fc_df[np.isfinite(fc_df["lfc_DENV"]) & np.isfinite(fc_df["lfc_ZIKV"])]

        r, p = pearsonr(fc_df["lfc_DENV"], fc_df["lfc_ZIKV"])
        rs, ps = spearmanr(fc_df["lfc_DENV"], fc_df["lfc_ZIKV"])

        # Count genes moving in same direction (concordant)
        concordant_up = ((fc_df["lfc_DENV"] > 0.5) & (fc_df["lfc_ZIKV"] > 0.5)).sum()
        concordant_dn = ((fc_df["lfc_DENV"] < -0.5) & (fc_df["lfc_ZIKV"] < -0.5)).sum()
        discordant    = ((fc_df["lfc_DENV"] > 0.5) & (fc_df["lfc_ZIKV"] < -0.5)).sum() + \
                        ((fc_df["lfc_DENV"] < -0.5) & (fc_df["lfc_ZIKV"] > 0.5)).sum()

        # Check shared genes specifically
        shared_in_tp = fc_df[fc_df.index.isin(shared_genes)]
        if len(shared_in_tp) > 1:
            r_shared, p_shared = pearsonr(shared_in_tp["lfc_DENV"], shared_in_tp["lfc_ZIKV"])
        else:
            r_shared, p_shared = np.nan, np.nan

        log(f"  {tp}: r={r:.4f} (p={p:.2e})  Spearman_r={rs:.4f}  "
            f"n={len(fc_df)}  concordant_up={concordant_up}  concordant_dn={concordant_dn}  "
            f"discordant={discordant}  r_shared={r_shared:.4f}")

        tp_results.append({
            "timepoint": tp,
            "timepoint_h": int(tp.replace("h", "")),
            "pearson_r": round(r, 4),
            "pearson_p": p,
            "spearman_r": round(rs, 4),
            "n_genes": len(fc_df),
            "concordant_up": int(concordant_up),
            "concordant_down": int(concordant_dn),
            "discordant": int(discordant),
            "r_shared_genes": round(r_shared, 4) if not np.isnan(r_shared) else None,
        })

        # Save per-timepoint log2FC table
        fc_df.to_csv(RES_DIR / f"log2FC_per_gene_{tp}.csv")

    summary_df = pd.DataFrame(tp_results)
    summary_df.to_csv(RES_DIR / "temporal_convergence_summary.csv", index=False)
    log(f"\nSummary saved → {RES_DIR}/temporal_convergence_summary.csv")

    # ─── Figure: Temporal convergence trajectory ───────────────────────────────
    log("\nGenerating figures ...")

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    tp_labels = summary_df["timepoint"].tolist()
    rs_all    = summary_df["pearson_r"].tolist()
    conc_up   = summary_df["concordant_up"].tolist()
    conc_dn   = summary_df["concordant_down"].tolist()

    # Panel A — Pearson r trajectory line + bars
    ax = axes[0]
    bar_colors = ["#1565C0" if r >= 0.4 else "#FB8C00" if r >= 0.25 else "#C62828" for r in rs_all]
    bars = ax.bar(tp_labels, rs_all, color=bar_colors, edgecolor="black", linewidth=0.8, width=0.5, alpha=0.85)
    ax.plot(tp_labels, rs_all, "ko-", linewidth=1.5, markersize=5, zorder=5)
    ax.axhline(0.4, color="green", linestyle="--", linewidth=1.3, label="G3 threshold (r=0.4)")
    ax.axhline(0.0, color="gray", linewidth=0.5)
    for bar, r in zip(bars, rs_all):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                f"{r:.3f}", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xlabel("Timepoint post-infection", fontsize=11)
    ax.set_ylabel("Pearson r (DENV vs ZIKV log2FC)", fontsize=11)
    ax.set_title("A  FC Correlation Trajectory", fontsize=12, fontweight="bold")
    ax.set_ylim(min(min(rs_all) - 0.05, -0.05), max(rs_all) + 0.12)
    ax.legend(fontsize=9)

    # Panel B — Concordant gene counts over time
    ax2 = axes[1]
    x = range(len(tp_labels))
    ax2.bar([i - 0.18 for i in x], conc_up, width=0.35, label="Concordant UP (both ↑)", color="#E53935", alpha=0.85, edgecolor="black")
    ax2.bar([i + 0.18 for i in x], conc_dn, width=0.35, label="Concordant DOWN (both ↓)", color="#1E88E5", alpha=0.85, edgecolor="black")
    ax2.set_xticks(list(x))
    ax2.set_xticklabels(tp_labels)
    ax2.set_xlabel("Timepoint post-infection", fontsize=11)
    ax2.set_ylabel("# Concordant Genes (|log2FC| > 0.5)", fontsize=11)
    ax2.set_title("B  Concordant Gene Counts", fontsize=12, fontweight="bold")
    ax2.legend(fontsize=9)

    # Panel C — FC scatter at the timepoint with highest r
    best_tp = tp_labels[rs_all.index(max(rs_all))]
    fc_best = pd.read_csv(RES_DIR / f"log2FC_per_gene_{best_tp}.csv", index_col=0)
    ax3 = axes[2]
    ax3.scatter(fc_best["lfc_DENV"], fc_best["lfc_ZIKV"],
                alpha=0.08, s=4, color="#9E9E9E", rasterized=True)
    # Highlight shared genes
    shared_in_best = fc_best[fc_best.index.isin(shared_genes)]
    if len(shared_in_best) > 0:
        ax3.scatter(shared_in_best["lfc_DENV"], shared_in_best["lfc_ZIKV"],
                    color="#E91E63", s=50, zorder=5, label=f"Shared DEGs (n={len(shared_in_best)})", edgecolors="black", linewidths=0.5)
    lim = max(abs(fc_best["lfc_DENV"].quantile(0.99)), abs(fc_best["lfc_ZIKV"].quantile(0.99))) * 1.1
    ax3.set_xlim(-lim, lim)
    ax3.set_ylim(-lim, lim)
    ax3.axhline(0, color="gray", linewidth=0.5)
    ax3.axvline(0, color="gray", linewidth=0.5)
    best_r = rs_all[tp_labels.index(best_tp)]
    ax3.text(0.05, 0.95, f"r = {best_r:.4f}", transform=ax3.transAxes,
             fontsize=11, fontweight="bold", va="top",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray"))
    ax3.set_xlabel(f"DENV log2FC ({best_tp})", fontsize=11)
    ax3.set_ylabel(f"ZIKV log2FC ({best_tp})", fontsize=11)
    ax3.set_title(f"C  FC Scatter at {best_tp} (peak r)", fontsize=12, fontweight="bold")
    ax3.legend(fontsize=9, loc="lower right")

    plt.suptitle("DENV–ZIKV Temporal FC Convergence — GSE110496", fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig_path = FIG_MAIN / "Figure_Temporal_Convergence.png"
    plt.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close()
    log(f"Figure saved → {fig_path}")

    # ─── Print final summary ──────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("STEP 07 COMPLETE — Temporal Convergence Summary")
    log("=" * 60)
    log(f"{'Timepoint':<10} {'Pearson r':>10} {'Spearman r':>12} {'Conc.UP':>9} {'Conc.DN':>9} {'G3 Status':>12}")
    log("-" * 65)
    for row in tp_results:
        g3 = "PASS" if row["pearson_r"] >= 0.4 else "MODERATE" if row["pearson_r"] >= 0.25 else "WEAK"
        log(f"  {row['timepoint']:<8} {row['pearson_r']:>10.4f} {row['spearman_r']:>12.4f} "
            f"{row['concordant_up']:>9} {row['concordant_down']:>9} {g3:>12}")

    max_r = max(rs_all)
    max_tp = tp_labels[rs_all.index(max_r)]
    log(f"\nPeak convergence: r = {max_r:.4f} at {max_tp}")
    log(f"Overall (all timepoints collapsed, from Step 03): r = 0.3569")

    ckpt["temporal_done"] = True
    ckpt["peak_r"] = max_r
    ckpt["peak_tp"] = max_tp
    ckpt["temporal_summary"] = {row["timepoint"]: row["pearson_r"] for row in tp_results}
    save_ckpt(ckpt)
    log("\nNext: run step08_prepare_literature_resources.py")


if __name__ == "__main__":
    main()
