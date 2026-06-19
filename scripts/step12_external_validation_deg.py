"""
Step 12: External Validation DEG Analysis (SOP Phase 8, Steps 8.1-8.4 + GATE G6)
Processes pre-downloaded expression data from three validation datasets:
  - GSE118305 (ZIKV Macrophages — FPKM matrix, needs DEA)
  - GSE94892  (DENV PBMCs — pre-computed DEG table from Cufflinks)
  - GSE78711  (ZIKV NPCs — pre-computed fold change table)

Then runs GATE G6: Fisher test for replication of discovery shared DEGs.
Checkpoint-based: safe to restart.
"""

import json, time, gzip, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import fisher_exact, pearsonr

warnings.filterwarnings("ignore")

BASE_DIR   = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
RAW_DIR    = BASE_DIR / "00_raw_data"
PROC_DIR   = BASE_DIR / "01_processed_data" / "deg_tables"
RES_DIR    = BASE_DIR / "03_results" / "phase8_validation"
FIG_MAIN   = BASE_DIR / "04_figures" / "main"
CKPT_DIR   = BASE_DIR / "checkpoints"
LOG_FILE   = BASE_DIR / "logs" / "step12_validation.log"

for d in [PROC_DIR, RES_DIR, FIG_MAIN, CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE  = CKPT_DIR / "step12_checkpoint.json"
FC_THRESH  = 1.0
P_THRESH   = 0.05

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


# ─── Parse GSE118305 (ZIKV Macrophages, FPKM matrix) ─────────────────────────
def parse_gse118305():
    """Parse HOMER FPKM output, aggregate to gene level, get ZIKV vs Mock DEGs."""
    log("Parsing GSE118305 (ZIKV Macrophages FPKM) ...")

    fpath = RAW_DIR / "GSE118305" / "GSE118305_expression.txt.gz"
    with gzip.open(fpath, "rt") as f:
        lines = f.readlines()

    # Header line 1 (comment with column info), skip it; header line 2 is actual columns
    # Find the data start
    header_line = None
    skip = 0
    for i, line in enumerate(lines):
        if line.startswith("Transcript"):
            header_line = i
            break

    df = pd.read_csv(fpath, sep="\t", skiprows=header_line, comment=None, low_memory=False)
    log(f"  Raw shape: {df.shape}, Columns (first 15): {list(df.columns[:15])}")

    # Column 1 is TranscriptID, Column 7 is gene name, rest are sample FPKMs
    # Find gene name column
    cols = df.columns.tolist()
    gene_col = None
    for c in cols:
        if "Gene Name" in c or "gene" in c.lower() and "id" not in c.lower():
            gene_col = c
            break
    if gene_col is None:
        gene_col = cols[6]  # typically 7th column is gene name

    # Sample columns are FPKM values — typically from column 9 onwards
    # Find FPKM columns (numeric, column names contain sample IDs)
    id_cols = [cols[0], gene_col]
    potential_id_cols = [c for c in cols[:9] if not c.startswith("Unnamed")]
    sample_cols = [c for c in cols if c not in potential_id_cols[:8] and
                   df[c].dtype in [np.float64, np.int64] and df[c].notna().sum() > 100]

    if len(sample_cols) == 0:
        # Try differently: look for columns starting with numeric-like names
        sample_cols = [c for c in cols[8:] if not c.startswith("Unnamed")]

    log(f"  Gene column: '{gene_col}'")
    log(f"  Sample columns ({len(sample_cols)}): {sample_cols[:10]}")

    # Aggregate to gene level (sum FPKM by gene)
    df_genes = df[[gene_col] + sample_cols].copy()
    df_genes.columns = ["gene"] + sample_cols
    df_genes = df_genes[df_genes["gene"].notna() & (df_genes["gene"] != "")]
    df_genes = df_genes.groupby("gene")[sample_cols].sum()
    log(f"  Gene-level matrix: {df_genes.shape}")

    # Get metadata to identify ZIKV vs Mock samples
    meta = pd.read_csv(RAW_DIR / "GSE118305" / "sample_metadata.csv")
    log(f"  Metadata columns: {list(meta.columns)}")
    log(f"  Unique viral_infection values: {meta['viral_infection'].unique() if 'viral_infection' in meta.columns else 'N/A'}")

    # Return what we have
    df_genes.to_csv(PROC_DIR / "GSE118305_fpkm_genes.csv")
    log(f"  Saved: GSE118305_fpkm_genes.csv")
    return df_genes, meta


# ─── Parse GSE94892 (DENV PBMCs, pre-computed Cufflinks DEG) ─────────────────
def parse_gse94892():
    """Parse pre-computed Cufflinks DEG results from GSE94892."""
    log("Parsing GSE94892 (DENV PBMCs — Cufflinks DEG table) ...")

    # Group A = primary comparison
    df = pd.read_csv(RAW_DIR / "GSE94892" / "GSE94892_grpA.txt.gz", sep="\t")
    log(f"  Shape: {df.shape}")
    log(f"  Columns: {list(df.columns)}")
    log(f"  Sample 1 values: {df['sample_1'].unique()[:5]}")
    log(f"  Sample 2 values: {df['sample_2'].unique()[:5]}")

    # Rename columns to standard format
    df = df.rename(columns={
        "gene": "symbol",
        "log2(fold_change)": "log2FoldChange",
        "p_value": "pvalue",
        "q_value": "padj",
        "significant": "sig_cufflinks",
    })

    # Filter to valid status (not NOTEST)
    df_valid = df[df["status"] != "NOTEST"].copy()
    df_valid = df_valid.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FoldChange"])

    # DEGs: significant by Cufflinks
    degs = df_valid[
        (df_valid["sig_cufflinks"] == "yes") &
        (df_valid["padj"] < P_THRESH) &
        (df_valid["log2FoldChange"].abs() >= FC_THRESH)
    ]
    log(f"  Valid entries: {len(df_valid)}")
    log(f"  Cufflinks-significant DEGs: {len(degs)}")
    log(f"  DEGs up (Dengue > Control): {(degs['log2FoldChange'] > 0).sum()}")
    log(f"  DEGs down: {(degs['log2FoldChange'] < 0).sum()}")

    df_valid.to_csv(PROC_DIR / "DEGs_DENV_PBMCs_GSE94892.csv", index=False)
    log(f"  Saved: DEGs_DENV_PBMCs_GSE94892.csv")
    return df_valid, degs


# ─── Parse GSE78711 (ZIKV NPCs, pre-computed fold change) ────────────────────
def parse_gse78711():
    """Parse pre-computed fold change table from GSE78711."""
    log("Parsing GSE78711 (ZIKV NPCs — pre-computed FC table) ...")

    df = pd.read_csv(RAW_DIR / "GSE78711" / "GSE78711_expression.txt.gz", sep="\t")
    log(f"  Shape: {df.shape}")
    log(f"  Columns: {list(df.columns)}")
    log(f"  Head:")
    log(f"  {df.head(3).to_string()}")

    # Rename to standard
    df = df.rename(columns={
        "gene": "symbol",
        "log2.fold_change.": "log2FoldChange",
        "p_value": "pvalue",
        "significant": "sig_flag",
    })
    df["padj"] = df["pvalue"]  # no FDR in this dataset

    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FoldChange"])

    degs = df[
        (df["sig_flag"] == "yes") &
        (df["log2FoldChange"].abs() >= FC_THRESH)
    ]
    log(f"  Valid entries: {len(df)}")
    log(f"  Significant DEGs: {len(degs)}")
    log(f"  DEGs up (ZIKV > Mock): {(degs['log2FoldChange'] > 0).sum()}")
    log(f"  DEGs down: {(degs['log2FoldChange'] < 0).sum()}")

    df.to_csv(PROC_DIR / "DEGs_ZIKV_NPCs_GSE78711.csv", index=False)
    log(f"  Saved: DEGs_ZIKV_NPCs_GSE78711.csv")
    return df, degs


# ─── GATE G6: Replication Fisher Test ────────────────────────────────────────
def gate_g6_replication(shared_up, val_degs_dict, val_bg_dict):
    """Test whether discovery shared DEGs replicate in validation datasets."""
    log("\n" + "=" * 55)
    log("GATE G6: Replication Testing")
    log("=" * 55)

    results = []
    for dataset, val_degs in val_degs_dict.items():
        bg_genes = val_bg_dict[dataset]
        val_up = set(val_degs[val_degs["log2FoldChange"] > 0]["symbol"].dropna().tolist())

        # Fisher test
        shared = set(shared_up)
        a = len(shared & val_up)
        b = len(shared - val_up)
        c = len((bg_genes - shared) & val_up)
        d = len((bg_genes - shared) - val_up)

        if a > 0:
            or_val, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        else:
            or_val, p = 0.0, 1.0

        expected = len(shared) * len(val_up) / max(len(bg_genes), 1)
        fe = a / expected if expected > 0 else 0
        rep_rate = round(100 * a / len(shared), 1) if len(shared) > 0 else 0
        replicated_genes = sorted(shared & val_up)

        log(f"\n  [{dataset}]")
        log(f"    Validation upregulated DEGs: {len(val_up)}")
        log(f"    Discovery shared up: {len(shared)}")
        log(f"    Overlap: {a}  Expected: {expected:.2f}  FE: {fe:.2f}×")
        log(f"    Replication rate: {rep_rate}%")
        log(f"    Fisher p = {p:.4e}")
        log(f"    Replicated genes: {replicated_genes}")

        results.append({
            "dataset": dataset,
            "val_up_degs": len(val_up),
            "discovery_shared": len(shared),
            "overlap": a,
            "expected": round(expected, 3),
            "fold_enrichment": round(fe, 2),
            "replication_pct": rep_rate,
            "fisher_p": p,
            "replicated_genes": "; ".join(replicated_genes),
        })

    return pd.DataFrame(results)


def main():
    ckpt = load_ckpt()

    log("=" * 60)
    log("Step 12: External Validation (Phase 8, GATE G6)")
    log("=" * 60)

    # ─── Load discovery shared DEGs ───────────────────────────────────────────
    shared_ann = pd.read_csv(BASE_DIR / "03_results" / "phase3_shared_degs" / "shared_DEGs_annotated.csv")
    shared_up = shared_ann[shared_ann["log2FC_DENV"] > 0]["symbol"].dropna().tolist()
    log(f"Discovery shared upregulated genes ({len(shared_up)}): {shared_up}")

    val_degs = {}
    val_bg   = {}

    # ─── GSE118305 ────────────────────────────────────────────────────────────
    if not ckpt.get("gse118305_parsed"):
        try:
            gse118305_expr, gse118305_meta = parse_gse118305()
            # For GSE118305, we don't have a clean DEG table yet
            # Use gene names from the expression matrix as background
            genes_118305 = set(gse118305_expr.index.str.upper().tolist())
            log(f"  GSE118305: {len(genes_118305)} genes in expression matrix")
            # We'll manually create a DEG-like table based on FPKM ratio
            # Identify ZIKV vs Mock columns from metadata
            meta = gse118305_meta.copy()
            log(f"  GSE118305 metadata columns: {list(meta.columns)}")
            if "viral_infection" in meta.columns:
                log(f"  Unique viral infections: {meta['viral_infection'].unique()}")
            # Save genes for background
            pd.DataFrame({"gene": sorted(genes_118305)}).to_csv(RES_DIR / "GSE118305_background_genes.csv", index=False)
            ckpt["gse118305_parsed"] = True
            ckpt["gse118305_n_genes"] = len(genes_118305)
            save_ckpt(ckpt)
        except Exception as e:
            log(f"  GSE118305 parse error: {e}")
            ckpt["gse118305_parsed"] = False
            save_ckpt(ckpt)

    # ─── GSE94892 ─────────────────────────────────────────────────────────────
    if not ckpt.get("gse94892_parsed"):
        try:
            gse94892_all, gse94892_degs = parse_gse94892()
            val_degs["GSE94892 (DENV PBMCs)"] = gse94892_degs
            val_bg["GSE94892 (DENV PBMCs)"] = set(gse94892_all["symbol"].dropna().str.upper().tolist())
            ckpt["gse94892_parsed"] = True
            ckpt["gse94892_n_degs"] = len(gse94892_degs)
            save_ckpt(ckpt)
        except Exception as e:
            log(f"  GSE94892 parse error: {e}")
            ckpt["gse94892_parsed"] = False
            save_ckpt(ckpt)
    else:
        log("GSE94892: loading cached ...")
        gse94892_all = pd.read_csv(PROC_DIR / "DEGs_DENV_PBMCs_GSE94892.csv")
        gse94892_degs = gse94892_all[
            (gse94892_all.get("sig_cufflinks", pd.Series(["no"] * len(gse94892_all))) == "yes") &
            (gse94892_all["log2FoldChange"].abs() >= FC_THRESH)
        ]
        val_degs["GSE94892 (DENV PBMCs)"] = gse94892_degs
        val_bg["GSE94892 (DENV PBMCs)"]   = set(gse94892_all["symbol"].dropna().str.upper())

    # ─── GSE78711 ─────────────────────────────────────────────────────────────
    if not ckpt.get("gse78711_parsed"):
        try:
            gse78711_all, gse78711_degs = parse_gse78711()
            val_degs["GSE78711 (ZIKV NPCs)"] = gse78711_degs
            val_bg["GSE78711 (ZIKV NPCs)"] = set(gse78711_all["symbol"].dropna().str.upper().tolist())
            ckpt["gse78711_parsed"] = True
            ckpt["gse78711_n_degs"] = len(gse78711_degs)
            save_ckpt(ckpt)
        except Exception as e:
            log(f"  GSE78711 parse error: {e}")
            ckpt["gse78711_parsed"] = False
            save_ckpt(ckpt)
    else:
        log("GSE78711: loading cached ...")
        gse78711_all = pd.read_csv(PROC_DIR / "DEGs_ZIKV_NPCs_GSE78711.csv")
        gse78711_degs = gse78711_all[
            (gse78711_all.get("sig_flag", pd.Series(["no"] * len(gse78711_all))) == "yes") &
            (gse78711_all["log2FoldChange"].abs() >= FC_THRESH)
        ]
        val_degs["GSE78711 (ZIKV NPCs)"] = gse78711_degs
        val_bg["GSE78711 (ZIKV NPCs)"]   = set(gse78711_all["symbol"].dropna().str.upper())

    # ─── GATE G6 ──────────────────────────────────────────────────────────────
    if len(val_degs) > 0:
        gate_df = gate_g6_replication(
            [g.upper() for g in shared_up],
            {k: v.copy().assign(symbol=v["symbol"].str.upper()) for k, v in val_degs.items()},
            val_bg
        )
        gate_df.to_csv(RES_DIR / "gate_g6_replication_results.csv", index=False)
        log(f"\nGATE G6 results saved → {RES_DIR}/gate_g6_replication_results.csv")

        # ─── Figure ───────────────────────────────────────────────────────────
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Panel A — Replication rate per dataset
        ax = axes[0]
        labels = [d.split(" ")[0] for d in gate_df["dataset"].tolist()]
        rates  = gate_df["replication_pct"].tolist()
        colors = ["#D32F2F" if r > 10 else "#607D8B" for r in rates]
        bars = ax.bar(labels, rates, color=colors, edgecolor="black", alpha=0.85, width=0.5)
        for bar, r in zip(bars, rates):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                    f"{r:.0f}%", ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.axhline(10, color="orange", linestyle="--", linewidth=1.2, label="10% threshold")
        ax.set_ylabel("Replication Rate (%)", fontsize=11)
        ax.set_title("A  GATE G6: Discovery Gene\nReplication Rate", fontsize=11, fontweight="bold")
        ax.legend(fontsize=9)
        ax.set_ylim(0, max(rates + [15]) * 1.3)

        # Panel B — Fisher p-values
        ax2 = axes[1]
        ps = gate_df["fisher_p"].tolist()
        neg_log_p = [-np.log10(p) if p > 0 else 10 for p in ps]
        bar_colors2 = ["#D32F2F" if p < 0.05 else "#607D8B" for p in ps]
        bars2 = ax2.bar(labels, neg_log_p, color=bar_colors2, edgecolor="black", alpha=0.85, width=0.5)
        ax2.axhline(-np.log10(0.05), color="green", linestyle="--", linewidth=1.2, label="p=0.05")
        for bar, p in zip(bars2, ps):
            ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
                     f"p={p:.3f}", ha="center", va="bottom", fontsize=9)
        ax2.set_ylabel("-log10(Fisher p)", fontsize=11)
        ax2.set_title("B  GATE G6: Replication\nFisher Test", fontsize=11, fontweight="bold")
        ax2.legend(fontsize=9)

        plt.suptitle("GATE G6 — External Validation Replication Analysis", fontsize=13, fontweight="bold")
        plt.tight_layout()
        fig_path = FIG_MAIN / "Figure_GATE_G6_Validation.png"
        plt.savefig(fig_path, dpi=200, bbox_inches="tight")
        plt.close()
        log(f"Figure saved → {fig_path}")

        # ─── Gate decision ────────────────────────────────────────────────────
        log("\n" + "=" * 60)
        log("GATE G6 DECISION")
        log("=" * 60)

        any_sig = any(gate_df["fisher_p"] < 0.05)
        any_mod = any(gate_df["replication_pct"] > 10)
        n_shared_pathways = 0  # (would need Phase 7 pathway overlap)

        if any_sig:
            log("  GATE G6: STRONG REPLICATION — Fisher p < 0.05 in at least one dataset")
            gate_status = "STRONG"
        elif any_mod:
            log("  GATE G6: MODERATE REPLICATION — >10% replication rate")
            gate_status = "MODERATE"
        else:
            log("  GATE G6: WEAK/NO REPLICATION — Gene-level replication < 10%")
            log("  → Per SOP: check pathway-level replication instead")
            gate_status = "WEAK"

        ckpt["gate_g6_done"] = True
        ckpt["gate_g6_status"] = gate_status
        save_ckpt(ckpt)

    log("\nNext: run step13_network_analysis.py")


if __name__ == "__main__":
    main()
