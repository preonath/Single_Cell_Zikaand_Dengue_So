"""
Step 02: Load all single-cell TSV files → build AnnData → QC filter → Normalize → UMAP
Dataset: GSE110496 — Zanini et al. eLife 2018
Checkpoint-based: safe to restart at any stage.
"""

import os, json, time, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import sparse

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

# ─── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR  = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
RAW_DIR   = BASE_DIR / "00_raw_data" / "GSE110496"
ANN_DIR   = BASE_DIR / "01_processed_data" / "anndata_objects"
RES_DIR   = BASE_DIR / "03_results" / "phase2_qc"
FIG_DIR   = BASE_DIR / "04_figures" / "supplementary"
CKPT_DIR  = BASE_DIR / "checkpoints"
LOG_FILE  = BASE_DIR / "logs" / "step02_qc.log"

for d in [ANN_DIR, RES_DIR, FIG_DIR, CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE = CKPT_DIR / "step02_checkpoint.json"

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


# ─── Step A: Parse metadata ───────────────────────────────────────────────────
def parse_metadata() -> pd.DataFrame:
    meta_raw = pd.read_csv(RAW_DIR / "sample_metadata.csv")

    def parse_chars(s):
        d = {}
        for part in str(s).split(";"):
            part = part.strip()
            if ":" in part:
                k, v = part.split(":", 1)
                d[k.strip()] = v.strip()
        return d

    parsed = meta_raw["characteristics_ch1"].apply(parse_chars)
    df = pd.DataFrame(list(parsed))
    df["gsm"]        = meta_raw["geo_accession"].values
    df["cell_title"] = meta_raw["title"].values

    # Numeric conversions
    df["moi"]               = pd.to_numeric(df["moi"], errors="coerce").fillna(0).astype(int)
    df["time_h"]            = pd.to_numeric(df["time[h]"], errors="coerce").fillna(0).astype(int)
    df["n_dengue_mol"]      = pd.to_numeric(df.get("n_dengue_molecules", 0), errors="coerce").fillna(0).astype(int)
    df["n_zika_mol"]        = pd.to_numeric(df.get("n_zika_molecules", 0), errors="coerce").fillna(0).astype(int)
    df["viral_molecules"]   = df["n_dengue_mol"] + df["n_zika_mol"]

    # Assign condition
    def assign_condition(row):
        if row["moi"] == 0:
            return "Control"
        elif row["virus"] == "dengue":
            return "DENV"
        elif row["virus"] == "zika":
            return "ZIKV"
        return "Unknown"

    df["condition"] = df.apply(assign_condition, axis=1)
    df["timepoint"] = df["time_h"].astype(str) + "h"

    # Build filename map: gsm → tsv path
    tsv_files = {f.name.split("_")[0]: f for f in RAW_DIR.glob("GSM*_counts.tsv")}
    df["tsv_path"] = df["gsm"].map(tsv_files)
    df = df[df["tsv_path"].notna()].copy()

    log(f"Metadata parsed: {len(df)} cells with matching TSV files")
    log(f"  Condition counts: {df['condition'].value_counts().to_dict()}")
    log(f"  Timepoint counts: {df['timepoint'].value_counts().to_dict()}")
    return df


# ─── Step B: Load count matrix ────────────────────────────────────────────────
def load_count_matrix(meta: pd.DataFrame) -> sc.AnnData:
    log(f"Loading {len(meta)} TSV files ...")
    all_counts = {}
    gene_ids   = None

    for i, row in meta.iterrows():
        tsv = pd.read_csv(row["tsv_path"], sep="\t", index_col=0)
        if gene_ids is None:
            gene_ids = tsv.index.tolist()
        all_counts[row["gsm"]] = tsv["count"].values

        if (list(meta.index).index(i) + 1) % 200 == 0:
            log(f"  Loaded {list(meta.index).index(i)+1}/{len(meta)} cells ...")

    log("Building count matrix ...")
    count_matrix = np.array([all_counts[gsm] for gsm in meta["gsm"]], dtype=np.float32)

    obs_df = meta.set_index("gsm").copy()
    obs_df["tsv_path"] = obs_df["tsv_path"].astype(str)  # Path → str for h5ad

    adata = sc.AnnData(
        X    = sparse.csr_matrix(count_matrix),
        obs  = obs_df,
        var  = pd.DataFrame(index=gene_ids)
    )
    adata.var.index.name = "gene_id"
    adata.var_names_make_unique()

    log(f"AnnData created: {adata.shape[0]} cells × {adata.shape[1]} genes")
    return adata


# ─── Step C: QC metrics + filtering ──────────────────────────────────────────
def run_qc(adata: sc.AnnData) -> sc.AnnData:
    log("Calculating QC metrics ...")

    # Mitochondrial genes (Ensembl IDs starting with MT genes)
    # For Ensembl IDs: annotate MT genes separately
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    log(f"  Median genes per cell : {adata.obs['n_genes_by_counts'].median():.0f}")
    log(f"  Median counts per cell: {adata.obs['total_counts'].median():.0f}")
    log(f"  Median MT%            : {adata.obs['pct_counts_mt'].median():.1f}%")

    # ── Pre-filter violin plots ────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, col, title in zip(axes,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        ["Genes per cell", "Total counts", "MT %"]):
        for cond, grp in adata.obs.groupby("condition"):
            ax.violinplot(grp[col].dropna(), positions=[["Control","DENV","ZIKV"].index(cond)],
                         showmedians=True, widths=0.7)
        ax.set_xticks([0,1,2]); ax.set_xticklabels(["Control","DENV","ZIKV"])
        ax.set_title(title); ax.set_ylabel(col)
    fig.suptitle("QC metrics — BEFORE filtering", fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(FIG_DIR / "QC_violin_prefilter.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "QC_violin_prefilter.png", bbox_inches="tight", dpi=150)
    plt.close()
    log("  Pre-filter violin saved")

    # ── Determine thresholds (GSE110496: Smart-seq2-like) ─────────────────────
    MIN_GENES  = 2000
    MAX_GENES  = 8000
    MAX_MT_PCT = 15

    before = adata.n_obs
    adata = adata[
        (adata.obs["n_genes_by_counts"] > MIN_GENES) &
        (adata.obs["n_genes_by_counts"] < MAX_GENES) &
        (adata.obs["pct_counts_mt"]     < MAX_MT_PCT)
    ].copy()
    after = adata.n_obs

    removal_rate = 100 * (1 - after / before)
    log(f"QC filtering: {before} → {after} cells (removed {before-after}, {removal_rate:.1f}%)")
    log(f"  Cells per condition after filter:")
    for cond, n in adata.obs["condition"].value_counts().items():
        log(f"    {cond}: {n}")
    log(f"  Cells per timepoint after filter:")
    for tp, n in adata.obs["timepoint"].value_counts().sort_index().items():
        log(f"    {tp}: {n}")

    # Condition × timepoint table
    ct_table = pd.crosstab(adata.obs["condition"], adata.obs["timepoint"])
    ct_table.to_csv(RES_DIR / "cell_count_per_condition_timepoint.csv")
    log(f"  Cross-tab saved → cell_count_per_condition_timepoint.csv")
    log(f"\n{ct_table.to_string()}\n")

    # ── Post-filter violin ─────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, col, title in zip(axes,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        ["Genes per cell", "Total counts", "MT %"]):
        for cond, grp in adata.obs.groupby("condition"):
            ax.violinplot(grp[col].dropna(), positions=[["Control","DENV","ZIKV"].index(cond)],
                         showmedians=True, widths=0.7)
        ax.set_xticks([0,1,2]); ax.set_xticklabels(["Control","DENV","ZIKV"])
        ax.set_title(title); ax.set_ylabel(col)
    fig.suptitle("QC metrics — AFTER filtering", fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(FIG_DIR / "QC_violin_postfilter.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "QC_violin_postfilter.png", bbox_inches="tight", dpi=150)
    plt.close()
    log("  Post-filter violin saved")

    assert removal_rate < 50, f"Removal rate {removal_rate:.1f}% too high — check thresholds!"

    return adata


# ─── Step D: Normalize + HVG + PCA + UMAP ────────────────────────────────────
def normalize_and_embed(adata: sc.AnnData) -> sc.AnnData:
    log("Normalizing ...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log_norm"] = adata.X.copy()

    log("Finding highly variable genes (HVGs) ...")
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
    n_hvg = adata.var["highly_variable"].sum()
    log(f"  HVGs selected: {n_hvg}")

    # Check ISGs in HVGs
    isgs = ["ISG15", "MX1", "OAS1", "IFIT1", "IFITM3", "IFI6"]
    # Map gene names if possible (we have Ensembl IDs)
    log(f"  (Note: genes are Ensembl IDs — ISG check after annotation)")

    log("Scaling (regress out MT%) ...")
    sc.pp.regress_out(adata, ["pct_counts_mt"])
    sc.pp.scale(adata, max_value=10)

    log("Running PCA (50 components) ...")
    sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

    # Elbow plot
    fig, ax = plt.subplots(figsize=(8, 4))
    variance_ratio = adata.uns["pca"]["variance_ratio"]
    ax.plot(range(1, len(variance_ratio)+1), variance_ratio, "o-", markersize=3)
    ax.axvline(x=20, color="red", linestyle="--", label="n_pcs=20")
    ax.set_xlabel("PC"); ax.set_ylabel("Variance ratio")
    ax.set_title("PCA Elbow Plot — GSE110496")
    ax.legend(); ax.set_xlim(1, 50)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "elbow_plot.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "elbow_plot.png", bbox_inches="tight", dpi=150)
    plt.close()
    log("  Elbow plot saved")

    # Choose PCs based on elbow (typically 15-25 for this dataset)
    n_pcs = 20
    log(f"Using {n_pcs} PCs for UMAP ...")

    log("Computing neighbors + UMAP ...")
    sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=15)
    sc.tl.umap(adata)

    log("UMAP done.")
    return adata


# ─── Step E: UMAP visualizations ─────────────────────────────────────────────
def make_umap_figures(adata: sc.AnnData):
    log("Generating UMAP figures ...")

    COLORS = {"Control": "#808080", "DENV": "#E41A1C", "ZIKV": "#377EB8"}
    TP_COLORS = {"4h": "#FEE5D9", "12h": "#FCAE91", "24h": "#FB6A4A", "48h": "#CB181D"}

    # Figure 1: Condition UMAP
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    for cond in ["Control", "DENV", "ZIKV"]:
        mask = adata.obs["condition"] == cond
        ax.scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c=COLORS[cond], s=4, alpha=0.6, label=cond, rasterized=True
        )
    ax.set_title("UMAP — Condition", fontsize=13, fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.legend(markerscale=3, framealpha=0.8)
    ax.set_aspect("equal")

    ax = axes[1]
    for tp in ["4h", "12h", "24h", "48h"]:
        mask = adata.obs["timepoint"] == tp
        ax.scatter(
            adata.obsm["X_umap"][mask, 0],
            adata.obsm["X_umap"][mask, 1],
            c=TP_COLORS[tp], s=4, alpha=0.6, label=tp, rasterized=True
        )
    ax.set_title("UMAP — Timepoint", fontsize=13, fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.legend(markerscale=3, framealpha=0.8)
    ax.set_aspect("equal")

    fig.suptitle("GSE110496 — Single-cell RNA-seq (Huh7 cells)", fontsize=14)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "UMAP_condition_timepoint.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "UMAP_condition_timepoint.png", bbox_inches="tight", dpi=150)
    plt.close()
    log("  Condition+Timepoint UMAP saved")

    # Figure 2: Viral load UMAP
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, col, title, cmap in zip(axes,
        ["n_dengue_mol", "n_zika_mol"],
        ["DENV Viral molecules per cell", "ZIKV Viral molecules per cell"],
        ["Reds", "Blues"]):

        vals = np.log1p(adata.obs[col].values.astype(float))
        sc_plot = ax.scatter(
            adata.obsm["X_umap"][:, 0],
            adata.obsm["X_umap"][:, 1],
            c=vals, cmap=cmap, s=3, alpha=0.7, rasterized=True
        )
        plt.colorbar(sc_plot, ax=ax, label="log1p(viral molecules)")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.set_aspect("equal")

    fig.suptitle("GSE110496 — Viral load per cell (log1p)", fontsize=14)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "UMAP_viral_load.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "UMAP_viral_load.png", bbox_inches="tight", dpi=150)
    plt.close()
    log("  Viral load UMAP saved")

    # Figure 3: Condition × Timepoint facet grid
    conditions = ["Control", "DENV", "ZIKV"]
    timepoints = ["4h", "12h", "24h", "48h"]
    fig, axes = plt.subplots(len(conditions), len(timepoints),
                             figsize=(16, 10), sharex=True, sharey=True)
    umap_bg = adata.obsm["X_umap"]

    for i, cond in enumerate(conditions):
        for j, tp in enumerate(timepoints):
            ax = axes[i][j]
            # Background (all cells, grey)
            ax.scatter(umap_bg[:, 0], umap_bg[:, 1],
                      c="lightgrey", s=1, alpha=0.3, rasterized=True)
            # Highlight condition+timepoint
            mask = (adata.obs["condition"] == cond) & (adata.obs["timepoint"] == tp)
            n = mask.sum()
            ax.scatter(umap_bg[mask, 0], umap_bg[mask, 1],
                      c=COLORS[cond], s=4, alpha=0.8, rasterized=True)
            ax.set_title(f"{cond} {tp}\n(n={n})", fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle("GSE110496 — Condition × Timepoint UMAP facets", fontsize=13)
    plt.tight_layout()
    fig.savefig(FIG_DIR / "UMAP_facet_condition_timepoint.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "UMAP_facet_condition_timepoint.png", bbox_inches="tight", dpi=150)
    plt.close()
    log("  Facet UMAP saved")


# ─── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    ckpt = load_ckpt()
    log("=" * 60)
    log("Step 02: Load → QC → Normalize → UMAP")
    log("=" * 60)

    # ── A: Load raw AnnData ────────────────────────────────────────────────────
    raw_path = ANN_DIR / "adata_raw.h5ad"
    if ckpt.get("raw_loaded") and raw_path.exists():
        log(f"✓ Raw AnnData exists — loading from {raw_path.name}")
        adata = sc.read_h5ad(raw_path)
        log(f"  Shape: {adata.shape}")
    else:
        meta = parse_metadata()
        adata = load_count_matrix(meta)
        adata.write_h5ad(raw_path)
        log(f"✓ Raw AnnData saved → {raw_path.name}")
        ckpt["raw_loaded"] = True
        save_ckpt(ckpt)

    # ── B: QC + filtering ──────────────────────────────────────────────────────
    filtered_path = ANN_DIR / "adata_filtered.h5ad"
    if ckpt.get("qc_done") and filtered_path.exists():
        log(f"✓ Filtered AnnData exists — loading from {filtered_path.name}")
        adata = sc.read_h5ad(filtered_path)
        log(f"  Shape: {adata.shape}")
    else:
        adata = run_qc(adata)
        adata.write_h5ad(filtered_path)
        log(f"✓ Filtered AnnData saved → {filtered_path.name}")
        ckpt["qc_done"] = True
        ckpt["n_cells_after_qc"] = int(adata.n_obs)
        save_ckpt(ckpt)

    # ── C: Normalize + embed ───────────────────────────────────────────────────
    processed_path = ANN_DIR / "adata_processed.h5ad"
    if ckpt.get("embedding_done") and processed_path.exists():
        log(f"✓ Processed AnnData exists — loading from {processed_path.name}")
        adata = sc.read_h5ad(processed_path)
        log(f"  Shape: {adata.shape}")
    else:
        adata = normalize_and_embed(adata)
        adata.write_h5ad(processed_path)
        log(f"✓ Processed AnnData saved → {processed_path.name}")
        ckpt["embedding_done"] = True
        save_ckpt(ckpt)

    # ── D: UMAP figures ────────────────────────────────────────────────────────
    if ckpt.get("figures_done"):
        log("✓ Figures already generated — skipping")
    else:
        make_umap_figures(adata)
        ckpt["figures_done"] = True
        save_ckpt(ckpt)

    # ── Summary ────────────────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("STEP 02 COMPLETE — Summary:")
    log(f"  Final cells     : {adata.n_obs}")
    log(f"  Final genes     : {adata.n_vars}")
    log(f"  Conditions      : {adata.obs['condition'].value_counts().to_dict()}")
    log(f"  Figures saved in: {FIG_DIR}")
    log(f"  AnnData files   : {ANN_DIR}")
    log("\nNext: run step03_pseudobulk_deg.py")
    log("=" * 60)


if __name__ == "__main__":
    main()
