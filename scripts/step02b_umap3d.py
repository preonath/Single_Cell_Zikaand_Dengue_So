"""
Step 02b: 3D UMAP (interactive HTML) for GSE110496 Huh7 scRNA-seq.

Re-uses the neighbors graph already stored in adata_processed.h5ad (built in step02)
and recomputes UMAP with n_components=3 into a SEPARATE embedding key
(obsm['X_umap_3d']) so the canonical 2D X_umap is never overwritten.

Outputs (04_figures/supplementary/):
  - UMAP_3D_condition.html      interactive, coloured by condition
  - UMAP_3D_timepoint.html      interactive, coloured by timepoint
  - UMAP_3D_viral_load.html     interactive, coloured by viral molecule count
  - umap_3d_coords.csv          the 3 UMAP coords + obs labels

Checkpoint: checkpoints/step02b_checkpoint.json  (early-return if done — delete to recompute).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import plotly.express as px

BASE_DIR = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
ANN_DIR  = BASE_DIR / "01_processed_data" / "anndata_objects"
FIG_DIR  = BASE_DIR / "04_figures" / "supplementary"
CKPT_DIR = BASE_DIR / "checkpoints"
CKPT_FILE = CKPT_DIR / "step02b_checkpoint.json"
for d in [FIG_DIR, CKPT_DIR]:
    d.mkdir(parents=True, exist_ok=True)


def load_ckpt():
    return json.load(open(CKPT_FILE)) if CKPT_FILE.exists() else {}


def save_ckpt(d):
    json.dump(d, open(CKPT_FILE, "w"), indent=2)


def log(m):
    print(m, flush=True)


COND_ORDER = ["Control", "DENV", "ZIKV"]
TP_ORDER   = ["4h", "12h", "24h", "48h"]
COND_COLORS = {"Control": "#7F7F7F", "DENV": "#D62728", "ZIKV": "#1F77B4"}


def compute_3d(adata):
    """Recompute UMAP in 3D on the existing neighbors graph -> obsm['X_umap_3d']."""
    if "neighbors" not in adata.uns:
        log("  neighbors graph missing -> rebuilding (PCA + neighbors) ...")
        sc.pp.neighbors(adata, n_pcs=min(30, adata.obsm["X_pca"].shape[1]),
                        n_neighbors=15)
    log("  Computing 3D UMAP (n_components=3) ...")
    sc.tl.umap(adata, n_components=3)
    adata.obsm["X_umap_3d"] = adata.obsm["X_umap"].copy()
    return adata


def build_df(adata):
    c = adata.obsm["X_umap_3d"]
    df = pd.DataFrame(c, columns=["UMAP1", "UMAP2", "UMAP3"], index=adata.obs_names)
    df["condition"] = pd.Categorical(adata.obs["condition"].values,
                                     categories=COND_ORDER, ordered=True)
    df["timepoint"] = pd.Categorical(adata.obs["timepoint"].astype(str).values,
                                     categories=TP_ORDER, ordered=True)
    df["viral_molecules"] = adata.obs["viral_molecules"].values
    df["log_viral"] = np.log10(df["viral_molecules"].clip(lower=0) + 1)
    return df


def save_html(df, color, fname, title, **kw):
    fig = px.scatter_3d(df, x="UMAP1", y="UMAP2", z="UMAP3", color=color,
                        title=title, opacity=0.8, **kw)
    fig.update_traces(marker=dict(size=3))
    fig.update_layout(legend=dict(itemsizing="constant"),
                      scene=dict(xaxis_title="UMAP 1", yaxis_title="UMAP 2",
                                 zaxis_title="UMAP 3"),
                      template="plotly_white")
    out = FIG_DIR / fname
    fig.write_html(out)
    log(f"  wrote {out.name}")


def main():
    ckpt = load_ckpt()
    if ckpt.get("html_done"):
        log("step02b already done (checkpoint present) — delete "
            "checkpoints/step02b_checkpoint.json to recompute.")
        return

    log("Loading adata_processed.h5ad ...")
    adata = sc.read_h5ad(ANN_DIR / "adata_processed.h5ad")
    log(f"  cells={adata.n_obs}, genes={adata.n_vars}")

    adata = compute_3d(adata)
    df = build_df(adata)
    df.to_csv(FIG_DIR / "umap_3d_coords.csv")
    log(f"  wrote umap_3d_coords.csv ({len(df)} cells)")

    base = "GSE110496 — Huh7 scRNA-seq — 3D UMAP"
    save_html(df, "condition", "UMAP_3D_condition.html",
              f"{base} (Condition)", color_discrete_map=COND_COLORS)
    save_html(df, "timepoint", "UMAP_3D_timepoint.html",
              f"{base} (Timepoint)",
              color_discrete_sequence=px.colors.sequential.YlOrRd)
    save_html(df, "log_viral", "UMAP_3D_viral_load.html",
              f"{base} (log10 viral molecules)",
              color_continuous_scale="Viridis")

    save_ckpt({"embedding_3d_done": True, "html_done": True,
               "n_cells": int(adata.n_obs)})
    log("step02b done.")


if __name__ == "__main__":
    main()
