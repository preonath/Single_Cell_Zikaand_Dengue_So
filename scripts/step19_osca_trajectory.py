"""
Step 19: OSCA-style trajectory analysis (scanpy equivalent of the Bioconductor OSCA
'trajectory-analysis' chapter — TSCAN/MST + pseudotime + genes-along-pseudotime).

OSCA chapter            ->  scanpy equivalent used here
  cluster cells         ->  sc.tl.leiden
  MST on cluster centroids (TSCAN)  ->  sc.tl.paga  (partition-based graph abstraction =
                                        abstracted connectivity backbone between clusters)
  root + map pseudotime ->  sc.tl.dpt rooted at the lowest viral-load cell (diffusion pseudotime)
  project / draw path   ->  PAGA edges drawn over the UMAP embedding
  genes vs pseudotime   ->  expression of the 15 shared genes ordered by DPT pseudotime

Infection anchoring: pseudotime is rooted at the lowest-`viral_molecules` cell and we report
pseudotime<->viral-load correlation, so the unsupervised trajectory is tied to infection progress.

Canonical adata_processed.h5ad is READ-ONLY; all clustering/PAGA/DPT run on a copy.
Outputs -> 04_figures/supplementary/ ; checkpoint checkpoints/step19_checkpoint.json.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
ANN  = BASE / "01_processed_data" / "anndata_objects" / "adata_processed.h5ad"
DEG  = BASE / "01_processed_data" / "deg_tables" / "DEGs_DENV_vs_Control_annotated_full.csv"
SHARED = BASE / "03_results" / "phase3_shared_degs" / "shared_DEGs_annotated.csv"
FIG  = BASE / "04_figures" / "supplementary"
RES  = BASE / "03_results" / "phase_trajectory"
CKPT = BASE / "checkpoints" / "step19_checkpoint.json"
for d in [FIG, RES, CKPT.parent]:
    d.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = FIG
sc.settings.verbosity = 1


def log(m): print(m, flush=True)


def main():
    if CKPT.exists() and json.load(open(CKPT)).get("done"):
        log("step19 already done — delete checkpoints/step19_checkpoint.json to recompute.")
        return

    log("Loading adata_processed.h5ad (read-only) -> working on a copy ...")
    a = sc.read_h5ad(ANN).copy()
    a.obs["log_viral"] = np.log10(a.obs["viral_molecules"].clip(lower=0) + 1)

    # --- 1. cluster (OSCA: cluster cells) ---
    log("Leiden clustering ...")
    sc.tl.leiden(a, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2,
                 directed=False)
    log(f"  {a.obs['leiden'].nunique()} clusters")

    # --- 2. PAGA: abstracted graph backbone (OSCA: MST on centroids / TSCAN) ---
    log("PAGA abstracted graph ...")
    sc.tl.paga(a, groups="leiden")
    sc.pl.paga(a, plot=False)   # compute layout positions (uns['paga']['pos']) without drawing

    # --- 3. embedding initialised from PAGA, so the backbone is drawable over it ---
    sc.tl.umap(a, init_pos="paga")

    # --- 4. root at lowest viral-load cell, diffusion pseudotime (OSCA: root + map pseudotime) ---
    log("Diffusion pseudotime rooted at lowest-viral-load cell ...")
    sc.tl.diffmap(a)
    a.uns["iroot"] = int(np.argmin(a.obs["viral_molecules"].values))
    sc.tl.dpt(a)
    rho = spearmanr(a.obs["viral_molecules"], a.obs["dpt_pseudotime"])[0]
    log(f"  pseudotime <-> viral load Spearman rho = {rho:.3f}")

    # --- 4b. ANNOTATE clusters (no cell types in a clonal line -> annotate by INFECTION STATE) ---
    log("Annotating clusters by infection state + dominant condition ...")
    g = a.obs.groupby("leiden", observed=True)
    cl = g.agg(n=("leiden", "size"),
               median_load=("viral_molecules", "median"),
               median_log_viral=("log_viral", "median"),
               median_dpt=("dpt_pseudotime", "median"),
               pct_infected=("viral_molecules", lambda s: float((s > 0).mean()))).copy()
    dom = g["condition"].agg(lambda s: s.value_counts().index[0])
    cl["dominant_condition"] = dom
    # infection-state label from median viral load (anchored, interpretable)
    def state(row):
        if row["pct_infected"] < 0.3 or row["median_load"] < 1:
            return "Uninfected / bystander"
        # terciles of dpt among infected-dominated clusters
        return None
    inf_clusters = cl[(cl["pct_infected"] >= 0.3) & (cl["median_load"] >= 1)].index
    dpt_inf = cl.loc[inf_clusters, "median_dpt"]
    q1, q2 = dpt_inf.quantile([0.34, 0.67]) if len(dpt_inf) >= 3 else (dpt_inf.min(), dpt_inf.max())
    labels = {}
    for c in cl.index:
        if c not in inf_clusters:
            labels[c] = "Uninfected / bystander"
        elif cl.loc[c, "median_dpt"] <= q1:
            labels[c] = "Early infection"
        elif cl.loc[c, "median_dpt"] <= q2:
            labels[c] = "Mid infection"
        else:
            labels[c] = "Late / high-load"
    cl["cell_state"] = pd.Series(labels)
    cl["cluster_label"] = [f"{cl.loc[c,'dominant_condition']}: {cl.loc[c,'cell_state']}" for c in cl.index]
    a.obs["cell_state"] = pd.Categorical(a.obs["leiden"].map(labels),
        categories=["Uninfected / bystander", "Early infection", "Mid infection", "Late / high-load"],
        ordered=True)
    a.obs["cluster_label"] = a.obs["leiden"].map(cl["cluster_label"].to_dict())
    cl.to_csv(RES / "osca_cluster_annotation.csv")
    log("  cluster annotation:")
    log(cl[["n", "dominant_condition", "pct_infected", "median_load", "median_dpt", "cell_state"]].to_string())

    # ---------- Figure 1: PAGA backbone ----------
    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.paga(a, color="leiden", ax=ax, show=False, fontsize=9,
               title="PAGA backbone (cluster connectivity = OSCA MST analog)")
    fig.savefig(FIG / "osca_paga_graph.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---------- Figure 2: UMAP panels (annotated state / condition / viral load / pseudotime) ----------
    state_pal = {"Uninfected / bystander": "#BDBDBD", "Early infection": "#FEE08B",
                 "Mid infection": "#FC8D59", "Late / high-load": "#D73027"}
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    sc.pl.umap(a, color="cell_state", ax=axes[0, 0], show=False,
               title="Annotated infection state (cluster annotation)", palette=state_pal)
    sc.pl.umap(a, color="condition", ax=axes[0, 1], show=False, title="Condition",
               palette={"Control": "#7F7F7F", "DENV": "#D62728", "ZIKV": "#1F77B4"})
    sc.pl.umap(a, color="log_viral", ax=axes[1, 0], show=False, title="log10 viral load", color_map="viridis")
    sc.pl.umap(a, color="dpt_pseudotime", ax=axes[1, 1], show=False,
               title=f"DPT pseudotime (rho vs load = {rho:.2f})", color_map="magma")
    fig.tight_layout()
    fig.savefig(FIG / "osca_umap_pseudotime_panels.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---------- Figure 3: trajectory backbone drawn over the UMAP ----------
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(a, color="dpt_pseudotime", ax=ax, show=False, color_map="magma",
               title="Trajectory (PAGA edges) over UMAP, coloured by pseudotime", size=25)
    sc.pl.paga(a, color="leiden", ax=ax, show=False, node_size_scale=1.2,
               edge_width_scale=0.7, fontsize=8, frameon=False, text_kwds={"alpha": 0})
    fig.savefig(FIG / "osca_trajectory_over_umap.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---------- Figure 4: 15 shared genes along pseudotime ----------
    deg = pd.read_csv(DEG); ens2sym = dict(zip(deg["gene_id"], deg["symbol"].astype(str)))
    shared = pd.read_csv(SHARED)
    ids = [g for g in shared["gene_id"] if g in set(a.var_names)]
    syms = [ens2sym.get(g, g) for g in ids]
    X = a.layers["log_norm"]; X = X.toarray() if sp.issparse(X) else np.asarray(X)
    vn = pd.Index(a.var_names)
    order = np.argsort(a.obs["dpt_pseudotime"].values)
    M = np.vstack([X[order, vn.get_loc(g)] for g in ids])
    Mz = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9)

    fig, (axh, axb) = plt.subplots(2, 1, figsize=(10, 6), height_ratios=[6, 1.2], sharex=True)
    im = axh.imshow(Mz, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2, interpolation="nearest")
    axh.set_yticks(range(len(syms))); axh.set_yticklabels(syms, fontsize=8)
    axh.set_title("15 shared genes ordered by DPT pseudotime (left=root/uninfected -> right=late/high-load)",
                  fontweight="bold", fontsize=11)
    fig.colorbar(im, ax=axh, fraction=0.012, pad=0.01, label="z-score")
    pt = a.obs["dpt_pseudotime"].values[order]; lv = a.obs["log_viral"].values[order]
    axb.plot(pt, c="k", lw=1, label="pseudotime")
    axb2 = axb.twinx(); axb2.plot(lv, c="#2ca02c", lw=1, alpha=0.7, label="log10 load")
    axb.set_ylabel("pseudotime", fontsize=8); axb2.set_ylabel("log10 load", fontsize=8, color="#2ca02c")
    axb.set_xlabel("cells ordered by pseudotime"); axb.spines[["top"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "osca_genes_along_pseudotime.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    json.dump({"done": True, "n_clusters": int(a.obs["leiden"].nunique()),
               "pseudotime_vs_viralload_rho": round(float(rho), 3)}, open(CKPT, "w"), indent=2)
    log("\nFigures written:")
    for f in ["osca_paga_graph", "osca_umap_pseudotime_panels", "osca_trajectory_over_umap",
              "osca_genes_along_pseudotime"]:
        log(f"  {f}.png")
    log("step19 done.")


if __name__ == "__main__":
    main()
