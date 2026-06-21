"""
Step 18: Viral-load infection-progression trajectory (single cell).

NOT a developmental pseudotime. Huh7 is a homogeneous line with no lineage tree, so we use a
SUPERVISED ordering of cells by measured intracellular viral RNA (`viral_molecules`) — the
"virus-inclusive scRNA-seq" approach of Zanini 2018 (the source study of this dataset). We then
ask how the 15 shared convergent genes change along that infection axis, in DENV and ZIKV.

Diagnostics that justify this (computed here, printed):
  - viral load vs PC1 (infection should be a dominant axis)
  - viral load vs sequencing depth   (confound check; expect ~0 / negative = host shutoff)
  - viral load vs cell-cycle S/G2M score (confound check)
  - DPT (rooted at lowest-load cell) vs viral load (does an unsupervised pseudotime recover load?)

Bystander handling: cells with viral_molecules==0 in an infected condition are exposed-but-uninfected
(paracrine state), NOT "early infection" — flagged and excluded from the trend fits (shown separately).

Outputs (04_figures/supplementary/ + 03_results/phase_trajectory/):
  - trajectory_convergence_genes.png/pdf      headline: genes rising in BOTH viruses vs load
  - trajectory_shared_genes_{DENV,ZIKV}.png   all 15 genes, lowess vs log10 viral load
  - trajectory_heatmap_{DENV,ZIKV}.png         cells ordered by load, 15-gene expression heatmap
  - viralload_trajectory_stats.csv             per gene per virus: rho, p, dpt_rho, cc/depth confound
Checkpoint: checkpoints/step18_checkpoint.json  (early-return if done — delete to recompute).
Canonical adata_processed.h5ad is read-only; DPT is computed on a copy.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
ANN  = BASE / "01_processed_data" / "anndata_objects" / "adata_processed.h5ad"
DEG  = BASE / "01_processed_data" / "deg_tables" / "DEGs_DENV_vs_Control_annotated_full.csv"
SHARED = BASE / "03_results" / "phase3_shared_degs" / "shared_DEGs_annotated.csv"
FIG  = BASE / "04_figures" / "supplementary"
RES  = BASE / "03_results" / "phase_trajectory"
CKPT = BASE / "checkpoints" / "step18_checkpoint.json"
for d in [FIG, RES, CKPT.parent]:
    d.mkdir(parents=True, exist_ok=True)

VIRUSES = ["DENV", "ZIKV"]
VCOL = {"DENV": "#D62728", "ZIKV": "#1F77B4"}
# the 4 genes that dose-respond in BOTH viruses (from the diagnostic) = convergence headline
CONVERGENT = ["CREBRF", "TSPYL2", "INHBE", "DUSP1"]


def log(m): print(m, flush=True)


def lognorm_matrix(adata):
    X = adata.layers["log_norm"]
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def cell_cycle_score(adata, sym2ens):
    """Score S/G2M using the standard regev list mapped to Ensembl; works on log_norm."""
    s_genes = ["MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7",
               "DTL","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76",
               "SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6",
               "EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8"]
    g2m_genes = ["HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2",
                 "CKS1B","MKI67","TMPO","CENPF","TACC3","SMC4","CCNB2","CKAP2L","CKAP2","AURKB",
                 "BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","CDC20",
                 "TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2",
                 "KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5"]
    tmp = adata.copy()
    tmp.X = tmp.layers["log_norm"].copy()
    tmp.var_names = [sym2ens.get(s, e) for s, e in zip(tmp.var.get("symbol", tmp.var_names), tmp.var_names)] \
        if "symbol" in tmp.var else tmp.var_names
    s_e   = [sym2ens[g] for g in s_genes   if g in sym2ens and sym2ens[g] in set(adata.var_names)]
    g2m_e = [sym2ens[g] for g in g2m_genes if g in sym2ens and sym2ens[g] in set(adata.var_names)]
    try:
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_e, g2m_genes=g2m_e, use_raw=False)
        # score_genes_cell_cycle uses adata.X (scaled) — fine for a relative confound check
        return True
    except Exception as e:
        log(f"  cell-cycle scoring skipped: {e}")
        adata.obs["S_score"] = np.nan
        adata.obs["G2M_score"] = np.nan
        return False


def main():
    if CKPT.exists() and json.load(open(CKPT)).get("done"):
        log("step18 already done — delete checkpoints/step18_checkpoint.json to recompute.")
        return

    log("Loading adata_processed.h5ad (read-only) ...")
    adata = sc.read_h5ad(ANN)

    # symbol<->ensembl from the full DEG table; map shared 15 genes to Ensembl
    deg = pd.read_csv(DEG)
    sym2ens = dict(zip(deg["symbol"].astype(str), deg["gene_id"]))
    ens2sym = dict(zip(deg["gene_id"], deg["symbol"].astype(str)))
    shared = pd.read_csv(SHARED)
    shared_ids = [g for g in shared["gene_id"] if g in set(adata.var_names)]
    shared_sym = [ens2sym.get(g, g) for g in shared_ids]
    log(f"  {len(shared_ids)}/15 shared genes present in object")

    # cell-cycle score (sensitivity)
    cell_cycle_score(adata, sym2ens)

    # DPT pseudotime on a COPY, rooted at lowest-load infected cell (per all cells)
    log("Computing DPT (rooted at lowest-load cell) on a copy ...")
    ad2 = adata.copy()
    sc.tl.diffmap(ad2)
    root = int(np.argmin(ad2.obs["viral_molecules"].values))
    ad2.uns["iroot"] = root
    sc.tl.dpt(ad2)
    adata.obs["dpt"] = ad2.obs["dpt_pseudotime"].values

    X = lognorm_matrix(adata)
    vn = pd.Index(adata.var_names)
    var_idx = {g: vn.get_loc(g) for g in shared_ids}

    stats_rows = []
    vload_all = adata.obs["viral_molecules"].values.astype(float)
    pcs = adata.obsm["X_pca"]

    for virus in VIRUSES:
        m = (adata.obs["condition"] == virus).values
        idx = np.where(m)[0]
        vl = vload_all[idx]
        infected = vl > 0                      # bystander handling
        depth = adata.obs["total_counts"].values[idx].astype(float)
        s_sc = adata.obs.get("S_score", pd.Series(np.nan, index=adata.obs_names)).values[idx]
        g2m  = adata.obs.get("G2M_score", pd.Series(np.nan, index=adata.obs_names)).values[idx]

        # axis diagnostics
        pc_r = [abs(spearmanr(vl, pcs[idx, k])[0]) for k in range(min(10, pcs.shape[1]))]
        top_pc = int(np.argmax(pc_r)) + 1
        rho_depth = spearmanr(vl, depth)[0]
        rho_cc = spearmanr(vl, np.nan_to_num(g2m))[0] if np.isfinite(g2m).any() else np.nan
        rho_dpt = spearmanr(vl, adata.obs["dpt"].values[idx])[0]
        log(f"\n[{virus}] n={m.sum()} infected={infected.sum()} ({100*infected.mean():.0f}%) | "
            f"viral↔PC{top_pc} rho={max(pc_r):.2f} | viral↔depth={rho_depth:.2f} | "
            f"viral↔G2M={rho_cc if rho_cc==rho_cc else float('nan'):.2f} | viral↔DPT={rho_dpt:.2f}")

        # per-gene Spearman vs load (infected cells only, to avoid bystander pull)
        for gid in shared_ids:
            expr = X[idx, var_idx[gid]]
            r_all, p_all = spearmanr(vl, expr)
            r_inf, p_inf = spearmanr(vl[infected], expr[infected]) if infected.sum() > 10 else (np.nan, np.nan)
            stats_rows.append(dict(virus=virus, gene=ens2sym.get(gid, gid),
                                   n=int(m.sum()), n_infected=int(infected.sum()),
                                   rho_load_all=round(r_all, 3), p_all=p_all,
                                   rho_load_infected=round(r_inf, 3), p_infected=p_inf,
                                   viral_vs_depth=round(rho_depth, 3),
                                   viral_vs_G2M=round(rho_cc, 3) if rho_cc == rho_cc else np.nan,
                                   viral_vs_dpt=round(rho_dpt, 3), top_pc=top_pc,
                                   top_pc_rho=round(max(pc_r), 2)))

    stats = pd.DataFrame(stats_rows)
    stats.to_csv(RES / "viralload_trajectory_stats.csv", index=False)
    log(f"\nWrote {RES/'viralload_trajectory_stats.csv'}")

    # ---------- Figures ----------
    def expr_of(gid, idx):
        return X[idx, var_idx[gid]]

    # (1) headline convergence figure: 4 genes rising in BOTH viruses
    conv_ids = [sym2ens[g] for g in CONVERGENT if g in sym2ens and sym2ens[g] in var_idx]
    fig, axes = plt.subplots(1, len(conv_ids), figsize=(4.2 * len(conv_ids), 4), squeeze=False)
    for ax, gid in zip(axes[0], conv_ids):
        sym = ens2sym.get(gid, gid)
        for virus in VIRUSES:
            idx = np.where((adata.obs["condition"] == virus).values)[0]
            vl = vload_all[idx]; inf = vl > 0
            x = np.log10(vl[inf] + 1); y = expr_of(gid, idx)[inf]
            ax.scatter(x, y, s=6, c=VCOL[virus], alpha=0.25)
            if inf.sum() > 15:
                lo = sm.nonparametric.lowess(y, x, frac=0.7)
                ax.plot(lo[:, 0], lo[:, 1], c=VCOL[virus], lw=2.5, label=virus)
        ax.set_title(sym, fontweight="bold"); ax.set_xlabel("log10(viral molecules+1)")
        ax.set_ylabel("expression (log-norm)"); ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Convergent shared genes rise with viral load in BOTH viruses\n"
                 "(infection-progression trajectory; infected cells only)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "trajectory_convergence_genes.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG / "trajectory_convergence_genes.pdf", bbox_inches="tight")
    plt.close(fig)
    log("  wrote trajectory_convergence_genes.png/.pdf")

    # (2) all 15 genes vs load, per virus (lowess panel) + (3) heatmap ordered by load
    for virus in VIRUSES:
        idx = np.where((adata.obs["condition"] == virus).values)[0]
        vl = vload_all[idx]; inf = vl > 0
        order = idx[np.argsort(vl)]                  # cells ordered by infection progression
        # lowess panel
        fig, ax = plt.subplots(figsize=(7, 6))
        for gid in shared_ids:
            x = np.log10(vl[inf] + 1); y = expr_of(gid, idx)[inf]
            if inf.sum() > 15:
                lo = sm.nonparametric.lowess(y, x, frac=0.8)
                rho = stats[(stats.virus == virus) & (stats.gene == ens2sym.get(gid, gid))]["rho_load_infected"].values
                lab = f"{ens2sym.get(gid,gid)} (ρ={rho[0]:.2f})" if len(rho) else ens2sym.get(gid, gid)
                ax.plot(lo[:, 0], lo[:, 1], lw=1.8, label=lab)
        ax.set_title(f"{virus}: 15 shared genes along viral-load trajectory", fontweight="bold")
        ax.set_xlabel("log10(viral molecules+1)"); ax.set_ylabel("expression (log-norm)")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(fontsize=7, ncol=2, frameon=False)
        fig.tight_layout(); fig.savefig(FIG / f"trajectory_shared_genes_{virus}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

        # heatmap: rows=genes, cols=cells ordered by load, z-scored per gene
        M = np.vstack([X[order, var_idx[g]] for g in shared_ids])
        Mz = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9)
        fig, (axh, axb) = plt.subplots(2, 1, figsize=(9, 5.5), height_ratios=[6, 1], sharex=True)
        im = axh.imshow(Mz, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2,
                        interpolation="nearest")
        axh.set_yticks(range(len(shared_sym))); axh.set_yticklabels(shared_sym, fontsize=8)
        axh.set_title(f"{virus}: shared-gene expression, cells ordered by viral load "
                      f"(left=low → right=high)", fontweight="bold", fontsize=11)
        fig.colorbar(im, ax=axh, fraction=0.012, pad=0.01, label="z-score")
        axb.plot(np.log10(np.sort(vl) + 1), c="k", lw=1)
        axb.set_ylabel("log10\nload", fontsize=8); axb.set_xlabel("cells (ordered by viral load)")
        axb.spines[["top", "right"]].set_visible(False)
        fig.tight_layout(); fig.savefig(FIG / f"trajectory_heatmap_{virus}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        log(f"  wrote trajectory_shared_genes_{virus}.png + trajectory_heatmap_{virus}.png")

    json.dump({"done": True, "n_cells": int(adata.n_obs),
               "convergent_genes": CONVERGENT}, open(CKPT, "w"), indent=2)
    log("\nstep18 done.")


if __name__ == "__main__":
    main()
