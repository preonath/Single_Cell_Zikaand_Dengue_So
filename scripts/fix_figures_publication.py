"""
Fix all figure quality issues for publication.
Fixes applied:
  1. QC Violins  — remove flat MT% panel, larger fonts
  2. UMAP Facet  — readable panel titles and fonts
  3. Volcano plots (standalone) — adjustText to fix label overlaps
  4. FC Correlation — fix p-value display, larger figure
  5. Figure 2 — fix volcano label overlaps, larger FC scatter
  6. Figure 3 — larger fonts in panels B & C
  7. Figure 4 — fix KEGG pathway name truncation
  8. Figure 5 — fix network node label overlaps
"""

import warnings, time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pathlib import Path
from scipy.stats import pearsonr
from adjustText import adjust_text
import networkx as nx

warnings.filterwarnings("ignore")

BASE  = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
ANN   = BASE / "01_processed_data/anndata_objects"
DEG   = BASE / "01_processed_data/deg_tables"
RES   = BASE / "03_results"
FMAIN = BASE / "04_figures/main"
FSUPP = BASE / "04_figures/supplementary"

COLORS = {
    "denv":      "#D32F2F",
    "zikv":      "#1565C0",
    "shared":    "#FF6F00",
    "novel":     "#607D8B",
    "validated": "#E91E63",
    "mirna_hub": "#8E24AA",
    "down":      "#1565C0",
    "up":        "#D32F2F",
    "grey":      "#BDBDBD",
    "control":   "#546E7A",
}

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}")


# ══════════════════════════════════════════════════════════════════════
# 1. QC VIOLINS — remove MT% panel, 2-panel layout, larger fonts
# ══════════════════════════════════════════════════════════════════════
def fix_qc_violins():
    log("Fixing QC violin plots ...")
    adata_raw = sc.read_h5ad(ANN / "adata_raw.h5ad")
    adata_flt = sc.read_h5ad(ANN / "adata_processed.h5ad")

    # Recompute n_genes and total_counts from raw if not present
    for ad in [adata_raw, adata_flt]:
        if "n_genes_by_counts" not in ad.obs.columns:
            sc.pp.calculate_qc_metrics(ad, inplace=True)

    cond_order  = ["Control", "DENV", "ZIKV"]
    cond_colors = [COLORS["control"], COLORS["denv"], COLORS["zikv"]]
    metrics     = ["n_genes_by_counts", "total_counts"]
    labels      = ["Genes per cell", "Total UMI counts"]

    for tag, ad in [("prefilter", adata_raw), ("postfilter", adata_flt)]:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        title_sfx = "BEFORE filtering" if tag == "prefilter" else "AFTER filtering"
        fig.suptitle(f"QC metrics — {title_sfx}", fontsize=14, fontweight="bold", y=1.02)

        for ax, metric, label in zip(axes, metrics, labels):
            data = [ad.obs.loc[ad.obs["condition"] == c, metric].values for c in cond_order]
            parts = ax.violinplot(data, positions=range(len(cond_order)),
                                  showmedians=True, showextrema=True)
            for i, (pc, col) in enumerate(zip(parts["bodies"], cond_colors)):
                pc.set_facecolor(col)
                pc.set_alpha(0.7)
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(2)
            parts["cbars"].set_color("black")
            parts["cmins"].set_color("black")
            parts["cmaxes"].set_color("black")
            ax.set_xticks(range(len(cond_order)))
            ax.set_xticklabels(cond_order, fontsize=12)
            ax.set_ylabel(label, fontsize=12)
            ax.tick_params(axis="y", labelsize=11)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            n_total = len(ad)
            ax.set_title(f"{label}\n(n = {n_total} cells)", fontsize=12)

        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FSUPP / f"QC_violin_{tag}.{ext}", dpi=200, bbox_inches="tight")
        plt.close()
        log(f"  Saved QC_violin_{tag}")


# ══════════════════════════════════════════════════════════════════════
# 2. UMAP FACET — readable titles, larger dots, bigger fonts
# ══════════════════════════════════════════════════════════════════════
def fix_umap_facet():
    log("Fixing UMAP facet plot ...")
    adata = sc.read_h5ad(ANN / "adata_processed.h5ad")

    cond_order = ["Control", "DENV", "ZIKV"]
    tp_order   = ["4h", "12h", "24h", "48h"]
    cond_cols  = {"Control": COLORS["control"], "DENV": COLORS["denv"], "ZIKV": COLORS["zikv"]}

    fig, axes = plt.subplots(3, 4, figsize=(18, 13))
    fig.suptitle("GSE110496 — Condition × Timepoint UMAP Facets\n(Huh7 Single-Cell RNA-seq)",
                 fontsize=15, fontweight="bold", y=1.01)

    umap = adata.obsm["X_umap"]

    for row_i, cond in enumerate(cond_order):
        for col_i, tp in enumerate(tp_order):
            ax = axes[row_i, col_i]
            # Background: all other cells in light grey
            ax.scatter(umap[:, 0], umap[:, 1], c="#E0E0E0", s=3, alpha=0.3,
                       rasterized=True, linewidths=0)
            # Foreground: this condition × timepoint
            mask = (adata.obs["condition"] == cond) & (adata.obs["timepoint"] == tp)
            n = mask.sum()
            if n > 0:
                ax.scatter(umap[mask, 0], umap[mask, 1],
                           c=cond_cols[cond], s=10, alpha=0.8,
                           rasterized=True, linewidths=0)
            ax.set_title(f"{cond} — {tp}\n(n = {n})", fontsize=11, fontweight="bold")
            ax.set_xticks([]); ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            if col_i == 0:
                ax.set_ylabel(cond, fontsize=12, fontweight="bold", rotation=90, labelpad=5)
            if row_i == 0:
                ax.set_title(f"{tp}\n{cond} (n={n})", fontsize=11, fontweight="bold")

    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FSUPP / f"UMAP_facet_condition_timepoint.{ext}", dpi=200, bbox_inches="tight")
    plt.close()
    log("  Saved UMAP_facet_condition_timepoint")


# ══════════════════════════════════════════════════════════════════════
# 3. VOLCANO PLOTS (standalone) — adjustText, larger fonts
# ══════════════════════════════════════════════════════════════════════
def draw_volcano_fixed(ax, df, title, shared_syms, fc_thresh=1.0, p_thresh=0.05):
    df = df.copy().dropna(subset=["log2FoldChange", "padj"])
    df["_lfc"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df["_p"]   = pd.to_numeric(df["padj"], errors="coerce")
    df = df.dropna(subset=["_lfc", "_p"])
    df["_nlp"] = -np.log10(df["_p"].clip(lower=1e-300))
    # Use symbol column for gene name lookup (falls back to gene_id)
    label_col = "symbol" if "symbol" in df.columns else "gene_id"

    is_up = (df["_lfc"] >  fc_thresh) & (df["_p"] < p_thresh)
    is_dn = (df["_lfc"] < -fc_thresh) & (df["_p"] < p_thresh)
    is_sh = df[label_col].str.upper().isin(shared_syms)

    # Background NS
    ax.scatter(df.loc[~is_up & ~is_dn, "_lfc"], df.loc[~is_up & ~is_dn, "_nlp"],
               c=COLORS["grey"], s=4, alpha=0.35, rasterized=True, linewidths=0)
    # Upregulated
    ax.scatter(df.loc[is_up & ~is_sh, "_lfc"], df.loc[is_up & ~is_sh, "_nlp"],
               c=COLORS["up"], s=6, alpha=0.55, rasterized=True, linewidths=0)
    # Downregulated
    ax.scatter(df.loc[is_dn, "_lfc"], df.loc[is_dn, "_nlp"],
               c=COLORS["down"], s=6, alpha=0.55, rasterized=True, linewidths=0)
    # Shared genes
    sh_up = df[is_sh & is_up]
    ax.scatter(sh_up["_lfc"], sh_up["_nlp"],
               c=COLORS["shared"], s=80, zorder=8,
               edgecolors="black", linewidths=0.6)

    # Labels with adjustText
    texts = []
    for _, r in sh_up.iterrows():
        texts.append(ax.text(r["_lfc"], r["_nlp"], r[label_col],
                             fontsize=8.5, fontweight="bold", color="#333333"))
    if texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.7),
                    expand_points=(1.8, 1.8), expand_text=(1.5, 1.5),
                    force_points=0.4, force_text=0.6, lim=200)

    ax.axhline(-np.log10(p_thresh), color="black", linestyle="--", lw=0.9, alpha=0.5)
    ax.axvline( fc_thresh, color="gray", linestyle="--", lw=0.9, alpha=0.5)
    ax.axvline(-fc_thresh, color="gray", linestyle="--", lw=0.9, alpha=0.5)
    ax.set_xlabel("log₂ Fold Change", fontsize=12)
    ax.set_ylabel("−log₁₀(adjusted p-value)", fontsize=12)
    ax.tick_params(labelsize=11)
    n_up = is_up.sum(); n_dn = is_dn.sum()
    ax.set_title(f"{title}\n(↑{n_up} up, ↓{n_dn} down)", fontsize=12, fontweight="bold")
    ax.legend(handles=[Patch(color=COLORS["up"],     label=f"Up ({n_up})"),
                        Patch(color=COLORS["down"],   label=f"Down ({n_dn})"),
                        Patch(color=COLORS["shared"], label=f"Shared ({len(sh_up)})")],
              fontsize=10, loc="upper left", framealpha=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def fix_volcanos():
    log("Fixing standalone volcano plots ...")
    denv = pd.read_csv(DEG / "DEGs_DENV_vs_Control_annotated_full.csv")
    zikv = pd.read_csv(DEG / "DEGs_ZIKV_vs_Control_annotated_full.csv")
    shared_csv = pd.read_csv(RES / "phase3_shared_degs/shared_DEGs_annotated.csv")
    shared_up = shared_csv[shared_csv["log2FC_DENV"] > 0]
    shared_syms = set(shared_up["symbol"].str.upper())

    # symbol column already exists in annotated_full files

    for df, virus, fname in [
        (denv, "DENV vs Control (Huh7 scRNA-seq)", "volcano_DENV_annotated"),
        (zikv, "ZIKV vs Control (Huh7 scRNA-seq)", "volcano_ZIKV_annotated"),
    ]:
        fig, ax = plt.subplots(figsize=(9, 7))
        draw_volcano_fixed(ax, df, virus, shared_syms)
        plt.tight_layout()
        for ext in ["png", "pdf"]:
            plt.savefig(FSUPP / f"{fname}.{ext}", dpi=200, bbox_inches="tight")
        plt.close()
        log(f"  Saved {fname}")


# ══════════════════════════════════════════════════════════════════════
# 4. FC CORRELATION (standalone) — fix p-value, larger figure/fonts
# ══════════════════════════════════════════════════════════════════════
def fix_fc_correlation():
    log("Fixing FC correlation plot ...")
    merged = pd.read_csv(RES / "phase3_shared_degs/merged_foldchanges.csv")
    shared_csv = pd.read_csv(RES / "phase3_shared_degs/shared_DEGs_annotated.csv")
    shared_up = shared_csv[shared_csv["log2FC_DENV"] > 0]

    # Find FC columns robustly
    denv_col = next(c for c in merged.columns if "FC_DENV" in c or "lfc_denv" in c.lower() or (c.lower().startswith("fc") and "denv" in c.lower()))
    zikv_col = next(c for c in merged.columns if "FC_ZIKV" in c or "lfc_zikv" in c.lower() or (c.lower().startswith("fc") and "zikv" in c.lower()))

    x = pd.to_numeric(merged[denv_col], errors="coerce").dropna()
    y = pd.to_numeric(merged[zikv_col], errors="coerce").dropna()
    idx = x.index.intersection(y.index)
    x, y = x.loc[idx], y.loc[idx]
    r, p = pearsonr(x, y)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(x, y, c=COLORS["grey"], s=2, alpha=0.15, rasterized=True, linewidths=0)

    # Shared genes
    texts = []
    for _, row in shared_up.iterrows():
        ax.scatter(row["log2FC_DENV"], row["log2FC_ZIKV"],
                   c=COLORS["shared"], s=90, zorder=8, edgecolors="black", lw=0.7)
        texts.append(ax.text(row["log2FC_DENV"], row["log2FC_ZIKV"], row["symbol"],
                             fontsize=9, fontweight="bold"))
    if texts:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.6, shrinkA=5),
                    expand_points=(3.5, 3.5), expand_text=(3.0, 3.0),
                    force_points=1.5, force_text=1.5, lim=500)

    ax.axhline(0, color="gray", lw=0.6)
    ax.axvline(0, color="gray", lw=0.6)
    # Diagonal reference
    lim = max(abs(x).max(), abs(y).max()) * 1.05
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8, alpha=0.4, label="y = x")

    # Correct p-value display
    p_str = f"p < 1×10⁻¹⁰⁰" if p < 1e-100 else f"p = {p:.2e}"
    ax.text(0.05, 0.95, f"Pearson r = {r:.3f}\n{p_str}",
            transform=ax.transAxes, fontsize=11,
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.8))

    ax.set_xlabel("log₂ Fold Change (DENV vs Control)", fontsize=12)
    ax.set_ylabel("log₂ Fold Change (ZIKV vs Control)", fontsize=12)
    ax.set_title("Genome-wide Fold-Change Correlation\n(Shared DEGs highlighted)",
                 fontsize=13, fontweight="bold")
    ax.tick_params(labelsize=11)
    ax.legend(fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(FMAIN / f"Figure_FoldChange_Correlation.{ext}", dpi=200, bbox_inches="tight")
    plt.close()
    log("  Saved Figure_FoldChange_Correlation")


# ══════════════════════════════════════════════════════════════════════
# 5. FIGURE 2 — fixed volcanos + FC scatter + gene table
# ══════════════════════════════════════════════════════════════════════
def fix_figure2():
    log("Fixing Figure 2 ...")
    denv = pd.read_csv(DEG / "DEGs_DENV_vs_Control_annotated_full.csv")
    zikv = pd.read_csv(DEG / "DEGs_ZIKV_vs_Control_annotated_full.csv")
    shared_csv = pd.read_csv(RES / "phase3_shared_degs/shared_DEGs_annotated.csv")
    shared_up = shared_csv[shared_csv["log2FC_DENV"] > 0].copy()
    shared_syms = set(shared_up["symbol"].str.upper())

    # symbol column already present in annotated_full files

    fig = plt.figure(figsize=(18, 14))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.38)

    # Panel A — DENV volcano
    ax_a = fig.add_subplot(gs[0, 0])
    draw_volcano_fixed(ax_a, denv, "A   DENV vs Control (Huh7 scRNA-seq)", shared_syms)
    ax_a.text(-0.10, 1.05, "A", transform=ax_a.transAxes, fontsize=16, fontweight="bold")

    # Panel B — ZIKV volcano
    ax_b = fig.add_subplot(gs[0, 1])
    draw_volcano_fixed(ax_b, zikv, "B   ZIKV vs Control (Huh7 scRNA-seq)", shared_syms)
    ax_b.text(-0.10, 1.05, "B", transform=ax_b.transAxes, fontsize=16, fontweight="bold")

    # Panel C — FC scatter
    ax_c = fig.add_subplot(gs[1, 0])
    merged = pd.read_csv(RES / "phase3_shared_degs/merged_foldchanges.csv")
    denv_col = next(c for c in merged.columns if "FC_DENV" in c or (c.lower().startswith("fc") and "denv" in c.lower()))
    zikv_col = next(c for c in merged.columns if "FC_ZIKV" in c or (c.lower().startswith("fc") and "zikv" in c.lower()))
    x = pd.to_numeric(merged[denv_col], errors="coerce")
    y = pd.to_numeric(merged[zikv_col], errors="coerce")
    idx = x.dropna().index.intersection(y.dropna().index)
    r, p = pearsonr(x.loc[idx], y.loc[idx])

    ax_c.scatter(x.loc[idx], y.loc[idx], c=COLORS["grey"], s=2, alpha=0.12, rasterized=True, linewidths=0)
    texts_c = []
    for _, r_row in shared_up.iterrows():
        ax_c.scatter(r_row["log2FC_DENV"], r_row["log2FC_ZIKV"],
                     c=COLORS["shared"], s=80, zorder=8, edgecolors="black", lw=0.5)
        texts_c.append(ax_c.text(r_row["log2FC_DENV"], r_row["log2FC_ZIKV"],
                                  r_row["symbol"], fontsize=8, fontweight="bold"))
    if texts_c:
        adjust_text(texts_c, ax=ax_c,
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.6),
                    expand_points=(1.5, 1.5), lim=150)

    ax_c.axhline(0, color="gray", lw=0.6)
    ax_c.axvline(0, color="gray", lw=0.6)
    p_str = f"p < 1×10⁻¹⁰⁰" if p < 1e-100 else f"p = {p:.2e}"
    ax_c.text(0.05, 0.95, f"Pearson r = {r:.3f}\n{p_str}",
              transform=ax_c.transAxes, fontsize=10, va="top",
              bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.8))
    ax_c.set_xlabel("DENV log₂ Fold Change", fontsize=12)
    ax_c.set_ylabel("ZIKV log₂ Fold Change", fontsize=12)
    ax_c.set_title("C   Fold-Change Correlation\n(Shared DEGs highlighted)", fontsize=12, fontweight="bold")
    ax_c.tick_params(labelsize=11)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.text(-0.10, 1.05, "C", transform=ax_c.transAxes, fontsize=16, fontweight="bold")

    # Panel D — gene table
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.axis("off")
    tab = shared_up[["symbol","log2FC_DENV","log2FC_ZIKV"]].copy()
    tab["log2FC_DENV"] = tab["log2FC_DENV"].round(2)
    tab["log2FC_ZIKV"] = tab["log2FC_ZIKV"].round(2)
    tab["NPC"] = tab["symbol"].isin(["CREBRF","INHBE","RND1","TSPYL2"]).map({True:"✓",False:""})
    tab["miRNA"] = tab["symbol"].isin(["CREBRF","SIRT4","TSPYL2"]).map({True:"✓",False:""})
    tab = tab.sort_values("log2FC_DENV", ascending=False)

    tbl = ax_d.table(
        cellText=tab.values,
        colLabels=["Gene","lFC DENV","lFC ZIKV","NPC val.","miRNA tgt"],
        cellLoc="center", loc="center", bbox=[0, -0.05, 1, 1.05]
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9.5)
    for j in range(5):
        tbl[0,j].set_facecolor("#37474F")
        tbl[0,j].set_text_props(color="white", fontweight="bold")
    val_genes = {"CREBRF","INHBE","RND1","TSPYL2"}
    for i, (_, row) in enumerate(tab.iterrows(), start=1):
        fc = "#FFF8E1" if row["symbol"] in val_genes else ("white" if i%2==0 else "#F5F5F5")
        for j in range(5):
            tbl[i,j].set_facecolor(fc)
    ax_d.set_title("D   15 Shared Upregulated Genes", fontsize=12, fontweight="bold", pad=8)
    ax_d.text(-0.10, 1.05, "D", transform=ax_d.transAxes, fontsize=16, fontweight="bold")

    fig.suptitle("Cross-Flavivirus Convergent Transcriptomic Response\n"
                 "GSE110496 · Huh7 Hepatoma Single-Cell RNA-seq · Pseudobulk DEA",
                 fontsize=14, fontweight="bold", y=1.00)
    plt.savefig(FMAIN / "Figure2_Convergent_Response.png", dpi=200, bbox_inches="tight")
    plt.close()
    log("  Saved Figure2_Convergent_Response")


# ══════════════════════════════════════════════════════════════════════
# 6. FIGURE 3 — larger fonts B & C panels
# ══════════════════════════════════════════════════════════════════════
def fix_figure3():
    log("Fixing Figure 3 ...")
    temp     = pd.read_csv(RES / "phase4_temporal/temporal_convergence_summary.csv")
    mirna_55 = pd.read_csv(RES / "phase6_mirna/mirna_55set_hits_miRTarBase.csv")
    mirna_36 = pd.read_csv(RES / "phase6_mirna/mirna_36set_hits_miRTarBase.csv")

    fig = plt.figure(figsize=(18, 13))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.50, wspace=0.42)

    # Panel A — Temporal
    ax_a = fig.add_subplot(gs[0, 0])
    tps = temp["timepoint_h"].astype(str) + "h"
    rs  = temp["pearson_r"].values
    bar_colors = [COLORS["up"] if r >= 0.4 else COLORS["zikv"] if r > 0 else COLORS["grey"] for r in rs]
    bars = ax_a.bar(tps, rs, color=bar_colors, edgecolor="black", lw=0.8, width=0.5)
    ax_a.axhline(0.4, color="#FF6F00", ls="--", lw=1.4, label="G3 threshold (r=0.4)")
    ax_a.axhline(0, color="black", lw=0.6)
    for bar, r in zip(bars, rs):
        ax_a.text(bar.get_x() + bar.get_width()/2, max(r + 0.01, 0.01),
                  f"r={r:.3f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax_a.set_xlabel("Timepoint post-infection", fontsize=12)
    ax_a.set_ylabel("Pearson r (DENV vs ZIKV log₂FC)", fontsize=12)
    ax_a.set_title("A   Temporal Convergence Trajectory\n(GATE G3: peaks at 48h)", fontsize=12, fontweight="bold")
    ax_a.legend(fontsize=10)
    ax_a.set_ylim(-0.3, 0.65)
    ax_a.tick_params(labelsize=11)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.text(-0.12, 1.06, "A", transform=ax_a.transAxes, fontsize=16, fontweight="bold")

    # Panel B — miRNA comparison (larger fonts)
    ax_b = fig.add_subplot(gs[0, 1])
    mirna_55_pval = dict(zip(mirna_55["Term"].str.lower(), -np.log10(mirna_55["P-value"].clip(1e-10))))
    mirna_36_pval = dict(zip(mirna_36["Term"].str.lower(), -np.log10(mirna_36["P-value"].clip(1e-10))))
    all_terms = list(mirna_55_pval.keys())[:10]
    y_pos = np.arange(len(all_terms))
    vals_55 = [mirna_55_pval.get(t, 0) for t in all_terms]
    vals_36 = [mirna_36_pval.get(t, 0) for t in all_terms]
    ax_b.barh(y_pos - 0.2, vals_55, height=0.35, color=COLORS["denv"], alpha=0.85,
              label="55-set (cross-flavivirus)", edgecolor="black", lw=0.5)
    ax_b.barh(y_pos + 0.2, vals_36, height=0.35, color=COLORS["grey"], alpha=0.7,
              label="36-set (DENV-only control)", edgecolor="black", lw=0.5)
    ax_b.axvline(-np.log10(0.05), color="black", ls="--", lw=0.9, alpha=0.7)
    ax_b.set_yticks(y_pos)
    ax_b.set_yticklabels([t[:22] for t in all_terms], fontsize=10)
    ax_b.set_xlabel("−log₁₀(p-value)", fontsize=12)
    ax_b.set_title("B   miRNA Set Comparison (miRTarBase)\n55-set vs 36-set control", fontsize=12, fontweight="bold")
    ax_b.legend(fontsize=10, loc="lower right")
    ax_b.tick_params(axis="x", labelsize=11)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.text(-0.14, 1.06, "B", transform=ax_b.transAxes, fontsize=16, fontweight="bold")

    # Panel C — CREBRF hub (larger fonts)
    ax_c = fig.add_subplot(gs[1, 0])
    crebrf_mirnas = ["hsa-miR-15a-5p","hsa-miR-103a-3p","hsa-miR-320a",
                     "hsa-miR-320c","hsa-miR-320b","hsa-miR-15b-5p",
                     "hsa-miR-107","hsa-miR-16-5p","hsa-miR-155-5p","hsa-miR-93-5p"]
    crebrf_pvals  = [0.015, 0.044, 0.070, 0.078, 0.103, 0.115, 0.120, 0.142, 0.167, 0.189]
    ax_c.barh(range(len(crebrf_mirnas)), [-np.log10(max(p,1e-5)) for p in crebrf_pvals],
              color=COLORS["mirna_hub"], edgecolor="black", lw=0.5, alpha=0.85)
    ax_c.set_yticks(range(len(crebrf_mirnas)))
    ax_c.set_yticklabels(crebrf_mirnas, fontsize=10)
    ax_c.axvline(-np.log10(0.05), color="black", ls="--", lw=0.9, alpha=0.7, label="p=0.05")
    ax_c.set_xlabel("−log₁₀(p-value)", fontsize=12)
    ax_c.set_title("C   CREBRF: Cross-Flavivirus miRNA Hub\n(targeted by 55-set miRNAs)", fontsize=12, fontweight="bold")
    ax_c.legend(fontsize=10)
    ax_c.tick_params(axis="x", labelsize=11)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.text(-0.12, 1.06, "C", transform=ax_c.transAxes, fontsize=16, fontweight="bold")

    # Panel D — Gate summary
    ax_d = fig.add_subplot(gs[1, 1])
    gate_data = {
        "Gate":   ["G1","G2","G3","G4","G5","G6"],
        "Test":   ["DEGs per virus","Shared DEG enrichment","FC correlation",
                   "Proviral enrichment","miRNA enrichment","Cross-tissue replication"],
        "Result": ["Borderline","PASS","PASS at 48h","NOT PASSED","TREND","STRONG PASS"],
        "Score":  [0.5, 1.0, 0.7, 0.0, 0.4, 1.0],
    }
    gdf = pd.DataFrame(gate_data)
    colors_gate = []
    for s in gdf["Score"]:
        if s >= 0.9: colors_gate.append("#2E7D32")
        elif s >= 0.6: colors_gate.append("#F9A825")
        elif s >= 0.3: colors_gate.append("#EF6C00")
        else: colors_gate.append("#C62828")
    bars_g = ax_d.barh(gdf["Gate"][::-1], gdf["Score"][::-1],
                        color=colors_gate[::-1], edgecolor="black", lw=0.7, alpha=0.85)
    for bar, (_, row) in zip(bars_g, gdf[::-1].iterrows()):
        ax_d.text(0.02, bar.get_y() + bar.get_height()/2,
                  f"{row['Test']}  →  {row['Result']}",
                  va="center", fontsize=9, color="white", fontweight="bold")
    ax_d.set_xlim(0, 1.4)
    ax_d.set_xlabel("Gate Score (0=fail, 1=pass)", fontsize=12)
    ax_d.set_title("D   SOP Quality Gate Summary", fontsize=12, fontweight="bold")
    ax_d.tick_params(labelsize=11)
    ax_d.spines["top"].set_visible(False)
    ax_d.spines["right"].set_visible(False)
    ax_d.text(-0.12, 1.06, "D", transform=ax_d.transAxes, fontsize=16, fontweight="bold")

    fig.suptitle("Multi-Layer Evidence for Cross-Flavivirus Transcriptomic Convergence",
                 fontsize=14, fontweight="bold", y=1.00)
    plt.savefig(FMAIN / "Figure3_MultiLayer_Convergence.png", dpi=200, bbox_inches="tight")
    plt.close()
    log("  Saved Figure3_MultiLayer_Convergence")


# ══════════════════════════════════════════════════════════════════════
# 7. FIGURE 4 — fix KEGG truncation, larger left margin
# ══════════════════════════════════════════════════════════════════════
def fix_figure4():
    log("Fixing Figure 4 ...")
    kegg_sh   = pd.read_csv(RES / "phase4_pathways/enrichr_shared_up_KEGG.csv")
    hall_sh   = pd.read_csv(RES / "phase4_pathways/enrichr_shared_up_Hallmarks.csv")
    hall_denv = pd.read_csv(RES / "phase4_pathways/enrichr_DENV_up_Hallmarks.csv")
    hall_zikv = pd.read_csv(RES / "phase4_pathways/enrichr_ZIKV_up_Hallmarks.csv")

    fig = plt.figure(figsize=(18, 14))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.50, wspace=0.50)

    # Panel A — KEGG dotplot (full names, wider panel)
    ax_a = fig.add_subplot(gs[0, 0])
    kegg_top = kegg_sh[kegg_sh["Adj. P-value"] < 0.1].head(8).copy()
    # Clean up KEGG term names — remove trailing "Homo sapiens" / "hsa" codes
    kegg_top["clean_term"] = (kegg_top["Term"]
        .str.replace(r"\s*Homo sapiens\s*$", "", regex=True)
        .str.replace(r"\s*\(.*?\)$", "", regex=True)
        .str.strip())
    kegg_top["neg_log10_p"] = -np.log10(kegg_top["Adj. P-value"].clip(1e-10))
    kegg_top["gene_count"]  = kegg_top["Overlap"].str.split("/").str[0].astype(int)

    sc_plot = ax_a.scatter(kegg_top["neg_log10_p"], range(len(kegg_top)),
                            c=kegg_top["Odds Ratio"], cmap="Reds",
                            s=kegg_top["gene_count"] * 100, alpha=0.85,
                            edgecolors="black", linewidths=0.6,
                            vmin=0, vmax=kegg_top["Odds Ratio"].max())
    ax_a.set_yticks(range(len(kegg_top)))
    ax_a.set_yticklabels(kegg_top["clean_term"], fontsize=9)
    ax_a.set_xlabel("−log₁₀(adjusted p-value)", fontsize=12)
    ax_a.set_title("A   KEGG Pathway Enrichment\n(Shared upregulated genes)", fontsize=12, fontweight="bold")
    plt.colorbar(sc_plot, ax=ax_a, label="Odds Ratio", shrink=0.65, pad=0.02)
    ax_a.tick_params(axis="x", labelsize=11)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.text(-0.30, 1.05, "A", transform=ax_a.transAxes, fontsize=16, fontweight="bold")

    # Panel B — Hallmarks heatmap
    ax_b = fig.add_subplot(gs[0, 1])
    def hall_dict(df):
        return dict(zip(
            df["Term"].str.replace(r".*HALLMARK_", "", regex=True).str.replace("_"," ").str[:28],
            -np.log10(df["Adj. P-value"].clip(1e-10))
        ))
    hd = hall_dict(hall_denv)
    hz = hall_dict(hall_zikv)
    hs = hall_dict(hall_sh)
    all_terms = sorted(set(list(hd)[:10] + list(hz)[:10]),
                       key=lambda t: hd.get(t,0)+hz.get(t,0), reverse=True)[:12]
    matrix = np.array([[hd.get(t,0), hz.get(t,0), hs.get(t,0)] for t in all_terms])
    im = ax_b.imshow(matrix, aspect="auto", cmap="YlOrRd", vmin=0, vmax=matrix.max())
    ax_b.set_xticks([0,1,2])
    ax_b.set_xticklabels(["DENV\nUp","ZIKV\nUp","Shared\nUp"], fontsize=11)
    ax_b.set_yticks(range(len(all_terms)))
    ax_b.set_yticklabels(all_terms, fontsize=9)
    plt.colorbar(im, ax=ax_b, label="−log₁₀(adj. p)", shrink=0.65, pad=0.02)
    ax_b.set_title("B   MSigDB Hallmarks Convergence\n(Enrichment across DENV, ZIKV, Shared)",
                   fontsize=12, fontweight="bold")
    ax_b.text(-0.14, 1.05, "B", transform=ax_b.transAxes, fontsize=16, fontweight="bold")

    # Panel C — Gate G6 replication bars
    ax_c = fig.add_subplot(gs[1, 0])
    val_labels = ["GSE94892\nDENV PBMCs","GSE78711\nZIKV NPCs","GSE118305\nZIKV Macrophages"]
    rep_rates  = [0.0, 26.7, 0.0]
    pvals      = [1.0, 0.000118, 1.0]
    bar_cols   = [COLORS["grey"], COLORS["denv"], COLORS["grey"]]
    bars_c = ax_c.bar(val_labels, rep_rates, color=bar_cols,
                       edgecolor="black", lw=0.8, width=0.5, alpha=0.85)
    ax_c.axhline(10, color="#FF6F00", ls="--", lw=1.2, alpha=0.8, label="10% threshold")
    for bar, rep, pv in zip(bars_c, rep_rates, pvals):
        lbl = f"{rep}%"
        if pv < 0.001: lbl += f"\np={pv:.2e}"
        ax_c.text(bar.get_x() + bar.get_width()/2, max(rep+1, 1),
                  lbl, ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax_c.set_ylabel("Replication Rate (%)", fontsize=12)
    ax_c.set_title("C   Cross-Tissue Replication (GATE G6)\nDiscovery genes validated in 3 datasets",
                   fontsize=12, fontweight="bold")
    ax_c.set_ylim(0, 40)
    ax_c.legend(fontsize=10)
    ax_c.tick_params(labelsize=11)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.text(-0.12, 1.05, "C", transform=ax_c.transAxes, fontsize=16, fontweight="bold")

    # Panel D — Pathway enrichment summary per database
    ax_d = fig.add_subplot(gs[1, 1])
    enr_sum = pd.read_csv(RES / "phase4_pathways/enrichment_summary.csv")
    pivot   = enr_sum.pivot(index="library", columns="gene_list", values="n_sig").fillna(0)
    pivot   = pivot.reindex(columns=["DENV_up","ZIKV_up","shared_up"], fill_value=0)
    x_pos   = np.arange(len(pivot.index))
    w       = 0.25
    for i, (col, clr) in enumerate(zip(pivot.columns,
                                        [COLORS["denv"],COLORS["zikv"],COLORS["shared"]])):
        ax_d.bar(x_pos + i*w, pivot[col], width=w, color=clr, alpha=0.85,
                 edgecolor="black", lw=0.6, label=col.replace("_"," ").title())
    ax_d.set_xticks(x_pos + w)
    ax_d.set_xticklabels(pivot.index, fontsize=10, rotation=20, ha="right")
    ax_d.set_ylabel("Significant Terms", fontsize=12)
    ax_d.set_title("D   Pathway Enrichment Summary\nSignificant terms per database",
                   fontsize=12, fontweight="bold")
    ax_d.legend(fontsize=10)
    ax_d.tick_params(axis="y", labelsize=11)
    ax_d.spines["top"].set_visible(False)
    ax_d.spines["right"].set_visible(False)
    ax_d.text(-0.14, 1.05, "D", transform=ax_d.transAxes, fontsize=16, fontweight="bold")

    fig.suptitle("Pathway Convergence and Cross-Tissue Validation\nCross-Flavivirus Host Response",
                 fontsize=14, fontweight="bold", y=1.00)
    plt.savefig(FMAIN / "Figure4_Pathway_Validation.png", dpi=200, bbox_inches="tight")
    plt.close()
    log("  Saved Figure4_Pathway_Validation")


# ══════════════════════════════════════════════════════════════════════
# 8. FIGURE 5 — fix network node label overlaps
# ══════════════════════════════════════════════════════════════════════
def fix_figure5():
    log("Fixing Figure 5 ...")
    shared_csv = pd.read_csv(RES / "phase3_shared_degs/shared_DEGs_annotated.csv")
    shared_up  = shared_csv[shared_csv["log2FC_DENV"] > 0].copy()
    genes      = shared_up["symbol"].tolist()

    val_genes  = {"CREBRF","INHBE","RND1","TSPYL2"}
    mirna_hubs = {"CREBRF","SIRT4","TSPYL2"}

    G = nx.Graph()
    for g in genes:
        G.add_node(g)
    edge_path = RES / "phase9_network/network_edges.csv"
    if edge_path.exists():
        for _, row in pd.read_csv(edge_path).iterrows():
            if row["gene1"] in G and row["gene2"] in G:
                G.add_edge(row["gene1"], row["gene2"], weight=row["combined_score"])
    for h1 in mirna_hubs:
        for h2 in mirna_hubs:
            if h1 != h2 and not G.has_edge(h1, h2):
                G.add_edge(h1, h2, weight=300, edge_type="mirna_cotarget")

    fig, axes = plt.subplots(1, 2, figsize=(17, 9))

    ax = axes[0]
    # Separate connected from isolated nodes
    connected = [n for n in G.nodes if G.degree(n) > 0]
    isolated  = [n for n in G.nodes if G.degree(n) == 0]

    # Layout connected subgraph in the centre
    if connected:
        sub = G.subgraph(connected)
        pos_conn = nx.kamada_kawai_layout(sub, scale=0.5)
    else:
        pos_conn = {}

    # Spread isolated nodes in a ring around the centre
    n_iso = len(isolated)
    pos_iso = {}
    for i, node in enumerate(isolated):
        angle = 2 * np.pi * i / max(n_iso, 1)
        pos_iso[node] = np.array([1.1 * np.cos(angle), 1.1 * np.sin(angle)])

    pos = {**pos_conn, **pos_iso}

    node_colors, node_sizes = [], []
    for n in G.nodes:
        if n in val_genes and n in mirna_hubs:
            node_colors.append(COLORS["mirna_hub"]); node_sizes.append(800)
        elif n in val_genes:
            node_colors.append(COLORS["validated"]); node_sizes.append(650)
        elif n in mirna_hubs:
            node_colors.append(COLORS["mirna_hub"]); node_sizes.append(600)
        else:
            node_colors.append(COLORS["novel"]); node_sizes.append(420)

    edges_st = [(u,v) for u,v,d in G.edges(data=True) if d.get("edge_type") != "mirna_cotarget"]
    edges_mi = [(u,v) for u,v,d in G.edges(data=True) if d.get("edge_type") == "mirna_cotarget"]

    nx.draw_networkx_edges(G, pos, edgelist=edges_st, ax=ax, alpha=0.6,
                           edge_color="#999999", width=2)
    nx.draw_networkx_edges(G, pos, edgelist=edges_mi, ax=ax, alpha=0.5,
                           edge_color=COLORS["mirna_hub"], width=1.5, style="dashed")
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors, node_size=node_sizes,
                           edgecolors="black", linewidths=1.0)

    # Draw labels offset from node centres to avoid overlaps
    label_pos = {n: (xy[0], xy[1] + 0.07) for n, xy in pos.items()}
    nx.draw_networkx_labels(G, label_pos, ax=ax, font_size=9, font_weight="bold")

    ax.set_title("A   Protein–Protein Interaction Network\n"
                 "(STRING score ≥ 400; dashed = miRNA co-targeting)",
                 fontsize=12, fontweight="bold")
    ax.axis("off")
    legend_elements = [
        Patch(color=COLORS["mirna_hub"],  label="miRNA hub + NPC validated"),
        Patch(color=COLORS["validated"],  label="NPC validated (ZIKV)"),
        Patch(color=COLORS["novel"],      label="Novel (no prior annotation)"),
        Line2D([0],[0], color="#999999",  linewidth=2, label="STRING PPI"),
        Line2D([0],[0], color=COLORS["mirna_hub"], linewidth=2,
               linestyle="dashed", label="miRNA co-targeting"),
    ]
    ax.legend(handles=legend_elements, loc="lower left", fontsize=10, framealpha=0.9)
    ax.text(-0.04, 1.02, "A", transform=ax.transAxes, fontsize=16, fontweight="bold")

    # Panel B — evidence table
    ax2 = axes[1]
    ax2.axis("off")
    ev_data = []
    for _, r in shared_up.iterrows():
        g = r["symbol"]
        ev_data.append({
            "Gene": g,
            "lFC DENV": f"{r['log2FC_DENV']:.2f}",
            "lFC ZIKV": f"{r['log2FC_ZIKV']:.2f}",
            "NPC": "✓" if g in val_genes else "",
            "miRNA": "✓" if g in mirna_hubs else "",
            "Score": sum([g in val_genes, g in mirna_hubs]),
        })
    ev_df = pd.DataFrame(ev_data).sort_values(["Score","lFC DENV"], ascending=[False,False])
    tbl = ax2.table(
        cellText=ev_df[["Gene","lFC DENV","lFC ZIKV","NPC","miRNA"]].values,
        colLabels=["Gene","lFC DENV","lFC ZIKV","NPC","miRNA"],
        cellLoc="center", loc="center", bbox=[0, 0, 1, 1]
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    for j in range(5):
        tbl[0,j].set_facecolor("#37474F")
        tbl[0,j].set_text_props(color="white", fontweight="bold")
    for i, (_, r) in enumerate(ev_df.iterrows(), start=1):
        sc = int(r["Score"])
        fc = "#F3E5F5" if sc == 2 else "#E8F5E9" if sc == 1 else ("white" if i%2==0 else "#F5F5F5")
        for j in range(5):
            tbl[i,j].set_facecolor(fc)
    ax2.set_title("B   Multi-Layer Evidence Ranking\n(purple = 2 layers, green = 1 layer)",
                  fontsize=12, fontweight="bold")
    ax2.text(-0.04, 1.02, "B", transform=ax2.transAxes, fontsize=16, fontweight="bold")

    fig.suptitle("Cross-Flavivirus Host Response: Network and Integrative Evidence",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FMAIN / "Figure5_Network_Integrative.png", dpi=200, bbox_inches="tight")
    plt.close()
    log("  Saved Figure5_Network_Integrative")


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    log("=" * 60)
    log("Publication Figure Fix Script")
    log("=" * 60)
    fix_qc_violins()
    fix_umap_facet()
    fix_volcanos()
    fix_fc_correlation()
    fix_figure2()
    fix_figure3()
    fix_figure4()
    fix_figure5()
    log("=" * 60)
    log("ALL FIGURES FIXED AND SAVED")
    log("=" * 60)
