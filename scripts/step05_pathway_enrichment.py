"""
Step 05: Pathway Enrichment Analysis (Enrichr via gseapy)
- GO BP / MF / CC, KEGG, Reactome, MSigDB Hallmarks
- Run for DENV up/down, ZIKV up/down, shared genes
- Generate bar plots, dot plots, comparison figures
Checkpoint-based: safe to restart.
"""

import json, time, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import gseapy as gp
from pathlib import Path

warnings.filterwarnings("ignore")

# ─── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR  = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
DEG_DIR   = BASE_DIR / "01_processed_data" / "deg_tables"
RES_DIR   = BASE_DIR / "03_results" / "phase4_pathways"
SHARED_DIR= BASE_DIR / "03_results" / "phase3_shared_degs"
FIG_MAIN  = BASE_DIR / "04_figures" / "main"
FIG_SUPP  = BASE_DIR / "04_figures" / "supplementary"
CKPT_DIR  = BASE_DIR / "checkpoints"
LOG_FILE  = BASE_DIR / "logs" / "step05_pathway.log"

for d in [RES_DIR, FIG_MAIN, FIG_SUPP, CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE = CKPT_DIR / "step05_checkpoint.json"

# ─── Enrichr libraries to query ───────────────────────────────────────────────
LIBRARIES = {
    "GO_BP"     : "GO_Biological_Process_2023",
    "GO_MF"     : "GO_Molecular_Function_2023",
    "GO_CC"     : "GO_Cellular_Component_2023",
    "KEGG"      : "KEGG_2021_Human",
    "Reactome"  : "Reactome_Pathways_2024",
    "Hallmarks" : "MSigDB_Hallmark_2020",
}

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


# ─── Run Enrichr for one gene list ────────────────────────────────────────────
def run_enrichr(gene_list: list, label: str, lib_key: str, lib_name: str,
                out_dir: Path) -> pd.DataFrame | None:
    out_csv = out_dir / f"enrichr_{label}_{lib_key}.csv"
    if out_csv.exists():
        log(f"  ✓ {label}/{lib_key} cached")
        return pd.read_csv(out_csv)

    if len(gene_list) < 3:
        log(f"  ! {label}/{lib_key} skipped — only {len(gene_list)} genes")
        return None

    try:
        enr = gp.enrichr(
            gene_list   = gene_list,
            gene_sets   = lib_name,
            organism    = "human",
            outdir      = None,
            verbose     = False,
        )
        df = enr.results
        if df is None or df.empty:
            log(f"  ! {label}/{lib_key} — no results")
            return None

        # Standardise column names across gseapy versions
        df.columns = [c.strip() for c in df.columns]
        if "Adjusted P-value" in df.columns:
            df = df.rename(columns={"Adjusted P-value": "Adj. P-value"})
        if "P-value" in df.columns:
            df = df.rename(columns={"P-value": "P.value"})

        df["neg_log10_padj"] = -np.log10(df["Adj. P-value"].clip(lower=1e-300))
        df["gene_set_lib"]   = lib_key
        df["query_label"]    = label
        df.to_csv(out_csv, index=False)
        log(f"  → {label}/{lib_key}: {len(df)} terms returned")
        return df

    except Exception as e:
        log(f"  ERROR {label}/{lib_key}: {e}")
        return None


# ─── Run all libraries for a gene list ────────────────────────────────────────
def enrich_all_libs(gene_list: list, label: str, out_dir: Path) -> dict:
    results = {}
    for lib_key, lib_name in LIBRARIES.items():
        df = run_enrichr(gene_list, label, lib_key, lib_name, out_dir)
        if df is not None and not df.empty:
            results[lib_key] = df
        time.sleep(0.4)   # be polite to Enrichr API
    return results


# ─── Bar plot: top N enriched terms ───────────────────────────────────────────
def plot_barplot(df: pd.DataFrame, title: str, out_path: Path,
                 color: str = "#4575B4", top_n: int = 15):
    df = df.copy()

    # Determine p-value column
    padj_col = next((c for c in ["Adj. P-value", "Adjusted P-value", "FDR"]
                     if c in df.columns), None)
    if padj_col is None:
        log(f"    ! No padj column found in {title}, skipping bar plot")
        return
    df = df[df[padj_col] < 0.05].copy()
    if df.empty:
        log(f"    ! No significant terms for {title}")
        return

    df["neg_log10_padj"] = -np.log10(df[padj_col].clip(lower=1e-300))

    # Shorten term names
    term_col = "Term" if "Term" in df.columns else df.columns[0]
    df[term_col] = df[term_col].str.replace(r"\(GO:\d+\)", "", regex=True).str.strip()
    df[term_col] = df[term_col].str[:60]

    df = df.nsmallest(top_n, padj_col)
    df = df.sort_values("neg_log10_padj", ascending=True)

    fig, ax = plt.subplots(figsize=(9, max(4, len(df) * 0.4)))
    bars = ax.barh(df[term_col], df["neg_log10_padj"], color=color, edgecolor="white")
    ax.axvline(x=-np.log10(0.05), color="red", linestyle="--", linewidth=1,
               label="FDR = 0.05")
    ax.set_xlabel("-log₁₀(Adj. P-value)", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close()
    log(f"    Bar plot → {out_path.name}")


# ─── Dot plot: compare DENV vs ZIKV for a library ────────────────────────────
def plot_comparison_dotplot(denv_df: pd.DataFrame, zikv_df: pd.DataFrame,
                            lib_key: str, out_path: Path, top_n: int = 20):
    padj_col = next((c for c in ["Adj. P-value", "Adjusted P-value", "FDR"]
                     if c in denv_df.columns), None)
    if padj_col is None:
        return
    term_col = "Term" if "Term" in denv_df.columns else denv_df.columns[0]

    def prep(df, virus):
        d = df[df[padj_col] < 0.05][[term_col, padj_col, "Overlap"]].copy()
        d[term_col] = d[term_col].str.replace(r"\(GO:\d+\)", "", regex=True).str.strip()
        d[term_col] = d[term_col].str[:55]
        d["neg_log10_padj"] = -np.log10(d[padj_col].clip(lower=1e-300))
        d["overlap_n"] = d["Overlap"].str.split("/").str[0].astype(int)
        d["virus"] = virus
        return d

    denv_p = prep(denv_df, "DENV")
    zikv_p = prep(zikv_df, "ZIKV")

    # Take top terms from each, then union
    top_terms = set(denv_p.nsmallest(top_n, padj_col)[term_col].tolist()) | \
                set(zikv_p.nsmallest(top_n, padj_col)[term_col].tolist())

    combined = pd.concat([denv_p, zikv_p])
    combined = combined[combined[term_col].isin(top_terms)]

    if combined.empty:
        return

    pivot_val  = combined.pivot_table(index=term_col, columns="virus",
                                       values="neg_log10_padj", aggfunc="max").fillna(0)
    pivot_size = combined.pivot_table(index=term_col, columns="virus",
                                       values="overlap_n", aggfunc="max").fillna(0)

    # Sort by DENV significance
    pivot_val = pivot_val.sort_values("DENV" if "DENV" in pivot_val.columns else pivot_val.columns[0],
                                      ascending=False).head(top_n)
    pivot_size = pivot_size.loc[pivot_val.index]

    viruses = [v for v in ["DENV", "ZIKV"] if v in pivot_val.columns]
    n_terms = len(pivot_val)

    fig, ax = plt.subplots(figsize=(7, max(5, n_terms * 0.5)))
    colors  = {"DENV": "#E41A1C", "ZIKV": "#377EB8"}
    x_pos   = {"DENV": 0, "ZIKV": 1}

    for virus in viruses:
        for i, term in enumerate(pivot_val.index):
            val  = pivot_val.loc[term, virus]
            size = pivot_size.loc[term, virus]
            if val > 0:
                ax.scatter(x_pos[virus], i,
                           s=size * 20, c=colors[virus], alpha=0.75,
                           edgecolors="grey", linewidths=0.5, zorder=3)
                ax.text(x_pos[virus], i, f"{val:.1f}",
                        ha="center", va="center", fontsize=6, color="white",
                        fontweight="bold")

    ax.set_yticks(range(n_terms))
    ax.set_yticklabels(pivot_val.index, fontsize=9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["DENV", "ZIKV"], fontsize=11, fontweight="bold")
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.5, n_terms - 0.5)
    ax.invert_yaxis()
    ax.set_title(f"{lib_key} enrichment — DENV vs ZIKV\n(bubble = overlap genes; text = -log₁₀ FDR)",
                 fontsize=11, fontweight="bold")
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    legend_patches = [mpatches.Patch(color=colors[v], label=v) for v in viruses]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close()
    log(f"    Comparison dot plot → {out_path.name}")


# ─── Hallmarks dot plot (single virus) ───────────────────────────────────────
def plot_hallmarks(df: pd.DataFrame, title: str, out_path: Path, top_n: int = 15):
    padj_col = next((c for c in ["Adj. P-value", "Adjusted P-value", "FDR"]
                     if c in df.columns), None)
    if padj_col is None:
        return
    term_col = "Term" if "Term" in df.columns else df.columns[0]

    df = df[df[padj_col] < 0.05].copy()
    if df.empty:
        log(f"    ! No significant Hallmarks for {title}")
        return

    df["neg_log10_padj"] = -np.log10(df[padj_col].clip(lower=1e-300))
    df["overlap_n"] = df["Overlap"].str.split("/").str[0].astype(int)
    df = df.nsmallest(top_n, padj_col).sort_values("neg_log10_padj")

    df[term_col] = df[term_col].str.replace("HALLMARK_", "").str.replace("_", " ").str.title()

    fig, ax = plt.subplots(figsize=(8, max(4, len(df) * 0.45)))
    scatter = ax.scatter(
        df["neg_log10_padj"], df[term_col],
        s=df["overlap_n"] * 15, c=df["neg_log10_padj"],
        cmap="Reds", edgecolors="grey", linewidths=0.5, zorder=3
    )
    ax.axvline(x=-np.log10(0.05), color="navy", linestyle="--", linewidth=1)
    ax.set_xlabel("-log₁₀(Adj. P-value)", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6)
    cbar.set_label("-log₁₀(FDR)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close()
    log(f"    Hallmarks dot plot → {out_path.name}")


# ─── Save combined result table ────────────────────────────────────────────────
def save_combined(results: dict, label: str, out_dir: Path):
    frames = [df.assign(library=k) for k, df in results.items()]
    if not frames:
        return
    combined = pd.concat(frames, ignore_index=True)
    out = out_dir / f"enrichment_all_{label}.csv"
    combined.to_csv(out, index=False)
    log(f"  Combined enrichment table → {out.name}")


# ─── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    ckpt = load_ckpt()

    log("=" * 60)
    log("Step 05: Pathway Enrichment Analysis (Enrichr)")
    log("=" * 60)

    # ── Load annotated DEG tables ──────────────────────────────────────────────
    denv_ann = pd.read_csv(DEG_DIR / "DEGs_DENV_vs_Control_annotated.csv")
    zikv_ann = pd.read_csv(DEG_DIR / "DEGs_ZIKV_vs_Control_annotated.csv")
    shared   = pd.read_csv(SHARED_DIR / "shared_DEGs_annotated.csv")

    log(f"DENV annotated DEGs  : {len(denv_ann)}")
    log(f"ZIKV annotated DEGs  : {len(zikv_ann)}")
    log(f"Shared upregulated   : {len(shared)}")

    # Gene lists
    fc_col = "log2FC" if "log2FC" in denv_ann.columns else "log2FoldChange"

    denv_up   = denv_ann[denv_ann[fc_col] > 0]["symbol"].dropna().unique().tolist()
    denv_dn   = denv_ann[denv_ann[fc_col] < 0]["symbol"].dropna().unique().tolist()
    zikv_up   = zikv_ann[zikv_ann[fc_col] > 0]["symbol"].dropna().unique().tolist()
    zikv_dn   = zikv_ann[zikv_ann[fc_col] < 0]["symbol"].dropna().unique().tolist()
    denv_all  = denv_ann["symbol"].dropna().unique().tolist()
    zikv_all  = zikv_ann["symbol"].dropna().unique().tolist()
    shared_up = shared["symbol"].dropna().unique().tolist()

    log(f"\nGene list sizes:")
    log(f"  DENV up={len(denv_up)}  dn={len(denv_dn)}  all={len(denv_all)}")
    log(f"  ZIKV up={len(zikv_up)}  dn={len(zikv_dn)}  all={len(zikv_all)}")
    log(f"  Shared up={len(shared_up)}")

    gene_lists = {
        "DENV_up"   : denv_up,
        "DENV_down" : denv_dn,
        "ZIKV_up"   : zikv_up,
        "ZIKV_down" : zikv_dn,
        "shared_up" : shared_up,
    }

    # ── A: Run Enrichr for all gene lists ─────────────────────────────────────
    log("\n── A: Running Enrichr ──────────────────────────────────────────")
    all_results = {}
    for label, genes in gene_lists.items():
        if ckpt.get(f"enrichr_{label}_done"):
            log(f"  ✓ {label} enrichment cached")
            # Reload from disk
            res = {}
            for lib_key in LIBRARIES:
                fp = RES_DIR / f"enrichr_{label}_{lib_key}.csv"
                if fp.exists():
                    res[lib_key] = pd.read_csv(fp)
            all_results[label] = res
        else:
            log(f"\n  Running Enrichr for {label} ({len(genes)} genes) ...")
            res = enrich_all_libs(genes, label, RES_DIR)
            all_results[label] = res
            save_combined(res, label, RES_DIR)
            ckpt[f"enrichr_{label}_done"] = True
            save_ckpt(ckpt)

    # ── B: Bar plots for each virus / direction / library ─────────────────────
    log("\n── B: Bar plots ────────────────────────────────────────────────")
    colors = {
        "DENV_up":   "#E41A1C",
        "DENV_down": "#4575B4",
        "ZIKV_up":   "#FF7F00",
        "ZIKV_down": "#984EA3",
        "shared_up": "#33A02C",
    }

    for label, res in all_results.items():
        for lib_key, df in res.items():
            if df is None or df.empty:
                continue
            title    = f"{label.replace('_',' ').title()} — {lib_key}"
            out_path = FIG_SUPP / f"barplot_{label}_{lib_key}"
            plot_barplot(df, title, out_path, color=colors.get(label, "#4575B4"))

    # ── C: DENV vs ZIKV comparison dot plots (upregulated genes) ─────────────
    log("\n── C: Comparison dot plots (DENV_up vs ZIKV_up) ────────────────")
    if ckpt.get("comparison_plots_done"):
        log("  ✓ Comparison plots cached")
    else:
        for lib_key in LIBRARIES:
            denv_res_lib = all_results.get("DENV_up", {}).get(lib_key)
            zikv_res_lib = all_results.get("ZIKV_up", {}).get(lib_key)
            if denv_res_lib is not None and zikv_res_lib is not None:
                out_path = FIG_MAIN / f"Figure_Enrichment_Comparison_{lib_key}"
                plot_comparison_dotplot(denv_res_lib, zikv_res_lib, lib_key, out_path)
        ckpt["comparison_plots_done"] = True
        save_ckpt(ckpt)

    # ── D: Hallmarks dot plots ────────────────────────────────────────────────
    log("\n── D: MSigDB Hallmarks dot plots ───────────────────────────────")
    if ckpt.get("hallmarks_done"):
        log("  ✓ Hallmarks plots cached")
    else:
        for label in ["DENV_up", "ZIKV_up", "shared_up"]:
            df = all_results.get(label, {}).get("Hallmarks")
            if df is not None and not df.empty:
                title    = f"MSigDB Hallmarks — {label.replace('_',' ').title()}"
                out_path = FIG_MAIN / f"Figure_Hallmarks_{label}"
                plot_hallmarks(df, title, out_path)
        ckpt["hallmarks_done"] = True
        save_ckpt(ckpt)

    # ── E: Shared gene enrichment summary ────────────────────────────────────
    log("\n── E: Shared gene enrichment summary ───────────────────────────")
    shared_res = all_results.get("shared_up", {})
    if shared_res:
        log(f"  Libraries with significant hits (FDR<0.05):")
        for lib_key, df in shared_res.items():
            padj_col = next((c for c in ["Adj. P-value","Adjusted P-value","FDR"]
                             if c in df.columns), None)
            if padj_col:
                n_sig = (df[padj_col] < 0.05).sum()
                log(f"    {lib_key}: {n_sig} significant terms")
                if n_sig > 0:
                    term_col = "Term" if "Term" in df.columns else df.columns[0]
                    top3 = df[df[padj_col]<0.05].nsmallest(3, padj_col)[term_col].tolist()
                    for t in top3:
                        log(f"      • {t}")

    # ── F: Cross-library summary table ────────────────────────────────────────
    log("\n── F: Building cross-library summary ───────────────────────────")
    summary_rows = []
    for label, res in all_results.items():
        for lib_key, df in res.items():
            if df is None or df.empty:
                continue
            padj_col = next((c for c in ["Adj. P-value","Adjusted P-value","FDR"]
                             if c in df.columns), None)
            if padj_col is None:
                continue
            n_sig  = (df[padj_col] < 0.05).sum()
            n_tot  = len(df)
            top1   = df.nsmallest(1, padj_col)
            term_col = "Term" if "Term" in df.columns else df.columns[0]
            top_term = top1[term_col].values[0] if not top1.empty else ""
            top_padj = top1[padj_col].values[0]  if not top1.empty else np.nan
            summary_rows.append({
                "gene_list" : label,
                "library"   : lib_key,
                "n_sig"     : n_sig,
                "n_total"   : n_tot,
                "top_term"  : top_term,
                "top_padj"  : top_padj,
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(RES_DIR / "enrichment_summary.csv", index=False)
    log(f"  Summary table → enrichment_summary.csv")

    # ── Summary ────────────────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("STEP 05 COMPLETE")
    log(f"  Enrichment results : {RES_DIR}")
    log(f"  Main figures       : {FIG_MAIN}")
    log(f"  Supplementary figs : {FIG_SUPP}")
    log(f"  Summary table      : enrichment_summary.csv")
    log("\nNext: run step06_network_analysis.py  (or review results)")
    log("=" * 60)


if __name__ == "__main__":
    main()
