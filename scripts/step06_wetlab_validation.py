"""
Step 06: Wetlab Validation
- Cross-reference computational DEGs against literature-curated
  antiviral / proviral gene lists from wetlab_results/
- Report overlaps, mismatches, and novel candidates
- Generate validation figures and summary tables
"""

import json, time, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from matplotlib_venn import venn2

warnings.filterwarnings("ignore")

# ─── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR    = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
DEG_DIR     = BASE_DIR / "01_processed_data" / "deg_tables"
SHARED_DIR  = BASE_DIR / "03_results" / "phase3_shared_degs"
WETLAB_DIR  = BASE_DIR / "wetlab_results"
RES_DIR     = BASE_DIR / "03_results" / "phase5_validation"
FIG_MAIN    = BASE_DIR / "04_figures" / "main"
FIG_SUPP    = BASE_DIR / "04_figures" / "supplementary"
CKPT_DIR    = BASE_DIR / "checkpoints"
LOG_FILE    = BASE_DIR / "logs" / "step06_validation.log"

for d in [RES_DIR, FIG_MAIN, FIG_SUPP, CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE = CKPT_DIR / "step06_checkpoint.json"

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

def gene_set(df: pd.DataFrame) -> set:
    col = next((c for c in ["Gene Symbol","gene_symbol","symbol","Gene"]
                if c in df.columns), df.columns[1])
    return set(df[col].dropna().str.strip().str.upper().tolist())


# ─── Load wetlab gene lists ────────────────────────────────────────────────────
def load_wetlab() -> dict:
    files = {
        "DENV_antiviral" : WETLAB_DIR / "DENV_Antiviral.xlsx",
        "DENV_proviral"  : WETLAB_DIR / "DENV_proviral.xlsx",
        "ZIKV_antiviral" : WETLAB_DIR / "Zika_Antiviral.xlsx",
        "ZIKV_proviral"  : WETLAB_DIR / "Zika_Proviral.xlsx",
        "medium_conf"    : WETLAB_DIR / "Medium Confidance.xlsx",
    }
    wl = {}
    for name, fp in files.items():
        df = pd.read_excel(fp)
        if name == "medium_conf":
            # First row is header
            df.columns = df.iloc[0]
            df = df.iloc[1:].reset_index(drop=True)
        wl[name] = df
        gs = gene_set(df)
        log(f"  Loaded {name}: {len(gs)} genes")
    return wl


# ─── Load computational results ───────────────────────────────────────────────
def load_comp() -> dict:
    fc_col = "log2FC"

    denv_ann = pd.read_csv(DEG_DIR / "DEGs_DENV_vs_Control_annotated.csv")
    zikv_ann = pd.read_csv(DEG_DIR / "DEGs_ZIKV_vs_Control_annotated.csv")
    shared   = pd.read_csv(SHARED_DIR / "shared_DEGs_annotated.csv")

    denv_ann["symbol_up"] = denv_ann["symbol"].str.upper()
    zikv_ann["symbol_up"] = zikv_ann["symbol"].str.upper()

    denv_up  = set(denv_ann[denv_ann[fc_col]>0]["symbol_up"].dropna())
    denv_dn  = set(denv_ann[denv_ann[fc_col]<0]["symbol_up"].dropna())
    zikv_up  = set(zikv_ann[zikv_ann[fc_col]>0]["symbol_up"].dropna())
    zikv_dn  = set(zikv_ann[zikv_ann[fc_col]<0]["symbol_up"].dropna())
    shared_s = set(shared["symbol"].dropna().str.upper())

    log(f"  Computational: DENV up={len(denv_up)} dn={len(denv_dn)}")
    log(f"  Computational: ZIKV up={len(zikv_up)} dn={len(zikv_dn)}")
    log(f"  Shared upregulated: {len(shared_s)}")

    return {
        "denv_up": denv_up, "denv_dn": denv_dn,
        "zikv_up": zikv_up, "zikv_dn": zikv_dn,
        "shared" : shared_s,
        "denv_ann": denv_ann, "zikv_ann": zikv_ann,
        "shared_df": shared,
    }


# ─── Cross-reference one pair ─────────────────────────────────────────────────
def cross_ref(comp_set: set, wetlab_df: pd.DataFrame,
              comp_label: str, wetlab_label: str) -> pd.DataFrame:
    wl_genes = gene_set(wetlab_df)
    overlap  = comp_set & wl_genes
    comp_only= comp_set - wl_genes
    wl_only  = wl_genes - comp_set

    log(f"  {comp_label} vs {wetlab_label}:")
    log(f"    wetlab={len(wl_genes)}  comp={len(comp_set)}  overlap={len(overlap)}")
    if overlap:
        log(f"    OVERLAP genes: {sorted(overlap)}")

    # Build detail table for overlap
    sym_col = next((c for c in ["Gene Symbol","symbol","Gene"]
                    if c in wetlab_df.columns), wetlab_df.columns[1])

    rows = []
    for g in sorted(overlap):
        wl_row = wetlab_df[wetlab_df[sym_col].str.upper().str.strip() == g]
        pathway = wl_row["Biological Pathway"].values[0] if "Biological Pathway" in wl_row.columns and len(wl_row) else ""
        evidence= wl_row["Evidence"].values[0]           if "Evidence"           in wl_row.columns and len(wl_row) else ""
        rows.append({
            "gene"          : g,
            "comp_direction": comp_label,
            "wetlab_role"   : wetlab_label,
            "pathway"       : pathway,
            "evidence"      : evidence,
        })
    return pd.DataFrame(rows)


# ─── Validation matrix heatmap ────────────────────────────────────────────────
def plot_validation_matrix(results: dict, out_path: Path):
    comp_keys  = ["denv_up","denv_dn","zikv_up","zikv_dn","shared"]
    wetlab_keys= ["DENV_antiviral","DENV_proviral","ZIKV_antiviral","ZIKV_proviral"]

    labels_c = ["DENV Up","DENV Down","ZIKV Up","ZIKV Down","Shared Up"]
    labels_w = ["DENV Antiviral","DENV Proviral","ZIKV Antiviral","ZIKV Proviral"]

    mat = np.zeros((len(comp_keys), len(wetlab_keys)), dtype=int)
    for i, ck in enumerate(comp_keys):
        for j, wk in enumerate(wetlab_keys):
            df = results.get(f"{ck}_vs_{wk}")
            mat[i, j] = len(df) if df is not None else 0

    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(mat, cmap="YlOrRd", aspect="auto")

    ax.set_xticks(range(len(labels_w)))
    ax.set_xticklabels(labels_w, rotation=30, ha="right", fontsize=10)
    ax.set_yticks(range(len(labels_c)))
    ax.set_yticklabels(labels_c, fontsize=10)

    for i in range(len(labels_c)):
        for j in range(len(labels_w)):
            v = mat[i, j]
            ax.text(j, i, str(v), ha="center", va="center",
                    fontsize=12, fontweight="bold",
                    color="white" if v > mat.max() * 0.6 else "black")

    plt.colorbar(im, ax=ax, label="# overlapping genes")
    ax.set_title("Computational DEGs vs Wetlab-validated gene lists\n(number of shared genes)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close()
    log(f"  Validation matrix → {out_path.name}")


# ─── Venn: shared genes vs all wetlab ─────────────────────────────────────────
def plot_shared_venn(shared_set: set, wetlab_all: set, out_path: Path):
    try:
        fig, ax = plt.subplots(figsize=(6, 5))
        v = venn2([shared_set, wetlab_all],
                  set_labels=("Shared DEGs\n(comp)", "Wetlab\nvalidated"),
                  ax=ax)
        ax.set_title("Shared upregulated genes vs all wetlab-validated factors",
                     fontsize=11, fontweight="bold")
        fig.tight_layout()
        fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
        fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
        plt.close()
        log(f"  Venn diagram → {out_path.name}")
    except ImportError:
        log("  ! matplotlib_venn not installed — skipping Venn")


# ─── Bar: overlap per wetlab category ────────────────────────────────────────
def plot_overlap_bars(results: dict, out_path: Path):
    rows = []
    for key, df in results.items():
        if df is not None and not df.empty:
            parts = key.split("_vs_")
            rows.append({"comparison": key.replace("_vs_", " vs "),
                         "n_overlap": len(df)})
    if not rows:
        return

    df_bar = pd.DataFrame(rows).sort_values("n_overlap", ascending=True)
    fig, ax = plt.subplots(figsize=(9, max(4, len(df_bar) * 0.45)))
    colors_map = {
        "antiviral": "#1F78B4",
        "proviral" : "#E31A1C",
        "conf"     : "#FF7F00",
    }
    bar_colors = [
        colors_map["antiviral"] if "antiviral" in r["comparison"].lower() else
        colors_map["conf"]      if "medium"    in r["comparison"].lower() else
        colors_map["proviral"]
        for _, r in df_bar.iterrows()
    ]
    ax.barh(df_bar["comparison"], df_bar["n_overlap"],
            color=bar_colors, edgecolor="white")
    ax.set_xlabel("Number of overlapping genes", fontsize=11)
    ax.set_title("Computational DEGs overlapping wetlab-validated lists",
                 fontsize=12, fontweight="bold")
    patches = [mpatches.Patch(color=c, label=l)
               for l, c in [("Antiviral","#1F78B4"),("Proviral","#E31A1C")]]
    ax.legend(handles=patches, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close()
    log(f"  Overlap bar chart → {out_path.name}")


# ─── Novel candidates: in comp but NOT in any wetlab list ────────────────────
def find_novel_candidates(comp: dict, wetlab: dict) -> pd.DataFrame:
    all_wetlab = set()
    for df in wetlab.values():
        all_wetlab |= gene_set(df)

    novel_shared = comp["shared"] - all_wetlab
    novel_denv   = comp["denv_up"] - all_wetlab
    novel_zikv   = comp["zikv_up"] - all_wetlab

    log(f"\n  Novel shared upregulated (not in any wetlab list): {len(novel_shared)}")
    for g in sorted(novel_shared):
        log(f"    • {g}")

    log(f"  Novel DENV up (not in any wetlab list): {len(novel_denv)} genes")
    log(f"  Novel ZIKV up (not in any wetlab list): {len(novel_zikv)} genes")

    rows = []
    for g in sorted(novel_shared):
        rows.append({"gene": g, "source": "shared_up",
                     "note": "shared DENV+ZIKV upregulated, not in any wetlab list"})
    return pd.DataFrame(rows)


# ─── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    ckpt = load_ckpt()
    log("=" * 60)
    log("Step 06: Wetlab Validation")
    log("=" * 60)

    # ── Load data ──────────────────────────────────────────────────────────────
    log("\nLoading wetlab gene lists ...")
    wetlab = load_wetlab()

    log("\nLoading computational results ...")
    comp   = load_comp()

    # ── A: Cross-reference all pairs ──────────────────────────────────────────
    log("\n── A: Cross-referencing computational DEGs vs wetlab lists ─────")

    pairs = {
        "denv_up" : ("DENV Up", ["DENV_antiviral","DENV_proviral",
                                  "ZIKV_antiviral","ZIKV_proviral"]),
        "denv_dn" : ("DENV Down", ["DENV_antiviral","DENV_proviral"]),
        "zikv_up" : ("ZIKV Up",  ["ZIKV_antiviral","ZIKV_proviral",
                                   "DENV_antiviral","DENV_proviral"]),
        "zikv_dn" : ("ZIKV Down", ["ZIKV_antiviral","ZIKV_proviral"]),
        "shared"  : ("Shared Up", ["DENV_antiviral","DENV_proviral",
                                   "ZIKV_antiviral","ZIKV_proviral",
                                   "medium_conf"]),
    }

    all_overlaps = {}
    all_rows     = []

    for comp_key, (comp_label, wl_keys) in pairs.items():
        for wk in wl_keys:
            key = f"{comp_key}_vs_{wk}"
            df  = cross_ref(comp[comp_key], wetlab[wk], comp_label, wk)
            all_overlaps[key] = df
            if not df.empty:
                all_rows.append(df)

    # Combined overlap table
    if all_rows:
        combined = pd.concat(all_rows, ignore_index=True)
        combined.to_csv(RES_DIR / "validation_overlap_all.csv", index=False)
        log(f"\n  Combined overlap table: {len(combined)} rows → validation_overlap_all.csv")
    else:
        log("\n  No overlaps found in any category")
        combined = pd.DataFrame()

    # ── B: Novel candidates ───────────────────────────────────────────────────
    log("\n── B: Novel candidates (computational, not in wetlab) ───────────")
    novel_df = find_novel_candidates(comp, wetlab)
    novel_df.to_csv(RES_DIR / "novel_candidates_shared.csv", index=False)

    # ── C: Validation matrix heatmap ──────────────────────────────────────────
    log("\n── C: Validation matrix heatmap ────────────────────────────────")
    plot_validation_matrix(all_overlaps, FIG_MAIN / "Figure_Validation_Matrix")

    # ── D: Overlap bar chart ──────────────────────────────────────────────────
    log("\n── D: Overlap bar chart ────────────────────────────────────────")
    plot_overlap_bars(all_overlaps, FIG_SUPP / "barplot_wetlab_overlaps")

    # ── E: Venn — shared genes vs all wetlab ──────────────────────────────────
    log("\n── E: Venn diagram ─────────────────────────────────────────────")
    all_wl_genes = set()
    for df in wetlab.values():
        all_wl_genes |= gene_set(df)
    plot_shared_venn(comp["shared"], all_wl_genes,
                     FIG_MAIN / "Figure_SharedGenes_vs_Wetlab_Venn")

    # ── F: Detailed report per gene ───────────────────────────────────────────
    log("\n── F: Detailed match report ────────────────────────────────────")

    # Shared genes: annotate each with wetlab status
    shared_df = comp["shared_df"].copy()
    shared_df["symbol_up"] = shared_df["symbol"].str.upper()

    antiviral_all = gene_set(wetlab["DENV_antiviral"]) | gene_set(wetlab["ZIKV_antiviral"])
    proviral_all  = gene_set(wetlab["DENV_proviral"])  | gene_set(wetlab["ZIKV_proviral"])
    medconf_all   = gene_set(wetlab["medium_conf"])

    shared_df["wetlab_antiviral"] = shared_df["symbol_up"].isin(antiviral_all)
    shared_df["wetlab_proviral"]  = shared_df["symbol_up"].isin(proviral_all)
    shared_df["wetlab_medconf"]   = shared_df["symbol_up"].isin(medconf_all)
    shared_df["wetlab_novel"]     = ~(shared_df["wetlab_antiviral"] |
                                      shared_df["wetlab_proviral"]  |
                                      shared_df["wetlab_medconf"])

    shared_df.to_csv(RES_DIR / "shared_DEGs_wetlab_annotated.csv", index=False)

    log("\n  Shared gene wetlab status:")
    log(f"  {'Gene':<12} {'FC_DENV':>8} {'FC_ZIKV':>8}  {'Wetlab'}  ")
    log(f"  {'-'*12} {'-'*8} {'-'*8}  {'-'*30}")
    for _, row in shared_df.sort_values("log2FC_DENV", ascending=False).iterrows():
        wl_tag = ("ANTIVIRAL" if row["wetlab_antiviral"] else
                  "PROVIRAL"  if row["wetlab_proviral"]  else
                  "MED-CONF"  if row["wetlab_medconf"]   else
                  "NOVEL")
        log(f"  {str(row.get('symbol','?')):<12} {row.get('log2FC_DENV',0):+8.2f} "
            f"{row.get('log2FC_ZIKV',0):+8.2f}  {wl_tag}")

    # ── G: Summary statistics ─────────────────────────────────────────────────
    log("\n── G: Summary ──────────────────────────────────────────────────")
    n_anti = shared_df["wetlab_antiviral"].sum()
    n_prov = shared_df["wetlab_proviral"].sum()
    n_conf = shared_df["wetlab_medconf"].sum()
    n_nov  = shared_df["wetlab_novel"].sum()
    log(f"  Shared genes confirmed as antiviral   : {n_anti}")
    log(f"  Shared genes confirmed as proviral    : {n_prov}")
    log(f"  Shared genes in medium-confidence     : {n_conf}")
    log(f"  Shared genes NOVEL (not in wetlab)    : {n_nov}")

    # DENV-specific wetlab validation
    log(f"\n  DENV upregulated vs DENV antiviral wetlab:")
    denv_anti_ov = comp["denv_up"] & gene_set(wetlab["DENV_antiviral"])
    log(f"    Overlap: {len(denv_anti_ov)} genes — {sorted(denv_anti_ov)}")

    log(f"\n  DENV upregulated vs DENV proviral wetlab:")
    denv_prov_ov = comp["denv_up"] & gene_set(wetlab["DENV_proviral"])
    log(f"    Overlap: {len(denv_prov_ov)} genes — {sorted(denv_prov_ov)}")

    log(f"\n  ZIKV upregulated vs ZIKV antiviral wetlab:")
    zikv_anti_ov = comp["zikv_up"] & gene_set(wetlab["ZIKV_antiviral"])
    log(f"    Overlap: {len(zikv_anti_ov)} genes — {sorted(zikv_anti_ov)}")

    log(f"\n  ZIKV upregulated vs ZIKV proviral wetlab:")
    zikv_prov_ov = comp["zikv_up"] & gene_set(wetlab["ZIKV_proviral"])
    log(f"    Overlap: {len(zikv_prov_ov)} genes — {sorted(zikv_prov_ov)}")

    # ── Final Summary ──────────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("STEP 06 COMPLETE")
    log(f"  Overlap table     : {RES_DIR}/validation_overlap_all.csv")
    log(f"  Novel candidates  : {RES_DIR}/novel_candidates_shared.csv")
    log(f"  Shared+wetlab annot: {RES_DIR}/shared_DEGs_wetlab_annotated.csv")
    log(f"  Figures           : {FIG_MAIN}")
    log("\nNext: review SESSION_LOG.md for summary + pipeline status")
    log("=" * 60)


if __name__ == "__main__":
    main()
