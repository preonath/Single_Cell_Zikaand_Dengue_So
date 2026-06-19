"""
Step 04: Gene annotation — Ensembl ID → HGNC symbol
- Annotate all DEG tables
- Identify the 15 shared genes by name
- Flag known ISGs, interferon pathway genes
- Generate annotated DEG summary tables
Checkpoint-based: safe to restart.
"""

import json, time, warnings
import numpy as np
import pandas as pd
import mygene
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path

warnings.filterwarnings("ignore")

# ─── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR  = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
DEG_DIR   = BASE_DIR / "01_processed_data" / "deg_tables"
RES_DIR   = BASE_DIR / "03_results" / "phase3_shared_degs"
FIG_MAIN  = BASE_DIR / "04_figures" / "main"
FIG_SUPP  = BASE_DIR / "04_figures" / "supplementary"
CKPT_DIR  = BASE_DIR / "checkpoints"
LOG_FILE  = BASE_DIR / "logs" / "step04_annotation.log"

for d in [CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE = CKPT_DIR / "step04_checkpoint.json"

# ─── Known gene sets for biological validation ────────────────────────────────
ISGS = {
    "ISG15","MX1","MX2","OAS1","OAS2","OAS3","OASL","IFIT1","IFIT2","IFIT3",
    "IFITM1","IFITM2","IFITM3","IFI6","IFI27","IFI44","IFI44L","RSAD2",
    "HERC5","HERC6","USP18","TRIM25","TRIM56","RIG-I","DDX58","IFIH1",
    "IRF7","IRF9","STAT1","STAT2","BST2","APOBEC3G","CXCL10","CCL5"
}

PROVIRAL_KNOWN = {
    "BCL2L1","BCL2","MCL1","BIRC2","BIRC3",          # anti-apoptotic
    "FASN","ACACA","HMGCR","LDLR","SREBF1","PPARA",   # lipid metabolism
    "ATG5","ATG7","ATG12","BECN1","ULK1",              # autophagy
    "HSP90AA1","HSP90AB1","HSPA8","HSPA1A",            # chaperones
    "EIF4A1","EIF4E","EIF4G1",                          # translation
    "RACK1","GNB2L1",                                   # scaffold
    "DZIP1","NPC1","NPC2",                              # cholesterol
    "SEC61A1","SEC61B","TRAM1",                         # ER translocation
    "RAB5A","RAB7A","RAB18",                            # endosomal trafficking
    "VAPA","VAPB"                                       # ER-membrane contact
}

ANTIVIRAL_KNOWN = ISGS | {
    "MAVS","STING1","CGAS","TBK1","IRF3","IRF5",
    "TRIM56","TRIM27","RNF135",
    "PKR","EIF2AK2","OAS1","RNASEL",
    "IFNA1","IFNB1","IFNG","IL6","TNF","CXCL10"
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


# ─── Step A: Query MyGene.info ─────────────────────────────────────────────────
def annotate_ensembl_ids(ensembl_ids: list) -> pd.DataFrame:
    log(f"Querying MyGene.info for {len(ensembl_ids)} Ensembl IDs ...")
    mg  = mygene.MyGeneInfo()

    # Query in batches of 1000
    results = []
    batch_size = 1000
    for i in range(0, len(ensembl_ids), batch_size):
        batch = ensembl_ids[i:i+batch_size]
        res   = mg.querymany(
            batch,
            scopes  = "ensembl.gene",
            fields  = "symbol,name,entrezgene,ensembl.gene,type_of_gene",
            species = "human",
            as_dataframe = True,
            returnall    = False
        )
        results.append(res)
        log(f"  Batch {i//batch_size + 1}: {len(batch)} queried")
        time.sleep(0.5)

    df = pd.concat(results)
    df = df.reset_index().rename(columns={"query":"gene_id"})

    # Keep best hit per Ensembl ID (highest score)
    if "symbol" in df.columns:
        df = df.sort_values("_score", ascending=False)
        df = df.drop_duplicates(subset="gene_id", keep="first")

    log(f"  Mapped: {df['symbol'].notna().sum()} / {len(ensembl_ids)} IDs")
    return df[["gene_id","symbol","name","entrezgene","type_of_gene"]].copy()


# ─── Step B: Annotate DEG tables ──────────────────────────────────────────────
def annotate_deg_table(deg_df: pd.DataFrame,
                       annot_df: pd.DataFrame,
                       label: str) -> pd.DataFrame:
    merged = deg_df.merge(
        annot_df[["gene_id","symbol","name","entrezgene"]],
        on="gene_id", how="left"
    )
    mapped = merged["symbol"].notna().sum()
    log(f"  {label}: {mapped}/{len(merged)} genes have symbols")

    # Biological category flags
    merged["is_ISG"]      = merged["symbol"].isin(ISGS)
    merged["is_proviral"] = merged["symbol"].isin(PROVIRAL_KNOWN)
    merged["is_antiviral"]= merged["symbol"].isin(ANTIVIRAL_KNOWN)

    # Significance flag
    merged["significant"] = (
        (merged["padj"] < 0.05) &
        (merged["log2FoldChange"].abs() >= 1.0)
    )
    return merged


# ─── Step C: Print shared gene details ────────────────────────────────────────
def report_shared_genes(shared_df: pd.DataFrame):
    log("\n" + "=" * 60)
    log("THE 15 SHARED UPREGULATED GENES (DENV & ZIKV vs Control):")
    log("=" * 60)

    cols = ["symbol","name","log2FC_DENV","padj_DENV","log2FC_ZIKV","padj_ZIKV",
            "is_ISG","is_proviral","is_antiviral"]
    available = [c for c in cols if c in shared_df.columns]

    for _, row in shared_df.sort_values("log2FC_DENV", ascending=False).iterrows():
        sym  = row.get("symbol","?")
        name = str(row.get("name",""))[:60]
        fc_d = row.get("log2FC_DENV", float("nan"))
        fc_z = row.get("log2FC_ZIKV", float("nan"))
        flags = []
        if row.get("is_ISG"):      flags.append("ISG")
        if row.get("is_proviral"): flags.append("PROVIRAL")
        if row.get("is_antiviral"):flags.append("ANTIVIRAL")
        flag_str = f"  [{', '.join(flags)}]" if flags else ""
        log(f"  {sym:<12} FC_DENV={fc_d:+.2f}  FC_ZIKV={fc_z:+.2f}  {name}{flag_str}")

    n_isg  = shared_df["is_ISG"].sum()
    n_prov = shared_df["is_proviral"].sum()
    n_anti = shared_df["is_antiviral"].sum()
    log(f"\n  ISGs in shared genes   : {n_isg}")
    log(f"  Proviral in shared     : {n_prov}")
    log(f"  Antiviral in shared    : {n_anti}")


# ─── Step D: Annotated volcano with gene names ────────────────────────────────
def make_annotated_volcano(deg_df: pd.DataFrame, title: str,
                           color_up: str, color_dn: str,
                           shared_genes: list, out_path: Path):
    df = deg_df.dropna(subset=["padj","log2FoldChange"]).copy()
    df["-log10padj"] = -np.log10(df["padj"].clip(lower=1e-300))

    sig_up   = (df["padj"]<0.05) & (df["log2FoldChange"]>=1.0)
    sig_down = (df["padj"]<0.05) & (df["log2FoldChange"]<=-1.0)
    is_shared= df["gene_id"].isin(shared_genes)
    ns       = ~(sig_up | sig_down)

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(df.loc[ns,  "log2FoldChange"], df.loc[ns,  "-log10padj"],
               c="lightgrey", s=5, alpha=0.4, rasterized=True, label="NS")
    ax.scatter(df.loc[sig_up & ~is_shared,  "log2FoldChange"],
               df.loc[sig_up & ~is_shared,  "-log10padj"],
               c=color_up,   s=10, alpha=0.8, rasterized=True,
               label=f"Up ({sig_up.sum()})")
    ax.scatter(df.loc[sig_down & ~is_shared, "log2FoldChange"],
               df.loc[sig_down & ~is_shared, "-log10padj"],
               c=color_dn,   s=10, alpha=0.8, rasterized=True,
               label=f"Down ({sig_down.sum()})")
    ax.scatter(df.loc[is_shared, "log2FoldChange"],
               df.loc[is_shared, "-log10padj"],
               c="#FF7F00", s=60, alpha=1.0, zorder=5,
               edgecolors="black", linewidths=0.5,
               label=f"Shared ({is_shared.sum()})")

    # Label shared genes by symbol
    for _, row in df[is_shared].iterrows():
        sym = row.get("symbol","")
        if pd.notna(sym) and sym:
            ax.annotate(sym,
                xy=(row["log2FoldChange"], row["-log10padj"]),
                xytext=(5, 3), textcoords="offset points",
                fontsize=7, fontweight="bold",
                arrowprops=dict(arrowstyle="-", lw=0.5))

    ax.axvline(x= 1.0, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(x=-1.0, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(y=-np.log10(0.05), color="black", linestyle=":", linewidth=0.8)
    ax.set_xlabel("log₂ Fold Change", fontsize=12)
    ax.set_ylabel("-log₁₀(padj)", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.legend(framealpha=0.8, fontsize=9)
    ax.set_xlim(-10, 10)
    fig.tight_layout()
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close()
    log(f"  Annotated volcano → {out_path.name}")


# ─── Step E: Shared gene heatmap ──────────────────────────────────────────────
def make_shared_heatmap(shared_df: pd.DataFrame):
    log("Generating shared gene heatmap ...")

    plot_df = shared_df[["symbol","log2FC_DENV","log2FC_ZIKV"]].copy()
    plot_df = plot_df.dropna(subset=["symbol"])
    plot_df = plot_df.set_index("symbol")
    plot_df.columns = ["DENV vs Control", "ZIKV vs Control"]
    plot_df = plot_df.sort_values("DENV vs Control", ascending=False)

    fig, ax = plt.subplots(figsize=(5, max(4, len(plot_df)*0.4)))
    sns.heatmap(
        plot_df, ax=ax,
        cmap="RdBu_r", center=0,
        annot=True, fmt=".2f", annot_kws={"size":9},
        linewidths=0.5, linecolor="white",
        cbar_kws={"label":"log₂ Fold Change"}
    )
    ax.set_title("Shared upregulated genes\n(DENV & ZIKV vs Control)",
                 fontsize=12, fontweight="bold")
    ax.set_xlabel(""); ax.set_ylabel("Gene symbol", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG_MAIN / "Figure_SharedGenes_Heatmap.pdf", bbox_inches="tight")
    fig.savefig(FIG_MAIN / "Figure_SharedGenes_Heatmap.png", bbox_inches="tight", dpi=150)
    plt.close()
    log("  Shared gene heatmap saved → FIG_MAIN")


# ─── Step F: Top DEG tables for paper ─────────────────────────────────────────
def make_top_deg_tables(denv_ann: pd.DataFrame, zikv_ann: pd.DataFrame):
    for df, label in [(denv_ann,"DENV"), (zikv_ann,"ZIKV")]:
        sig = df[df["significant"]].copy()
        sig = sig.sort_values("padj")
        top = sig[["symbol","name","log2FoldChange","baseMean","stat",
                    "pvalue","padj","is_ISG","is_proviral","is_antiviral"]]
        top = top.rename(columns={"log2FoldChange":"log2FC"})
        out = DEG_DIR / f"DEGs_{label}_vs_Control_annotated.csv"
        top.to_csv(out, index=False)
        log(f"  Annotated DEG table saved → {out.name} ({len(top)} genes)")


# ─── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    ckpt = load_ckpt()
    log("=" * 60)
    log("Step 04: Gene Annotation (Ensembl → Symbol)")
    log("=" * 60)

    # Load DEG tables
    denv_res = pd.read_csv(DEG_DIR / "DEGs_DENV_vs_Control.csv")
    zikv_res = pd.read_csv(DEG_DIR / "DEGs_ZIKV_vs_Control.csv")
    shared_up_ids = pd.read_csv(RES_DIR / "shared_DEGs_shared_up.csv")["gene_id"].tolist()

    log(f"Loaded DENV DEGs: {len(denv_res)}, ZIKV DEGs: {len(zikv_res)}")
    log(f"Shared upregulated IDs: {len(shared_up_ids)}")

    # ── A: Query annotation ────────────────────────────────────────────────────
    annot_path = BASE_DIR / "02_literature_resources" / "ensembl_annotation.csv"
    if ckpt.get("annotation_done") and annot_path.exists():
        log("✓ Annotation already done — loading ...")
        annot_df = pd.read_csv(annot_path)
        log(f"  {annot_df['symbol'].notna().sum()}/{len(annot_df)} IDs mapped")
    else:
        all_ids  = list(set(denv_res["gene_id"].tolist() +
                            zikv_res["gene_id"].tolist()))
        annot_df = annotate_ensembl_ids(all_ids)
        annot_path.parent.mkdir(parents=True, exist_ok=True)
        annot_df.to_csv(annot_path, index=False)
        log(f"✓ Annotation saved → {annot_path.name}")
        ckpt["annotation_done"] = True
        save_ckpt(ckpt)

    # ── B: Annotate DEG tables ─────────────────────────────────────────────────
    log("\nAnnotating DEG tables ...")
    denv_ann = annotate_deg_table(denv_res, annot_df, "DENV")
    zikv_ann = annotate_deg_table(zikv_res, annot_df, "ZIKV")

    denv_ann.to_csv(DEG_DIR / "DEGs_DENV_vs_Control_annotated_full.csv", index=False)
    zikv_ann.to_csv(DEG_DIR / "DEGs_ZIKV_vs_Control_annotated_full.csv", index=False)

    # ── C: Build shared gene detail table ─────────────────────────────────────
    log("\nBuilding shared gene table ...")
    denv_sig = denv_ann[["gene_id","symbol","name",
                          "log2FoldChange","padj",
                          "is_ISG","is_proviral","is_antiviral"]].rename(
        columns={"log2FoldChange":"log2FC_DENV","padj":"padj_DENV"})
    zikv_sig = zikv_ann[["gene_id","log2FoldChange","padj"]].rename(
        columns={"log2FoldChange":"log2FC_ZIKV","padj":"padj_ZIKV"})

    shared_df = denv_sig[denv_sig["gene_id"].isin(shared_up_ids)].merge(
        zikv_sig, on="gene_id", how="left")

    shared_df.to_csv(RES_DIR / "shared_DEGs_annotated.csv", index=False)

    # ── D: Report shared genes ─────────────────────────────────────────────────
    report_shared_genes(shared_df)

    # ── E: Annotated volcano plots ─────────────────────────────────────────────
    if not ckpt.get("annotated_volcano_done"):
        denv_ann_merged = denv_ann.copy()
        zikv_ann_merged = zikv_ann.copy()

        make_annotated_volcano(
            denv_ann_merged, "DENV vs Control (shared genes highlighted)",
            "#E41A1C","#4575B4", shared_up_ids,
            FIG_SUPP / "volcano_DENV_annotated"
        )
        make_annotated_volcano(
            zikv_ann_merged, "ZIKV vs Control (shared genes highlighted)",
            "#E41A1C","#4575B4", shared_up_ids,
            FIG_SUPP / "volcano_ZIKV_annotated"
        )
        ckpt["annotated_volcano_done"] = True
        save_ckpt(ckpt)
    else:
        log("✓ Annotated volcanos already saved")

    # ── F: Shared gene heatmap ─────────────────────────────────────────────────
    if not ckpt.get("heatmap_done"):
        make_shared_heatmap(shared_df)
        ckpt["heatmap_done"] = True
        save_ckpt(ckpt)
    else:
        log("✓ Heatmap already saved")

    # ── G: Top DEG tables ─────────────────────────────────────────────────────
    make_top_deg_tables(denv_ann, zikv_ann)

    # ── H: ISG enrichment check ───────────────────────────────────────────────
    log("\nISG enrichment in shared upregulated genes:")
    denv_sig_genes = set(denv_ann[denv_ann["significant"] &
                                   (denv_ann["log2FoldChange"]>0)]["symbol"].dropna())
    zikv_sig_genes = set(zikv_ann[zikv_ann["significant"] &
                                   (zikv_ann["log2FoldChange"]>0)]["symbol"].dropna())
    shared_sym     = set(shared_df["symbol"].dropna())

    isg_in_denv   = denv_sig_genes & ISGS
    isg_in_zikv   = zikv_sig_genes & ISGS
    isg_in_shared = shared_sym & ISGS

    log(f"  ISGs upregulated in DENV : {len(isg_in_denv)} — {sorted(isg_in_denv)}")
    log(f"  ISGs upregulated in ZIKV : {len(isg_in_zikv)} — {sorted(isg_in_zikv)}")
    log(f"  ISGs in shared genes     : {len(isg_in_shared)} — {sorted(isg_in_shared)}")

    # ── Summary ────────────────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("STEP 04 COMPLETE")
    log(f"  Annotation file   : {annot_path.name}")
    log(f"  Annotated DEG DENV: DEGs_DENV_vs_Control_annotated.csv")
    log(f"  Annotated DEG ZIKV: DEGs_ZIKV_vs_Control_annotated.csv")
    log(f"  Shared genes table: shared_DEGs_annotated.csv")
    log(f"  Heatmap + volcanos: {FIG_MAIN} / {FIG_SUPP}")
    log("\nNext: run step05_pathway_enrichment.py")
    log("=" * 60)


if __name__ == "__main__":
    main()
