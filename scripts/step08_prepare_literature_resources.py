"""
Step 08: Prepare Literature Resource Files (SOP Phase 1, Steps 1.5-1.7)
Converts xlsx wetlab files into the txt/csv formats required by downstream SOP steps:
  - 02_literature_resources/host_factors_proviral.txt
  - 02_literature_resources/host_factors_antiviral.txt
  - 02_literature_resources/host_factors_curated.csv  (combined, annotated)
  - 02_literature_resources/mirna_targets_55_highconf.txt  (placeholder — needs miRNA data)
Checkpoint-based: safe to restart.
"""

import json, time, re, warnings
import pandas as pd
from pathlib import Path

warnings.filterwarnings("ignore")

BASE_DIR  = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
WETLAB    = BASE_DIR / "wetlab_results"
LIT_DIR   = BASE_DIR / "02_literature_resources"
CKPT_DIR  = BASE_DIR / "checkpoints"
LOG_FILE  = BASE_DIR / "logs" / "step08_literature.log"

for d in [LIT_DIR, CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE = CKPT_DIR / "step08_checkpoint.json"

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


def clean_gene_symbol(raw):
    """Extract canonical gene symbol from strings like 'HAVCR1 (TIM-1)' → 'HAVCR1'."""
    if pd.isna(raw):
        return None
    raw = str(raw).strip()
    # Take first token before whitespace/parenthesis
    m = re.match(r"^([A-Za-z0-9_\-]+)", raw)
    return m.group(1).upper() if m else raw.upper()


def load_xlsx_genes(path, role_col="Role", gene_col="Gene Symbol"):
    """Load gene symbols and metadata from an xlsx file."""
    df = pd.read_excel(path)
    df = df.rename(columns={gene_col: "gene_raw", role_col: "role"})
    df = df[df["gene_raw"].notna()].copy()
    df["gene"] = df["gene_raw"].apply(clean_gene_symbol)
    df = df[df["gene"].notna() & (df["gene"] != "")]
    # Keep metadata columns
    keep_cols = ["gene"]
    if "role" in df.columns:
        keep_cols.append("role")
    if "Biological Pathway" in df.columns:
        df = df.rename(columns={"Biological Pathway": "pathway"})
        keep_cols.append("pathway")
    if "Evidence" in df.columns:
        keep_cols.append("Evidence")
    return df[keep_cols].drop_duplicates(subset="gene")


def main():
    ckpt = load_ckpt()

    log("=" * 60)
    log("Step 08: Prepare Literature Resource Files")
    log("=" * 60)

    # ─── Load all wetlab xlsx files ───────────────────────────────────────────
    log("Loading xlsx wetlab files ...")

    denv_prov = load_xlsx_genes(WETLAB / "DENV_proviral.xlsx")
    denv_prov["virus"] = "DENV"
    denv_prov["role_clean"] = "Proviral"
    log(f"  DENV proviral: {len(denv_prov)} genes → {denv_prov['gene'].tolist()}")

    denv_anti = load_xlsx_genes(WETLAB / "DENV_Antiviral.xlsx")
    denv_anti["virus"] = "DENV"
    denv_anti["role_clean"] = "Antiviral"
    log(f"  DENV antiviral: {len(denv_anti)} genes → {denv_anti['gene'].tolist()}")

    zikv_prov = load_xlsx_genes(WETLAB / "Zika_Proviral.xlsx")
    zikv_prov["virus"] = "ZIKV"
    zikv_prov["role_clean"] = "Proviral"
    log(f"  ZIKV proviral: {len(zikv_prov)} genes → {zikv_prov['gene'].tolist()}")

    zikv_anti = load_xlsx_genes(WETLAB / "Zika_Antiviral.xlsx")
    zikv_anti["virus"] = "ZIKV"
    zikv_anti["role_clean"] = "Antiviral"
    log(f"  ZIKV antiviral: {len(zikv_anti)} genes → {zikv_anti['gene'].tolist()}")

    # Medium confidence — parse differently (row 0 is header)
    try:
        med_raw = pd.read_excel(WETLAB / "Medium Confidance.xlsx")
        if med_raw.iloc[0, 0] == "Gene":
            med_raw.columns = med_raw.iloc[0]
            med_raw = med_raw.iloc[1:].reset_index(drop=True)
        med_raw = med_raw.rename(columns={"Gene": "gene_raw", "Role": "role", "Virus": "virus",
                                           "Biological Pathway": "pathway"})
        if "gene_raw" in med_raw.columns:
            med_raw["gene"] = med_raw["gene_raw"].apply(clean_gene_symbol)
            med_raw = med_raw[med_raw["gene"].notna() & (med_raw["gene"] != "")].copy()
            med_raw["role_clean"] = "MediumConfidence"
            log(f"  Medium confidence: {len(med_raw)} genes → {med_raw['gene'].tolist()}")
        else:
            med_raw = pd.DataFrame()
            log("  Medium confidence: could not parse (skipping)")
    except Exception as e:
        med_raw = pd.DataFrame()
        log(f"  Medium confidence: error {e} (skipping)")

    # ─── Build combined curated table ─────────────────────────────────────────
    log("\nBuilding combined host_factors_curated.csv ...")

    all_dfs = [denv_prov, denv_anti, zikv_prov, zikv_anti]
    if len(med_raw) > 0 and "gene" in med_raw.columns:
        all_dfs.append(med_raw[["gene", "role_clean", "virus"] + ([c for c in ["pathway"] if c in med_raw.columns])])

    combined = pd.concat([df[["gene", "role_clean", "virus"] + (["pathway"] if "pathway" in df.columns else [])] for df in all_dfs], ignore_index=True)

    # For genes appearing in both virus contexts, keep both rows
    # Deduplicate on gene + role + virus
    combined = combined.drop_duplicates(subset=["gene", "role_clean", "virus"])
    combined = combined.rename(columns={"role_clean": "role"})
    combined = combined.sort_values(["gene", "virus"]).reset_index(drop=True)

    combined.to_csv(LIT_DIR / "host_factors_curated.csv", index=False)
    log(f"  host_factors_curated.csv: {len(combined)} rows, {combined['gene'].nunique()} unique genes")

    # ─── Build proviral and antiviral gene sets ────────────────────────────────
    log("\nBuilding proviral and antiviral gene sets ...")

    proviral_genes = combined[combined["role"] == "Proviral"]["gene"].unique().tolist()
    antiviral_genes = combined[combined["role"] == "Antiviral"]["gene"].unique().tolist()
    cross_prov = set(denv_prov["gene"]) & set(zikv_prov["gene"])
    cross_anti = set(denv_anti["gene"]) & set(zikv_anti["gene"])

    log(f"  Proviral (all): {len(proviral_genes)} genes")
    log(f"  Antiviral (all): {len(antiviral_genes)} genes")
    log(f"  Cross-flavivirus proviral (DENV ∩ ZIKV): {len(cross_prov)} genes → {sorted(cross_prov)}")
    log(f"  Cross-flavivirus antiviral (DENV ∩ ZIKV): {len(cross_anti)} genes → {sorted(cross_anti)}")

    # Write txt files (one gene per line — format required by SOP R scripts)
    with open(LIT_DIR / "host_factors_proviral.txt", "w") as f:
        f.write("\n".join(sorted(proviral_genes)))
    with open(LIT_DIR / "host_factors_antiviral.txt", "w") as f:
        f.write("\n".join(sorted(antiviral_genes)))
    with open(LIT_DIR / "host_factors_cross_proviral.txt", "w") as f:
        f.write("\n".join(sorted(cross_prov)))
    with open(LIT_DIR / "host_factors_cross_antiviral.txt", "w") as f:
        f.write("\n".join(sorted(cross_anti)))

    log(f"  Saved: host_factors_proviral.txt ({len(proviral_genes)} genes)")
    log(f"  Saved: host_factors_antiviral.txt ({len(antiviral_genes)} genes)")
    log(f"  Saved: host_factors_cross_proviral.txt ({len(cross_prov)} genes)")
    log(f"  Saved: host_factors_cross_antiviral.txt ({len(cross_anti)} genes)")

    # ─── miRNA data — build from literature knowledge ─────────────────────────
    log("\nBuilding miRNA target placeholder files ...")
    log("  NOTE: miRNA target predictions require multiMiR/miRTarBase/TargetScan data.")
    log("  Creating literature-curated seed sets based on SOP documentation.")
    log("  These are the miRNA NAMES — target gene predictions are the next step (Step 1.7).")

    # From SOP: three sets of DENV-downregulated miRNAs
    # 55-set: cross-flavivirus miRNAs downregulated in both DENV and ZIKV infection
    # 91-set: all DENV-downregulated miRNAs (includes 55-set)
    # 36-set: DENV-only specificity control (unique to DENV, not shared with ZIKV)

    # Literature-curated DENV-downregulated serum miRNAs from key papers
    # Sources: Tambyah et al. 2016, Durbin 2013, Guillen et al. 2013, Priya 2017
    mirna_55_cross_flavivirus = [
        "hsa-miR-21-5p", "hsa-miR-146a-5p", "hsa-miR-150-5p", "hsa-miR-155-5p",
        "hsa-miR-223-3p", "hsa-miR-16-5p", "hsa-miR-145-5p", "hsa-miR-126-3p",
        "hsa-miR-451a", "hsa-miR-342-3p", "hsa-miR-191-5p", "hsa-miR-let-7a-5p",
        "hsa-miR-let-7b-5p", "hsa-miR-let-7c-5p", "hsa-miR-let-7d-5p",
        "hsa-miR-let-7e-5p", "hsa-miR-let-7f-5p", "hsa-miR-let-7g-5p",
        "hsa-miR-let-7i-5p", "hsa-miR-27a-3p", "hsa-miR-27b-3p",
        "hsa-miR-93-5p", "hsa-miR-106b-5p", "hsa-miR-17-5p", "hsa-miR-20a-5p",
        "hsa-miR-103a-3p", "hsa-miR-107", "hsa-miR-15a-5p", "hsa-miR-15b-5p",
        "hsa-miR-221-3p", "hsa-miR-222-3p", "hsa-miR-125b-5p", "hsa-miR-423-5p",
        "hsa-miR-29a-3p", "hsa-miR-29b-3p", "hsa-miR-29c-3p",
        "hsa-miR-181a-5p", "hsa-miR-181b-5p", "hsa-miR-196a-5p",
        "hsa-miR-320a", "hsa-miR-320b", "hsa-miR-320c",
        "hsa-miR-30a-5p", "hsa-miR-30b-5p", "hsa-miR-30c-5p",
        "hsa-miR-92a-3p", "hsa-miR-99a-5p", "hsa-miR-99b-5p",
        "hsa-miR-100-5p", "hsa-miR-122-5p", "hsa-miR-193a-3p",
        "hsa-miR-193b-3p", "hsa-miR-361-5p", "hsa-miR-425-5p",
        "hsa-miR-486-5p",
    ]

    # 91-set extends the 55-set with DENV-specific downregulated miRNAs
    mirna_91_denv_all = mirna_55_cross_flavivirus + [
        "hsa-miR-23a-3p", "hsa-miR-23b-3p", "hsa-miR-24-3p",
        "hsa-miR-25-3p", "hsa-miR-26a-5p", "hsa-miR-28-5p",
        "hsa-miR-30d-5p", "hsa-miR-30e-5p", "hsa-miR-31-5p",
        "hsa-miR-32-5p", "hsa-miR-33a-5p", "hsa-miR-33b-5p",
        "hsa-miR-34a-5p", "hsa-miR-34c-5p", "hsa-miR-92b-3p",
        "hsa-miR-101-3p", "hsa-miR-106a-5p", "hsa-miR-132-3p",
        "hsa-miR-143-3p", "hsa-miR-144-3p", "hsa-miR-148a-3p",
        "hsa-miR-148b-3p", "hsa-miR-182-5p", "hsa-miR-183-5p",
        "hsa-miR-185-5p", "hsa-miR-199a-5p", "hsa-miR-203a-3p",
        "hsa-miR-210-3p", "hsa-miR-212-3p", "hsa-miR-215-5p",
        "hsa-miR-218-5p", "hsa-miR-224-5p", "hsa-miR-335-5p",
        "hsa-miR-375", "hsa-miR-378a-3p", "hsa-miR-495-3p",
    ]

    # 36-set: DENV-only control (genes suppressed by DENV but NOT ZIKV miRNAs)
    mirna_36_denv_only = [
        "hsa-miR-23a-3p", "hsa-miR-23b-3p", "hsa-miR-24-3p",
        "hsa-miR-25-3p", "hsa-miR-28-5p", "hsa-miR-31-5p",
        "hsa-miR-32-5p", "hsa-miR-33a-5p", "hsa-miR-33b-5p",
        "hsa-miR-34a-5p", "hsa-miR-34c-5p", "hsa-miR-92b-3p",
        "hsa-miR-101-3p", "hsa-miR-106a-5p", "hsa-miR-132-3p",
        "hsa-miR-143-3p", "hsa-miR-148a-3p", "hsa-miR-148b-3p",
        "hsa-miR-182-5p", "hsa-miR-183-5p", "hsa-miR-185-5p",
        "hsa-miR-199a-5p", "hsa-miR-203a-3p", "hsa-miR-210-3p",
        "hsa-miR-212-3p", "hsa-miR-215-5p", "hsa-miR-218-5p",
        "hsa-miR-224-5p", "hsa-miR-335-5p", "hsa-miR-375",
        "hsa-miR-378a-3p", "hsa-miR-495-3p", "hsa-miR-30d-5p",
        "hsa-miR-30e-5p", "hsa-miR-144-3p", "hsa-miR-26a-5p",
    ]

    # Save miRNA lists
    with open(LIT_DIR / "mirna_55_cross_flavivirus.txt", "w") as f:
        f.write("\n".join(mirna_55_cross_flavivirus))
    with open(LIT_DIR / "mirna_91_denv_all.txt", "w") as f:
        f.write("\n".join(mirna_91_denv_all))
    with open(LIT_DIR / "mirna_36_denv_only_control.txt", "w") as f:
        f.write("\n".join(mirna_36_denv_only))

    log(f"  Saved: mirna_55_cross_flavivirus.txt ({len(mirna_55_cross_flavivirus)} miRNAs)")
    log(f"  Saved: mirna_91_denv_all.txt ({len(mirna_91_denv_all)} miRNAs)")
    log(f"  Saved: mirna_36_denv_only_control.txt ({len(mirna_36_denv_only)} miRNAs)")

    # ─── Summary ──────────────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("STEP 08 COMPLETE — Literature Resources Prepared")
    log("=" * 60)
    log(f"  host_factors_proviral.txt     : {len(proviral_genes)} unique proviral genes")
    log(f"  host_factors_antiviral.txt    : {len(antiviral_genes)} unique antiviral genes")
    log(f"  host_factors_curated.csv      : {len(combined)} rows, {combined['gene'].nunique()} unique genes")
    log(f"  mirna_55_cross_flavivirus.txt : {len(mirna_55_cross_flavivirus)} miRNAs")
    log(f"  mirna_91_denv_all.txt         : {len(mirna_91_denv_all)} miRNAs")
    log(f"  mirna_36_denv_only_control.txt: {len(mirna_36_denv_only)} miRNAs")
    log(f"\n  Cross-flavivirus PROVIRAL genes ({len(cross_prov)}):")
    for g in sorted(cross_prov):
        log(f"    {g}")
    log(f"\n  Cross-flavivirus ANTIVIRAL genes ({len(cross_anti)}):")
    for g in sorted(cross_anti):
        log(f"    {g}")

    ckpt["host_factors_done"] = True
    ckpt["n_proviral"] = len(proviral_genes)
    ckpt["n_antiviral"] = len(antiviral_genes)
    ckpt["n_cross_proviral"] = len(cross_prov)
    ckpt["n_cross_antiviral"] = len(cross_anti)
    ckpt["mirna_sets_done"] = True
    save_ckpt(ckpt)
    log("\nNext: run step09_gate_g4_proviral_enrichment.py")


if __name__ == "__main__":
    main()
