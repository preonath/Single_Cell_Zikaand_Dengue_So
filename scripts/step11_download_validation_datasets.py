"""
Step 11: Download External Validation Datasets (SOP Phase 1, Steps 1.2-1.4 & Phase 8)
Downloads:
  - GSE118305 (ZIKV macrophages — validation dataset 1)
  - GSE94892  (DENV patient PBMCs — validation dataset 2)
  - GSE78711  (ZIKV neural progenitor cells — neural extension)

For each dataset, extracts expression matrix and sample metadata.
Checkpoint-based: safe to restart.
"""

import json, time, warnings, gzip, shutil, os
import numpy as np
import pandas as pd
import GEOparse
from pathlib import Path

warnings.filterwarnings("ignore")

BASE_DIR   = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
RAW_DIR    = BASE_DIR / "00_raw_data"
CKPT_DIR   = BASE_DIR / "checkpoints"
LOG_FILE   = BASE_DIR / "logs" / "step11_validation_download.log"

for d in [CKPT_DIR, LOG_FILE.parent]:
    d.mkdir(parents=True, exist_ok=True)

CKPT_FILE = CKPT_DIR / "step11_checkpoint.json"

DATASETS = {
    "GSE118305": {
        "description": "ZIKV Macrophage RNA-seq (validation 1)",
        "virus": "ZIKV",
        "cell_type": "macrophages",
        "condition_key": ["zikv", "mock", "uninfected", "infected", "control"],
    },
    "GSE94892": {
        "description": "DENV Patient PBMCs RNA-seq (validation 2)",
        "virus": "DENV",
        "cell_type": "PBMCs",
        "condition_key": ["dengue", "denv", "healthy", "control", "patient", "fever"],
    },
    "GSE78711": {
        "description": "ZIKV Neural Progenitor Cells RNA-seq (neural extension)",
        "virus": "ZIKV",
        "cell_type": "NPCs",
        "condition_key": ["zikv", "mock", "infected", "uninfected"],
    },
}

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


def download_and_parse_gse(gse_id, out_dir):
    """Download a GEO series and extract expression matrix + metadata."""
    out_dir.mkdir(parents=True, exist_ok=True)

    log(f"  Downloading {gse_id} from GEO ...")
    try:
        gse = GEOparse.get_GEO(geo=gse_id, destdir=str(out_dir), silent=True)
    except Exception as e:
        log(f"  ERROR downloading {gse_id}: {e}")
        return None, None

    log(f"  {gse_id}: {len(gse.gsms)} samples")

    # ─── Extract metadata ──────────────────────────────────────────────────
    meta_rows = []
    for gsm_id, gsm in gse.gsms.items():
        row = {
            "sample_id": gsm_id,
            "title": gsm.metadata.get("title", [""])[0],
            "source": gsm.metadata.get("source_name_ch1", [""])[0],
            "organism": gsm.metadata.get("organism_ch1", [""])[0],
            "platform": gsm.metadata.get("platform_id", [""])[0],
        }
        # Characteristics
        for char in gsm.metadata.get("characteristics_ch1", []):
            if ":" in char:
                k, v = char.split(":", 1)
                row[k.strip().lower().replace(" ", "_")] = v.strip()
        meta_rows.append(row)

    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(out_dir / "sample_metadata.csv", index=False)
    log(f"  Metadata: {len(meta_df)} samples, columns: {list(meta_df.columns)}")

    # ─── Extract expression table ──────────────────────────────────────────
    # Try to build expression matrix from GSM tables
    expr_dict = {}
    for gsm_id, gsm in gse.gsms.items():
        try:
            tbl = gsm.table
            if tbl is not None and len(tbl) > 0:
                # Look for value column (count or expression)
                val_cols = [c for c in tbl.columns if c.upper() not in ["ID_REF", "IDENTIFIER"]]
                if val_cols:
                    expr_dict[gsm_id] = tbl.set_index(tbl.columns[0])[val_cols[0]]
        except Exception:
            pass

    if expr_dict:
        expr_df = pd.DataFrame(expr_dict)
        expr_df.index.name = "gene_id"
        expr_df.to_csv(out_dir / "expression_matrix.csv")
        log(f"  Expression matrix: {expr_df.shape[0]} features × {expr_df.shape[1]} samples")
        return meta_df, expr_df
    else:
        log(f"  WARNING: Could not extract expression table from GSM records")
        log(f"  Raw SOFT file saved in {out_dir} — may need manual processing")
        return meta_df, None


def infer_condition(meta_df, info):
    """Try to assign condition labels (infected vs control) from metadata."""
    condition_keywords = info["condition_key"]
    virus = info["virus"].lower()

    # Try to find condition column
    cond_col = None
    for col in meta_df.columns:
        col_lower = col.lower()
        if any(k in col_lower for k in ["condition", "treatment", "status", "group", "infect", "disease"]):
            cond_col = col
            break

    if cond_col is None:
        # Use title
        cond_col = "title"

    log(f"  Using column '{cond_col}' for condition assignment")
    meta_df["condition"] = "Unknown"
    for idx, row in meta_df.iterrows():
        val = str(row.get(cond_col, "")).lower()
        if virus in val or "infect" in val or "positive" in val or "case" in val or "patient" in val:
            if "mock" not in val and "uninfect" not in val and "healthy" not in val and "control" not in val:
                meta_df.at[idx, "condition"] = "Infected"
        elif "mock" in val or "uninfect" in val or "healthy" in val or "control" in val or "normal" in val:
            meta_df.at[idx, "condition"] = "Control"

    cond_counts = meta_df["condition"].value_counts()
    log(f"  Condition assignment: {cond_counts.to_dict()}")
    meta_df.to_csv(meta_df.index.name and "sample_metadata_labeled.csv" or "sample_metadata.csv", index=False)
    return meta_df


def main():
    ckpt = load_ckpt()

    log("=" * 60)
    log("Step 11: Download Validation Datasets (Phase 8)")
    log("=" * 60)

    for gse_id, info in DATASETS.items():
        if ckpt.get(f"{gse_id}_done"):
            log(f"\n{gse_id}: Already downloaded — skipping")
            continue

        log(f"\n{'=' * 50}")
        log(f"Downloading {gse_id}: {info['description']}")
        log(f"{'=' * 50}")

        out_dir = RAW_DIR / gse_id
        meta_df, expr_df = download_and_parse_gse(gse_id, out_dir)

        if meta_df is not None:
            meta_df = infer_condition(meta_df, info)
            meta_df.to_csv(out_dir / "sample_metadata.csv", index=False)

            ckpt[f"{gse_id}_done"] = True
            ckpt[f"{gse_id}_n_samples"] = len(meta_df)
            ckpt[f"{gse_id}_has_expression"] = expr_df is not None
            save_ckpt(ckpt)
            log(f"  {gse_id}: Download complete")
        else:
            log(f"  {gse_id}: Download FAILED")
            ckpt[f"{gse_id}_done"] = False
            save_ckpt(ckpt)

    # ─── Summary ──────────────────────────────────────────────────────────────
    log("\n" + "=" * 60)
    log("STEP 11 COMPLETE — Validation Dataset Summary")
    log("=" * 60)
    for gse_id in DATASETS:
        done = ckpt.get(f"{gse_id}_done", False)
        n = ckpt.get(f"{gse_id}_n_samples", 0)
        has_expr = ckpt.get(f"{gse_id}_has_expression", False)
        status = "OK" if done else "FAILED"
        log(f"  {gse_id}: {status}  n_samples={n}  expression={'yes' if has_expr else 'no (manual needed)'}")

    log("\nNext: run step12_external_validation_deg.py")


if __name__ == "__main__":
    main()
