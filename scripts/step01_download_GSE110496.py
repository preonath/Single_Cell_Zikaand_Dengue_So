"""
Step 01: Download GSE110496 from GEO
Dataset: Zanini et al., eLife 2018 — DENV-2 + ZIKV single-cell RNA-seq in Huh7 cells
Checkpoint-based: safe to restart if interrupted
"""

import os
import json
import time
import requests
import gzip
import shutil
import subprocess
from pathlib import Path
import GEOparse

# ─── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR   = Path("/home/preonath/Desktop/Preonath_Project/Zika_Dengue")
RAW_DIR    = BASE_DIR / "00_raw_data" / "GSE110496"
CKPT_DIR   = BASE_DIR / "checkpoints"
LOG_FILE   = BASE_DIR / "logs" / "step01_download.log"

RAW_DIR.mkdir(parents=True, exist_ok=True)
CKPT_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

CKPT_FILE  = CKPT_DIR / "step01_checkpoint.json"

# ─── Checkpoint helpers ────────────────────────────────────────────────────────
def load_checkpoint():
    if CKPT_FILE.exists():
        with open(CKPT_FILE) as f:
            return json.load(f)
    return {}

def save_checkpoint(data: dict):
    with open(CKPT_FILE, "w") as f:
        json.dump(data, f, indent=2)

def log(msg: str):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

# ─── Main ──────────────────────────────────────────────────────────────────────
def main():
    ckpt = load_checkpoint()
    log("=" * 60)
    log("Step 01: Download GSE110496")
    log("=" * 60)

    # ── 1. Fetch GEO metadata ──────────────────────────────────────────────────
    if ckpt.get("geo_metadata_done"):
        log("✓ GEO metadata already fetched — skipping")
    else:
        log("Fetching GEO metadata for GSE110496 ...")
        try:
            gse = GEOparse.get_GEO("GSE110496", destdir=str(RAW_DIR), silent=False)
            log(f"  GSE title : {gse.metadata.get('title', ['?'])[0]}")
            log(f"  Organism  : {gse.metadata.get('organism_ch1', ['?'])[0] if 'organism_ch1' in gse.metadata else 'see samples'}")
            log(f"  N samples : {len(gse.gsms)}")

            # Save sample metadata table
            sample_rows = []
            for gsm_name, gsm in gse.gsms.items():
                row = {"gsm": gsm_name}
                for k, v in gsm.metadata.items():
                    row[k] = "; ".join(v) if isinstance(v, list) else v
                sample_rows.append(row)

            import pandas as pd
            meta_df = pd.DataFrame(sample_rows)
            meta_df.to_csv(RAW_DIR / "sample_metadata.csv", index=False)
            log(f"  Metadata saved → sample_metadata.csv ({len(meta_df)} rows)")

            ckpt["geo_metadata_done"] = True
            ckpt["n_samples"] = len(gse.gsms)
            save_checkpoint(ckpt)
        except Exception as e:
            log(f"ERROR fetching metadata: {e}")
            raise

    # ── 2. Download supplementary files ───────────────────────────────────────
    if ckpt.get("suppl_download_done"):
        log("✓ Supplementary files already downloaded — skipping")
    else:
        log("Downloading supplementary files from GEO ...")
        suppl_base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110496/suppl/"

        # Try to list available supplementary files via GEO FTP
        try:
            import urllib.request
            import re
            ftp_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110496/suppl/"
            with urllib.request.urlopen(ftp_url) as response:
                html = response.read().decode()
            # Extract filenames
            files = re.findall(r'href="(GSE110496[^"]+)"', html)
            log(f"  Found {len(files)} supplementary files: {files}")
        except Exception as e:
            log(f"  Could not list FTP directory: {e}")
            # Known files from paper
            files = [
                "GSE110496_RAW.tar",
            ]
            log(f"  Using known filenames: {files}")

        downloaded = ckpt.get("downloaded_files", [])
        for fname in files:
            if fname in downloaded:
                log(f"  ✓ Already downloaded: {fname}")
                continue
            url = suppl_base + fname
            dest = RAW_DIR / fname
            log(f"  Downloading: {fname} ...")
            try:
                _download_file(url, dest)
                downloaded.append(fname)
                ckpt["downloaded_files"] = downloaded
                save_checkpoint(ckpt)
                log(f"  ✓ Saved: {dest} ({dest.stat().st_size / 1e6:.1f} MB)")
            except Exception as e:
                log(f"  ERROR downloading {fname}: {e}")

        ckpt["suppl_download_done"] = True
        save_checkpoint(ckpt)

    # ── 3. Extract TAR if present ──────────────────────────────────────────────
    tar_file = RAW_DIR / "GSE110496_RAW.tar"
    if tar_file.exists() and not ckpt.get("tar_extracted"):
        log(f"Extracting {tar_file.name} ...")
        import tarfile
        with tarfile.open(tar_file) as tar:
            tar.extractall(path=RAW_DIR)
        log("  ✓ Extraction complete")
        ckpt["tar_extracted"] = True
        save_checkpoint(ckpt)
    elif ckpt.get("tar_extracted"):
        log("✓ TAR already extracted — skipping")

    # ── 4. Decompress .gz files ────────────────────────────────────────────────
    if not ckpt.get("gz_extracted"):
        gz_files = list(RAW_DIR.glob("*.gz"))
        log(f"Found {len(gz_files)} .gz files to decompress")
        for gz_path in gz_files:
            out_path = gz_path.with_suffix("")
            if out_path.exists():
                log(f"  ✓ Already decompressed: {out_path.name}")
                continue
            log(f"  Decompressing: {gz_path.name} ...")
            with gzip.open(gz_path, "rb") as f_in:
                with open(out_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            log(f"    → {out_path.name} ({out_path.stat().st_size / 1e6:.1f} MB)")
        ckpt["gz_extracted"] = True
        save_checkpoint(ckpt)
    else:
        log("✓ .gz files already decompressed — skipping")

    # ── 5. List downloaded files ───────────────────────────────────────────────
    log("\nFiles in 00_raw_data/GSE110496/:")
    all_files = sorted(RAW_DIR.iterdir())
    for f in all_files:
        size_mb = f.stat().st_size / 1e6
        log(f"  {f.name:<60} {size_mb:8.2f} MB")

    # ── 6. Quick inspection of data format ────────────────────────────────────
    log("\nInspecting data format ...")
    _inspect_files(RAW_DIR)

    log("\n✓ Step 01 complete. Checkpoint saved.")
    log(f"  Checkpoint file : {CKPT_FILE}")
    log(f"  Log file        : {LOG_FILE}")
    log("\nNext: run step02_load_and_qc.py")


def _download_file(url: str, dest: Path, chunk_size: int = 1024 * 1024):
    """Stream download with progress."""
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        downloaded = 0
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = 100 * downloaded / total
                        print(f"\r    {pct:.1f}% ({downloaded/1e6:.1f}/{total/1e6:.1f} MB)", end="", flush=True)
        print()


def _inspect_files(raw_dir: Path):
    """Print first few lines of text files and shape of matrix files."""
    import pandas as pd

    for fpath in sorted(raw_dir.iterdir()):
        suffix = fpath.suffix.lower()
        name   = fpath.name.lower()

        if suffix in (".txt", ".tsv", ".csv") and fpath.stat().st_size < 500 * 1e6:
            try:
                sep = "\t" if suffix in (".txt", ".tsv") else ","
                df  = pd.read_csv(fpath, sep=sep, nrows=3, index_col=0)
                log(f"  {fpath.name}: shape preview → {df.shape} | cols: {list(df.columns[:5])}")
            except Exception as e:
                log(f"  {fpath.name}: could not parse ({e})")

        elif suffix == ".h5" or ".h5ad" in name:
            try:
                import anndata as ad
                adata = ad.read_h5ad(fpath)
                log(f"  {fpath.name}: AnnData {adata.shape}")
            except Exception as e:
                log(f"  {fpath.name}: not h5ad ({e})")

        elif suffix == ".mtx":
            log(f"  {fpath.name}: MEX format matrix — will load in step02")


if __name__ == "__main__":
    main()
