# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A linear bioinformatics pipeline (not an app) testing one hypothesis: **do DENV-2 and ZIKV
activate the same host genes when infecting the same human cells?** The primary dataset
(GSE110496, Zanini 2018) is single-cell RNA-seq of Huh7 liver cells. The pipeline derives
**15 shared upregulated genes**, then stacks orthogonal evidence (pathways, miRNA targets,
external datasets, PPI network) and scores everything against **6 quality gates (G1–G6)**.
See `README.md` for the exhaustive per-step walkthrough; this file is the orientation layer.

## Environment & running

- **Conda env is `scipy`** (misleadingly named — it's the full bioinformatics stack:
  scanpy, pydeseq2, gseapy, mygene, GEOparse, networkx, statsmodels, openpyxl, matplotlib-venn).
- **The system `python3` has no pandas.** Use the env interpreter explicitly for any ad-hoc
  inspection: `/home/preonath/miniconda3/envs/scipy/bin/python`, or `conda activate scipy` first.
- Run a single stage directly: `python scripts/stepNN_*.py`. Each script is self-contained
  and hardcodes the absolute `BASE_DIR`.
- There is no test suite, linter, or build. "Correctness" is the gate results.

## Architecture — the two things to understand first

**1. Two parallel representations of the same pipeline (scripts ⇄ notebook).**
`scripts/step01..step17` is the canonical, ordered source of truth. `Zika_Dengue_Complete_Analysis.ipynb`
embeds the same code but **cell order ≠ step number** (e.g. Step 11 runs at cell#5, Step 15 at cell#27),
and some notebook cells were never executed (`execution_count: None`) even though their on-disk
outputs exist — those outputs came from the *script* runs. **When results and code seem to disagree,
trust the scripts and the dated files in `03_results/`, not the notebook's inline outputs.**

**2. Checkpoint-based resumability.** Every step reads/writes `checkpoints/stepNN_checkpoint.json`
with boolean keys (e.g. `qc_done`, `shared_done`). A step with its checkpoint present **skips and
returns early.** To force a stage to recompute, delete its checkpoint JSON. This is why re-running a
script may appear to "do nothing."

## Data flow (numbered folders are pipeline stages, in order)

```
00_raw_data/        GEO downloads (GSE110496 primary + 3 external validation GSEs)
01_processed_data/  anndata_objects/ (.h5ad), pseudobulk/, deg_tables/ (canonical DEG CSVs)
02_literature_resources/  curated host-factor + miRNA gene lists (.txt/.csv) feeding the gates
03_results/         phaseN_*/ — one subdir per analysis phase; *.csv gate outputs live here
04_figures/         main/ (manuscript figs) + supplementary/ + wetlab_validation/
05_manuscript_tables/ Tables 1–4 + styled Supplementary_Tables.xlsx
wetlab_results/      5 curated .xlsx host-factor sheets (lab input, not generated)
```

Stage N consumes stage N-1's folder outputs; nothing writes upstream. `checkpoints/`, `logs/`,
and `SESSION_LOG.md` are cross-cutting run records.

## Gotchas specific to this project

- **"Shared DEGs" has two distinct meanings — don't conflate them.** (a) The strict
  cross-flavivirus intersection = the **15 novel genes** (BIRC3, CCL4, CXCL1, DUSP1, INHBE,
  CREBRF, …) used by gates G4/G5/G6; these have **zero** overlap with the curated wetlab lists
  (that's the finding, not a bug — G4 is *expected* to be non-significant). (b) The per-virus
  DEG lists (`DEGs_DENV_vs_Control.csv`, etc.) **do** overlap known ISGs (IFI6, IFITM1/3, MX1,
  ISG15). Validation figures showing overlaps refer to (b); gate G4=0 refers to (a).
- **DEG tables key on Ensembl IDs (`gene_id`); only the `_annotated` tables carry HGNC `symbol`.**
  Map IDs→symbols before comparing against any gene-name list.
- **MT% is 0.0 for all cells** in GSE110496 (non-standard MT gene naming), so the mito QC filter
  is a no-op here. This is documented/expected — don't "fix" it.
- The `wetlab_results/Medium Confidance.xlsx` sheet has its **real header on row 2** (row 1 is
  `Column1..Column8`) and each `.xlsx` carries an empty trailing `Sheet1`. Parse with `header=1`.
- `02_literature_resources/host_factors_*.txt` (fed to the Fisher gates) are **high-confidence
  only**; the 25 medium-confidence genes live in `host_factors_curated.csv` but are deliberately
  excluded from the gate `.txt` lists.

## External APIs the pipeline calls (network-dependent, rate-limited)

MyGene.info (ID annotation, batched 1000 + 0.5s sleep), Enrichr via gseapy (pathway & miRNA-target
enrichment), STRINGdb REST (PPI network), NCBI GEO via GEOparse (all dataset downloads). A failure
mid-pipeline is usually one of these timing out, not a logic bug — re-run the step (checkpoints
preserve prior progress).

---

## The findings (what the pipeline actually concluded)

The headline result is **15 shared upregulated genes** induced by *both* DENV and ZIKV in the same
hepatocytes, and the claim that they are **novel convergent host-response genes** rather than
rediscoveries of known antiviral/proviral factors. The 6 gates are the evidence chain for that claim:

| Gate | Question it answers | Result | Why it matters |
|------|---------------------|--------|----------------|
| **G1** | Did each virus produce enough DEGs to analyze? | DENV 527 / ZIKV 176 — BORDERLINE | Sanity floor; ZIKV is weaker but usable |
| **G2** | Is the *overlap* between DENV and ZIKV DEGs more than chance? | p=2.26×10⁻¹⁵, 18.6× — **PASS** | This is the core "convergence is real" test |
| **G3** | Do the two viruses' genome-wide fold-changes correlate over time? | r=0.462 at 48h — **PASS at 48h** | Convergence *emerges late*, not early (4–24h are flat) |
| **G4** | Are the 15 shared genes known host factors? | 0/15, p=1.0 — NOT PASSED | **A "failed" gate that is the point: novelty.** Not a bug. |
| **G5** | Are they targets of flavivirus-suppressed miRNAs? | miR-15a-5p p=0.015; CREBRF hit by 10 miRNAs — TREND | Proposes a *mechanism* (miRNA de-repression) for the convergence |
| **G6** | Do they replicate in an independent ZIKV dataset? | ZIKV NPCs 4/15, p=1.18×10⁻⁴ — **STRONG PASS** | Replication in a different cell type rules out a Huh7 artifact |

**How to read this honestly (and how the project frames it):** the strong gates (G2, G3@48h, G6)
support genuine, replicating convergence. G4=0 is *reframed as the discovery* — these genes sit
outside curated host-factor databases, so they are candidates, not confirmations. G1 borderline and
G5 only-a-trend are the honest weak points; nothing here is overstated to "significant." The lead
candidate is **CREBRF** (an ER-stress / ATF6α regulator): it replicates in NPCs (G6) *and* is the
top miRNA hub (G5), so it carries two independent lines of orthogonal evidence.

**Why the negative validations (GSE94892 PBMCs, GSE118305 macrophages → 0% replication) are kept,
not hidden:** they are cell-type-specificity controls. The signature replicates in NPCs (another
infectable parenchymal/neural cell) but not in blood/myeloid cells, which is the *expected* pattern
if this is a hepatocyte-like intrinsic response rather than a generic immune-cell program.

## Why these methods were chosen (rationale behind the approach)

- **Pseudobulk + pyDESeq2 instead of single-cell DE.** Counts are summed per
  (condition × timepoint) into 12 pseudobulk samples before differential testing. Single-cell DE
  tests treat each cell as a replicate and massively inflate significance (pseudoreplication);
  pseudobulk restores the real experimental unit (the sample) so the negative-binomial model and
  p-values are statistically valid. This is the current best-practice for designs with real
  replicates.
- **Fisher exact test for every "overlap" question (G2, G4, G6).** Each is the same shape: "is the
  intersection of two gene sets bigger than expected given the transcriptome background?" Fisher on a
  2×2 (in-A/not-A × in-B/not-B) is the exact, assumption-light answer and avoids χ² approximation
  issues at small overlaps. The background size (~22,369 genes) is the denominator that makes
  "enrichment" meaningful.
- **MOI=1 fairness re-run (step03b).** DENV was assayed at multiple MOIs but ZIKV only at MOI=1.
  Re-deriving DENV DEGs from MOI=1 cells only removes "different viral dose" as an alternative
  explanation for the shared signature — a built-in robustness control, not a separate analysis.
- **The 36-miRNA "DENV-only" set in G5 is a deliberate negative control.** G5 isn't just "do
  flavivirus miRNAs target these genes" — it's "do the *cross-flavivirus* (55) miRNAs target them
  *more* than DENV-specific (36) miRNAs that can't bind the ZIKV genome." The control set is what
  turns a correlation into a specificity argument.
- **CPM + pseudocount for the temporal trajectory (step07), but pyDESeq2 for the headline DEGs.**
  G3 needs a *per-timepoint, genome-wide* fold-change at 4 timepoints where full DESeq2 modeling
  per slice is overkill and unstable; a CPM ratio with `+1` pseudocount (to avoid log(0)) is the
  lightweight, appropriate tool for a correlation. The pseudocount choice trades a little bias for
  numerical stability — standard for trajectory-style comparisons.
- **STRINGdb with a ring layout for isolated nodes (step13).** 11/15 genes have *no* STRING edges —
  which is itself consistent with novelty (uncharacterized genes lack curated interactions). The
  ring layout is a presentation choice so "no known interactions" is shown honestly rather than
  hidden by a force-directed layout that would scatter them ambiguously.

## What the recurring code idioms mean

- **The checkpoint pattern (`load_ckpt`/`save_ckpt` + early-return on a boolean key).** These
  scripts call slow, rate-limited external APIs (GEO, MyGene, Enrichr, STRING). The checkpoint makes
  every step *idempotent and resumable*: a crash on step 12 doesn't cost you steps 1–11. The mental
  model is "make targets," done in plain JSON. **Consequence for you: a script that "does nothing" is
  skipping because its checkpoint says done — delete the JSON to force recompute** (this is exactly
  what the step08/09 refresh required).
- **`alternative="greater"` in every `fisher_exact` call.** The biological hypothesis is *enrichment*
  (more overlap than chance), a one-sided question. A two-sided test would waste power testing for
  depletion, which is not the hypothesis.
- **`design_factors=["timepoint", "condition"]` in DeseqDataSet.** Timepoint is a nuisance/blocking
  variable; putting it in the model *adjusts out* time-driven expression so the `condition` contrast
  (DENV vs Control) isn't confounded by when the sample was taken. The contrast
  `["condition","DENV","Control"]` then reads as "virus effect, holding timepoint constant."
- **Batched MyGene queries with `time.sleep(0.5)`.** Politeness/rate-limit avoidance against a public
  API; batching 1000 IDs per call is the throughput knob, the sleep is the don't-get-blocked knob.
- **`regress_out(["pct_counts_mt"])` then `scale(max_value=10)` in QC.** Standard scanpy hygiene:
  remove mito-fraction as a technical covariate, then z-score and clip extreme values so a few
  outlier cells don't dominate PCA/UMAP. Here MT% is all-zero so the regress-out is inert (see
  gotchas) — the line stays for pipeline uniformity, not effect.
- **Two annotated DEG table variants (`_annotated` vs `_annotated_full`).** `full` keeps every gene
  (for background sets and volcano plots); the non-full is filtered to significant DEGs (for
  reporting). Pick by whether you need the universe or the hits.

---

## Session additions — 2026-06-21 (wetlab↔computational validation, 3D UMAP, trajectory)

This block records analyses added after the original pipeline. There is a **second working copy /
git repo at `/home/preonath/Desktop/Preonath_Project/Part_Zika_To_git/`** that mirrors selected
folders (`wetlab_validation/`, `03_results/`, `04_figures/`, `scripts/`…). New artifacts are produced
in the canonical `Zika_Dengue/` tree and **copied** into `Part_Zika_To_git/` for committing. As of
this session those copies are staged-but-**not committed**.

### Reconciling the headline DEG counts (DENV 527 / ZIKV 176)
The official per-virus DEG count = **`padj < 0.05` AND `|log2FC| ≥ 1`** (constants `P_THRESHOLD`,
`FC_THRESHOLD` in `step03_pseudobulk_deg.py`). This is what `volcano_DENV.png`/`volcano_ZIKV.png`
(supplementary) and Gate G1 report.
- DENV = **527** (300 up / 227 down). Raw `padj<0.05` *without* the FC cutoff = 1442 — don't quote that.
- ZIKV = **176** (66 up / 110 down) — down-skewed, the borderline dataset.
Result files: `01_processed_data/deg_tables/DEGs_{DENV,ZIKV}_vs_Control.csv`. The `log2FoldChange`
column is named `log2FoldChange` in the `_annotated_full` tables (the older `_annotated` ZIKV table
uses `log2FC` — watch the rename).

### Category-aware wetlab validation (reframes G4) — `03_results/phase5_validation/`
G4's "0/15, p=1.0" lumps all wetlab categories together, which **buries the antiviral signal under
proviral noise**. The correct test is **direction-aware and category-stratified**, because DEGs and
host-factor screens measure different things: DEGs detect *transcriptional induction* (antiviral/ISGs
are induced → detectable) whereas proviral factors are constitutively-expressed *functional*
requirements (CRISPR/RNAi hits) invisible to a DEG test. Splitting by category:
- **Antiviral/ISG: OR = 5.76, p = 4.4×10⁻⁴ (ENRICHED ✅)** — pipeline recovers IFI6, IFITM1/3, ISG15,
  MX1, TRIM56, RETREG1.
- **Proviral: OR = 1.87, p = 0.14 (ns)** — category mismatch, *not* a failure; transcriptomics can't
  measure functional dependency.
Outputs: `wetlab_vs_computational_pergene.csv` (per-gene DENV/ZIKV log2FC+padj+verdict),
`_summary.csv`, `_fisher.csv`. Figure: `04_figures/wetlab_validation/Figure_Wetlab_Computational_Validation.png`
(3 panels). This is the honest reframe: agreement on *known* antiviral genes validates the method;
the 15 shared genes being *outside* the lists is the novelty claim, not a bug.

### DEG ∩ wetlab overlap volcanoes — `04_figures/wetlab_validation/`
Volcano = full per-virus DEG cloud with genes that **also appear in a wetlab list** highlighted
(red = antiviral, blue = proviral) + labelled. Two families:
- **vs the virus's own list:** `volcano_{DENV,ZIKV}_wetlab_overlap.png` → DENV 5 genes, ZIKV 1 (RRBP1).
- **vs the full 128-gene wetlab union:** `volcano_{DENV,ZIKV}_vs_wetlab128_overlap.png` → DENV **9**
  (adds MX1, ISG15, RRBP1 because the 128 union includes genes curated on the *other* virus's sheet),
  ZIKV **1** (RRBP1). Companion `common_*_DEG_vs_wetlab*.csv`. All antiviral overlaps are UP.

### Where the "128" wetlab union comes from (Venn `Figure_SharedGenes_vs_Wetlab_Venn.png`)
**128 = union of all 5 wetlab sheets = 106 (4 high-confidence) + 24 (Medium-Confidence, dedup).**
CRITICAL DRIFT: the Venn (Jun 13) parsed Medium with `header=1` (correct → 24 genes → 128). The
**current `step06_wetlab_validation.py` reads Medium with default `header=0`**, which hits the wrong
column (per the Medium-header-on-row-2 gotcha) and finds only **6** → re-running step06 today yields
**112, not 128**. The figure is right; the script drifted. Do NOT confuse the 128 with
`02_literature_resources/host_factors_curated.csv` (123 genes — a different curated superset).
When matching wetlab symbols to DEGs, **expand aliases** (`BST2 (Tetherin)` → BST2 + TETHERIN;
128 raw → 145 match keys) or you'll miss real hits like MX1/ISG15.

### 3D UMAP (interactive HTML) — `scripts/step02b_umap3d.py`
Recomputes UMAP with `n_components=3` on the **existing neighbors graph** from step02 (so it's
consistent with the 2D embedding), into a **separate** key `obsm['X_umap_3d']` — the canonical
`adata_processed.h5ad` and its 2D `X_umap` are **never overwritten**. Outputs to
`04_figures/supplementary/`: `UMAP_3D_{condition,timepoint,viral_load}.html` (Plotly, rotatable) +
`umap_3d_coords.csv`. Checkpoint `checkpoints/step02b_checkpoint.json` (early-return when
`html_done`). **Plotly was pip-installed into the `scipy` env this session** (was missing).

### Trajectory analysis verdict — viral-load trajectory, NOT developmental pseudotime
Asked whether single-cell trajectory analysis is valid here. Deep answer, grounded in a diagnostic on
the 2075 cells:
- **Infection is the dominant axis**: in DENV cells `viral_molecules` ↔ PC1 Spearman **ρ = 0.65**.
- **No depth confound**: viral_load vs `total_counts` ρ = −0.08 (DENV) / −0.22 (ZIKV) — the negative
  sign is real **viral host-shutoff**, not library-size artifact.
- **The 15 shared genes dose-respond to viral load in BOTH viruses**: CREBRF, TSPYL2, INHBE, DUSP1 all
  rise with load in DENV *and* ZIKV (DENV 6/15, ZIKV 8/15 at ρ>0.2).
**Guidance:** do a **supervised viral-load-ordered "infection-progression" trajectory** (cite Zanini
2018, same dataset), regress out cell-cycle/depth, handle bystander (load=0) cells separately, treat
DENV as the powered analysis and ZIKV as supportive (median ZIKV load = 4 mol, 60% positive). Do NOT
run unsupervised Monocle/DPT branching framed as differentiation — Huh7 is a homogeneous line with no
lineage tree. This per-cell, dose-resolved convergence is *stronger* evidence than the pseudobulk G3
time-course and directly rebuts the "late bulk artifact" critique.

**BUILT — `scripts/step18_viral_load_trajectory.py`** (checkpoint `step18_checkpoint.json`; canonical
h5ad read-only, DPT on a copy). Supervised ordering of cells by `viral_molecules`, infected cells
only (bystanders = load-0 excluded), lowess trends of the 15 shared genes vs log10 load.
Cross-checks (all confirm infection is the axis, not artifact): viral↔PC1 ρ=0.65 (DENV) / PC2 ρ=0.72
(ZIKV); an unsupervised **DPT rooted at the lowest-load cell recovers load (ρ=0.59 DENV / 0.72 ZIKV)**;
depth confound ρ=−0.08/−0.22 (negative = host shutoff); cell-cycle G2M ρ=−0.28/−0.32 (infected cells
arrest, not a positive confound). **Genes rising with load in BOTH viruses (ρ>0.2, infected cells):
CREBRF, INHBE, TSPYL2** (DUSP1 borderline 0.20/0.38). CREBRF — the lead candidate — rises in both
(0.23 DENV / 0.35 ZIKV), now carrying a THIRD orthogonal line of evidence (G5 miRNA hub + G6 NPC
replication + per-cell viral-load dose-response). Stats: `03_results/phase_trajectory/viralload_trajectory_stats.csv`.
Figures (`04_figures/supplementary/`): `trajectory_convergence_genes.png` (headline, 4 genes ×
2 viruses), `trajectory_shared_genes_{DENV,ZIKV}.png`, `trajectory_heatmap_{DENV,ZIKV}.png`.
