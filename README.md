# Reference LD for multi-ancestry fine-mapping: an empirical benchmark

Code for *An empirical benchmark of reference linkage disequilibrium for multi-ancestry fine-mapping*
(He, Yuan, Zhao, Huang, Neale, Cui & Chen).

We fine-mapped 250 PGC3 schizophrenia loci with reference LD instead of in-sample LD, and scored the results against the
published in-sample SuSiEx analysis of Yuan *et al.* (2024).

- **Ancestry-aware methods, with per-ancestry LD:** SuSiEx, MultiSuSiE and MESuSiE. European LD comes from UK Biobank
  (GCTB, 10,000 participants) and East Asian LD from All of Us v9.
- **Merged-LD methods:** SuSiE-RSS, FINEMAP, RSparsePro and CARMA, run on the trans-ancestry meta-analysis with one
  sample-size-weighted merged LD matrix.

## Repository layout

| Path | What it holds |
|---|---|
| `pipeline/` | the cluster pipeline, in run order (see `RUN_ORDER.md`): master variant table (`build_master_table.py`), LD builds (`build_ld_*.py`), harmonisation (`harmonize.py`, `ld_io.py`), the fit driver (`run_finemap.py`, every parameter in `run_params.py`) and the export (`export_results.py`) |
| `pipeline/methods/` | per-method runners called by the fit driver (SuSiE-RSS, MESuSiE, MultiSuSiE, CARMA) |
| `analyses/` | supplementary cluster analyses, one script per question |
| `benchmark/` | building the published benchmark tables the results are scored against |
| `results/` | every committed output the figures are built from (below) |
| `figures/` | figure and table code; `make_figures.sh` rebuilds them all into `figures/output/`; `hand_made/` holds the two hand-drawn panels (Fig 1a, Supp Fig 1a); legends are in the manuscript |
| `inputs/` | locus lists and coordinates, the master variant table, the GCTB block file and the PolyFun windows; `inputs/external/` holds the benchmark output and three published tables (the tables for review only; see its README) |
| `environment/` | conda environment spec, lock file, setup script, and tool versions (`VERSIONS.md`) |
| `RUN_ORDER.md` | every step as it was run, with commands |

### Results and the figures they feed

| Folder | Contents | Used for |
|---|---|---|
| `results/finemap/main/` | PIPs and credible sets, nine methods × 250 loci, strand-ambiguous variants dropped (the primary analysis) | Figs 1–2; Supp Figs 1, 2, 4–8, 11; Supp Tables 1, 4 |
| `results/finemap/with_ambiguous/` | the same fits with A/T and C/G variants kept | Supp Fig 2b, c |
| `results/finemap/L_sweep/` | each method at L = 1, 3, 5 (and 10) | Supp Fig 3 |
| `results/finemap/coverage_95/` | SuSiE-RSS and FINEMAP with merged LD at 95% coverage | Supp Fig 5 |
| `results/susiex_eur_panel/` | SuSiEx with All of Us instead of UK Biobank European LD | Supp Fig 9 |
| `results/eur_panels/` | European-only SuSiE-RSS under three European LD panels | Supp Fig 10 |
| `results/fig2_loci/` | association statistics, r² and gene tracks for the Fig. 2 loci | Fig 2a–c |
| `results/kriging/` | susieR kriging diagnostics at the Fig. 2b locus | Supp Fig 8a, b |
| `results/s_statistic/` | z/LD consistency s (susieR `estimate_s_rss`), per ancestry and per European panel | Supp Table 2; Supp Figs 8c, 10 |
| `results/lambda_merged/` | SuSiE-RSS λ on the merged LD, per locus | Supp Tables 2, 3 |
| `results/carma_outliers/` | variants CARMA flagged as inconsistent with the merged LD | Supp Table 3 |
| `results/benchmark/` | published PIPs, the r² > 0.95 proxy map and the Mawson *et al.* tables and the full benchmark posterior | scoring throughout; Supp Figs 5, 6 |

## Pipeline design

- **One variant key.** Each locus has a master variant table, keyed on `chr:bp:A1:A2` (GRCh37, alleles sorted). It covers
  every variant with a z-score in the European, East Asian or trans-ancestry summary statistics. Every source joins to it
  by this key, never by rsID, and every z and LD row is signed to the master's reference allele.
- **LD stored once, in the source's own coding.** For each source (UK Biobank GCTB, All of Us EAS, All of Us EUR,
  PolyFun UKBB) we store one matrix per locus, cut to the master variants, together with a per-variant `sign`.
  `ld_io.py` is the only code that applies the sign or touches allele orientation.
- **Everything else is derived at run time.** The merged LD, each method's variant set and the strand-ambiguous filter
  are built when a run starts and are never stored. So is each tool's input format.
- **Two variants of every fit.** `noambig` drops A/T and C/G variants and is the primary analysis, matching the
  benchmark. `primary` keeps them and is a sensitivity analysis.

## Rebuilding the figures and tables (laptop)

Requirements: Python ≥ 3.10 with numpy, pandas, scipy, matplotlib, Pillow and openpyxl.

```
cd figures
OUT=/tmp/fig_check bash make_figures.sh     # writes to /tmp/fig_check; omit OUT to write figures/output
```

This rebuilds every figure panel, composite and supplementary table from the committed exports, without refitting. Most
outputs need the benchmark, and a few read the published tables; both are in `inputs/external/` (see its README).
Built that way, every output matched the figures and tables in the paper exactly (PNG pixels, CSV, XLSX cells).

## Rerunning the fine-mapping (All of Us Researcher Workbench)

The LD builds and fits ran on the All of Us Researcher Workbench (Dataproc, with Hail for the LD extraction). Set up the
environment with `bash environment/setup_env.sh`. That installs the conda environment plus SuSiEx, FINEMAP,
RSparsePro and CARMA at the versions in `VERSIONS.md`. Then follow `RUN_ORDER.md` step by step.

Scripts read the workspace bucket from `WORKSPACE_BUCKET`, which the Workbench sets. Local scratch defaults to `~/work/`.
Some steps were run by hand or from recorded commands rather than by a single script; `RUN_ORDER.md` says which, and it
lists the few commands whose exact arguments were not recorded.

## Data

| Input | Availability |
|---|---|
| PGC3 schizophrenia summary statistics (EUR, EAS and core meta-analysis) | public, from the [Psychiatric Genomics Consortium](https://pgc.unc.edu/for-researchers/download-results/) |
| UK Biobank European LD (GCTB, 10,000 participants) | public, [distributed with GCTB](https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/) |
| PolyFun UK Biobank LD (sensitivity analysis) | public, `https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/` |
| All of Us v9 genotypes and the LD derived from them | Controlled Tier; registered researchers only. Neither genotypes nor derived LD can leave the Workbench, so neither is included. The code to rebuild the LD is. |
| Benchmark: Yuan *et al.* in-sample SuSiEx output, and the per-variant posterior rebuilt from it | included, with the authors' permission (`inputs/external/yuan2024_susiex_output.zip`, `results/benchmark/`) |
| Published tables (Trubetskoy 2022 ST11, Yuan 2024 supplementary tables, Mawson 2026 media-2) | included for review in `inputs/external/`; from the publishers (links in its README), licences to be checked before release |

No participant identifiers or individual-level data are included.

## Licence

MIT; see `LICENSE`.
