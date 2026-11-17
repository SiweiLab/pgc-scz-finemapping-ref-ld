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
| `multi_ancestry/aou_ld/` | the cluster pipeline: master variant table, LD builds, harmonisation (`harmonize_v2.py`, `v2_io.py`), the fit driver (`run_primary_2026.py`, parameters in `run_params.py`), exports, and the supplementary analyses |
| `multi_ancestry/aou_ld/task_v2_out/` | per-method PIP and credible-set exports for all 250 loci, which every figure is built from |
| `multi_ancestry/aou_ld/task_*_out/` | exports from the supplementary analyses (see `RUN_ORDER.md`) |
| `multi_ancestry/aou_ld/v2_inputs/` | the master variant table, locus list and coordinates (the GRCh38 target lists are rebuilt by `build_v2_hg38_targets.py`; `RUN_ORDER.md` step 3a) |
| `multi_ancestry/metrics/` | concordance metrics, figures and supplementary tables; `regenerate_v2.sh` rebuilds them all |
| `multi_ancestry/metrics/figures/v2/` | the four hand-made figure inputs; `regenerate_v2.sh` writes every figure and table here |
| `multi_ancestry/multisusie/`, `multi_ancestry/mesusie/`, `susie/` | per-method runners called by the fit driver |
| `multi_ancestry/envs/` | conda environment spec, lock file and setup script |
| `multi_ancestry/inputs/` | where external inputs go (none are redistributed; see its README) |
| `ld/`, `data_prep/`, `gcs_io.py` | GCTB LD reader, sumstats preparation, and storage I/O helpers |
| `RUN_ORDER.md` | every step as it was run, with commands |
| `multi_ancestry/VERSIONS.md` | tool versions, and the defaults each driver leaves unpassed |

## Pipeline design

- **One variant key.** Each locus has a master variant table, keyed on `chr:bp:A1:A2` (GRCh37, alleles sorted). It covers
  every variant with a z-score in the European, East Asian or trans-ancestry summary statistics. Every source joins to it
  by this key, never by rsID, and every z and LD row is signed to the master's reference allele.
- **LD stored once, in the source's own coding.** For each source (UK Biobank GCTB, All of Us EAS, All of Us EUR,
  PolyFun UKBB) we store one matrix per locus, cut to the master variants, together with a per-variant `sign`.
  `v2_io.py` is the only code that applies the sign or touches allele orientation.
- **Everything else is derived at run time.** The merged LD, each method's variant set and the strand-ambiguous filter
  are built when a run starts and are never stored. So is each tool's input format.
- **Two variants of every fit.** `noambig` drops A/T and C/G variants and is the primary analysis, matching the
  benchmark. `primary` keeps them and is a sensitivity analysis.

## Rebuilding the figures and tables (laptop)

Requirements: Python ≥ 3.10 with numpy, pandas, scipy, matplotlib, Pillow and openpyxl.

```
cd multi_ancestry/metrics
OUT=/tmp/regen_check bash regenerate_v2.sh     # writes to /tmp/regen_check; omit OUT to write figures/v2
```

This rebuilds every figure panel, composite and supplementary table from the committed exports, without refitting. Most
outputs need the benchmark (below), and a few read the third-party tables listed in `multi_ancestry/inputs/README.md`.
Built that way, every output matched the figures and tables in the paper exactly (PNG pixels, CSV, XLSX cells).

## Rerunning the fine-mapping (All of Us Researcher Workbench)

The LD builds and fits ran on the All of Us Researcher Workbench (Dataproc, with Hail for the LD extraction). Set up the
environment with `bash multi_ancestry/envs/setup_ilae4_env.sh`. That installs the conda environment plus SuSiEx, FINEMAP,
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
| Benchmark: Yuan *et al.* in-sample SuSiEx output, and the per-variant posterior rebuilt from it | available from the authors of Yuan *et al.* on request; not included |
| Published tables (Trubetskoy 2022 ST11, Yuan 2024 supplementary tables, Mawson 2026 media-2) | from the publishers; links in `multi_ancestry/inputs/README.md` |

No participant identifiers or individual-level data are included.

## Licence

MIT; see `LICENSE`.
