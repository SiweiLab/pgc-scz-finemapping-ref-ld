# Run order

The steps that produced the paper's numbers and figures, in the order they were run. For the design, see
[README.md](README.md#pipeline-design). Paths are relative to the repository root; cluster commands run from `pipeline/` (supplementary runs from `analyses/`).

- **Where.** Steps 1–7 and the supplementary runs ran on the All of Us Researcher Workbench (Dataproc; Hail in step 3 only).
  Steps 3a and 8 ran on a laptop. `~/work/` is the cluster scratch root.
- **Scratch names.** Scratch folders under `~/work/` and checkpoint names are given generically here; as run they
  carried internal batch labels.
- **Buckets.** `${WORKSPACE_BUCKET}` is the workspace bucket (set by the Workbench; every script reads it from the
  environment). Below, `<root>` = `${WORKSPACE_BUCKET}/manuscript_final_2026-09`.
- **Controlled-tier data** (AoU genotypes, the participant exclusion list) stay in the Workbench. No participant ID is in
  this repository or in any output.
- **Software.** Versions, and the defaults each driver leaves unpassed: [`environment/VERSIONS.md`](environment/VERSIONS.md).
  Every fit parameter: `pipeline/run_params.py`. Conda env: `environment/finemap_env.yml` (lock:
  `finemap_env.lock.yml`), built by `bash environment/setup_env.sh` (prefix `~/envs/finemap_env`; also
  installs SuSiEx, FINEMAP and RSparsePro, and pins CARMA at `856f0fb`).
- **Two variants of every fit:** `noambig` (A/T and C/G variants dropped; the primary analysis, matching the benchmark's
  domain; results in `results/finemap/main/`) and `primary` (kept; sensitivity; `results/finemap/with_ambiguous/`).
  They share LD and differ only by `--drop-ambiguous`.

## 1. Master variant table (cluster)

**Purpose.** One table per locus, one row per variant in any of the EUR, EAS or core PGC3 sumstats, keyed on
`vkey = chr:bp:A:B` (GRCh37, sorted alleles); every later join uses this key and signs to its `ref_a1`.
```
python build_master_table.py SS_ROOT ../inputs/loci_250.txt OUT_DIR        # (inferred: the script's usage line)
```
- **In:** `<root>/PGC_SCZ_sumstats/{eur,eas,core}_sumstats/` (copied to `~/work/PGC_SCZ_sumstats/`). **Out:** `<root>/master_v2/`,
  `inputs/master_variants_250.tsv.gz` (committed).
- **Check:** identical row for row to `inputs/master_expected_laptop.tsv.gz` (built on the laptop from the public PGC3
  files); 272,841 rows; non-palindromic set = Yuan et al.'s 231,272.

## 2. EUR LD from the GCTB UKB 13M FullLDM (cluster, no Hail)

**Purpose.** Range-fetch only the 229 GCTB blocks the loci touch, and cut each locus over the master variants GCTB carries,
joined on position and unordered alleles.
```
python3 fetch_gctb_blocks.py --out <local dir>     # resumable; --list-only prints the plan
python build_ld_eur_gctb.py MASTER_TSV GCTB_DIR OUT_DIR [--workers N]   # (inferred: usage line)
```
- The fetch was split by hand over 6 parallel workers; rerunning the plain command afterwards fetched nothing.
- **In:** `https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/ukbEUR_13M_FullLDM.zip`,
  `inputs/ref_b37_1588blocks.pos`. **Out:** `ld_v2/eur_gctb/{locus}/{R.npz,snps.tsv.gz}` (R in GCTB coding, `sign`
  column), `gctb_presence.tsv.gz`; block cache in `${WORKSPACE_BUCKET}/source_cache/gctb_ldm13M/`.
- **Check:** three server members md5-identical to the copy the v1 build used; 235,941 variants; over v1's variants, signed
  R = v1's to float32.

### 2b. EUR LD from PolyFun/UKBB-337k (cluster, supplement only)

**Purpose.** The same cut from the Broad UKBB LD windows, streamed one window at a time.
```
build_ld_eur_polyfun.py ../inputs/master_variants_250.tsv.gz ~/work/ld_eur_polyfun ~/work/ld_eur_polyfun_scratch --workers 2
```
- Ran in two sessions; the second reran the same command, which skips finished windows.
- **In:** `https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/{window}.{npz,gz}`, windows from
  `inputs/polyfun_windows.tsv`; 249 loci (`chr3_48182866_53553354` is wider than any window). **Out:** `ld_v2/eur_polyfun/`.
- **Check:** signed rows correlate positively with GCTB's; coverage and s as expected; stored as built, **not PSD** at most
  loci.

## 3a. GRCh38 target list (laptop)

**Purpose.** Lift every master row hg19 → hg38 (`pyliftover`; the cluster has no egress for the UCSC chain), with `vkey`
as the Hail row ID and alleles complemented where the segment is inverted.
```
python3 build_hg38_targets.py --master ../inputs/master_expected_laptop.tsv.gz \
    --out ../inputs/hg38_targets
```
- **Out:** `inputs/hg38_targets/{locus}_hg38_targets.tsv`, `intervals.tsv`, `liftover_all.tsv.gz`. Not committed: this command
  rebuilds all 252 files byte-identically from the committed master table (needs `pyliftover`, which downloads the UCSC chain).
- **Check:** 272,758 of 272,841 lifted to their own chromosome, 250/250 loci, 0 collisions.

## 3. AoU v9 EAS and EUR LD (cluster, Hail)

**Purpose.** One Hail pass over the AoU v9 ACAF split MatrixTable at the target positions, then per-locus ALT-dosage
correlation in numpy (replicating `hl.row_correlation`) and alignment of each row to the master.
```
hailctl config set gcs_requester_pays/project $GOOGLE_PROJECT
python build_ld_aou.py --out ~/work/ld_aou --ckpt <root>/ld_v2_extract/eas_targets.mt          # (inferred)
python build_ld_aou.py --out ~/work/ld_aou --ckpt <root>/ld_v2_extract/eas_targets.mt --skip-a # (inferred)
python build_ld_aou.py --ancestry-code eur --expected-n 267101 --pilot-check eur_gctb \
    --ckpt <root>/ld_v2_extract/eur_targets.mt --out <dir>
```
- EAS ran twice: the second run (`--skip-a`, reusing the Stage A checkpoint) adds the filtered-entries fix.
  EUR used the recorded flags above. Defaults applied: call rate ≥ 0.95, `--n-partitions 1000`,
  `--collect-gb 2`. Stage D calls `align_ld_aou.py` itself. The cluster autoscaled from 2 to about 21 nodes.
- **Genotypes:** `gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt`.
  **Samples:** CDR v9 `ancestry_pred_other` minus the exclusion list: EAS 12,374, EUR 267,101.
- **Exclusion list** (`<root>/aou_samples/combined_exclude_samples.txt`, workspace only): the union (44,206 IDs) of the CDR v9
  `aux` files `relatedness_flagged_samples.tsv`, `qc/flagged_samples.tsv` and `qc/wgs_sex_concordance_exception_manifest.tsv`.
- **Target list** (`<root>/ld_v2_extract/all_targets.tsv`): the 250 committed step-3a files in one table (inferred).
- **Out:** `ld_v2/{eas,eur}_aou/{locus}/` (`R.npz` in Hail's coding; `snps.tsv.gz` with `sign`, `strand`, AF, call rate).
- **Check:** at 5 pilot loci R equals the per-locus Hail build to ≤ 1e-6; PSD at 250/250; AoU EUR vs GCTB median row correlation 0.999.

**Promotion.** Each staged tree was copied to `<root>/ld_matrices/ld_per_anc_{eur_gctb,eas_aou,eur_aou,eur_polyfun}/`:
```
gsutil -m rsync -r <staging>/ <dest>/
gsutil -m rsync -r -n -c <staging>/ <dest>/       # must print zero operations
```
The promoted trees and `PGC_SCZ_sumstats/` were then synced by hand to `~/work/` (`ld_matrices/`, `PGC_SCZ_sumstats/`); that copy is the `LD_ROOT` /
`SS_ROOT` of every later step. **Check:** `python test_ld_io.py` (7/7); every source × locus loads through `ld_io`; round
trip against v1's signed EUR LD gives max |Δ| = 0.

## 4–6. Harmonise every arm's inputs (cluster)

**Purpose.** `harmonize.py` reads only through `ld_io.py` (the only code that flips alleles), builds each arm's domain as
a set expression over the master flags, computes the ESS-weighted merged LD at run time, and writes tool inputs to scratch.
```
python test_ld_io.py && python test_harmonize.py          # 7 + 6
bash run_harmonize.sh noambig
bash run_harmonize.sh primary
```
- The wrapper runs `harmonize.py --loci <L> --out ~/work/v2_harmonized/<variant> [--drop-ambiguous]` per locus, 6 in
  parallel, with `LD_ROOT=~/work/ld_matrices SS_ROOT=~/work/PGC_SCZ_sumstats`.
- **MultiSuSiE censored inputs**: the paper's MultiSuSiE reads `union_censored_bothld/` (both-LD domain; a
  missing z gets z = 0 and a zero R row/column in that ancestry). First built alone with
  `python harmonize.py --loci-file ../inputs/loci_250.txt --out ~/work/v2_harmonized/noambig --arms multisusie --drop-ambiguous`
  (and `primary` without the flag), then both trees were rebuilt with all arms by the wrapper above.
- **Out:** `~/work/v2_harmonized/{noambig,primary}/{union,union_censored_bothld,intersect,merged,eur,_logs}/`.
- **Check:** `python check_harmonize.py --run-dir ~/work/v2_harmonized/primary --v1-merged <local v1 ld_eur_eas_merged> --ld-root <local ld_matrices> --loci-file ../inputs/loci_250.txt --out ~/work/check_harmonize`
  (and `--noambig`): every domain equals its set expression at 250/250, the run-time merge reproduces v1's merged LD, and
  the all-arm rebuild is byte-identical at 5 loci per tree.

## 7. Fit: nine arms × two variants (cluster)

**Purpose.** `run_finemap.py` fits SuSiEx, MultiSuSiE (censored, both-LD), MESuSiE, SuSiE-merged, FINEMAP-merged,
RSparsePro, CARMA, EUR-SuSiE-RSS and EUR-FINEMAP with `run_params.py`; failures are logged per (arm, locus), never retried.
```
export HARMONIZED_INPUTS=~/work/v2_harmonized/noambig
python run_finemap.py --dry-run --loci chr12_110129434_111489434        # read every argv
python run_finemap.py --loci chr12_110129434_111489434 --no-backup      # smoke test, 9 arms
python run_finemap.py --jobs <cores-2>                                  # full run
```
Then the same for `primary`. `WORK` = `~/work/primary_v2_<variant>`; backup to `<root>/fine_mapping_runs/primary_v2_2026-09/<variant>/`;
loci file: as run, the PGC3 locus table `~/work/ref_inputs/loci_250_pgc3scz.tsv`; the release default
`inputs/loci_250_coords.tsv` holds the same 250 windows (chr, start, stop).

**Seeded CARMA**. CARMA is stochastic, so the first, unseeded pass is not used:
`run_params.SEED_CARMA = 42` was added and CARMA rerun for both variants, launched detached:
```
for V in noambig primary; do
  HARMONIZED_INPUTS=~/work/v2_harmonized/$V setsid nohup python3 run_finemap.py --arms CARMA \
      --loci chr3_48182866_53553354 --jobs 1 &> ~/work/carma_chr3_$V.log &
done
HARMONIZED_INPUTS=~/work/v2_harmonized/$V setsid nohup python3 run_finemap.py --arms CARMA --loci <comma list> --jobs 6 &> ~/work/carma_resume_$V.log &
```
The driver skips loci with output, so `<comma list>` is the loci still missing, chr3 excluded. **Check:** every argv ends
`... 0.99 10 10 42`; two seeded runs of one locus give byte-identical PIPs.

**L sweep**. L = 1, 3, 5 (plus L = 10 for SuSiEx with mult-step off, and for FINEMAP-merged);
the primary run supplies every other arm's own L. CARMA ran on a second cluster from inputs staged in `<root>/run_inputs_v2/`:
```
B=<root>/fine_mapping_runs/l_sweep_v2_2026-10
for V in noambig primary; do for L in 1 3 5; do
  HARMONIZED_INPUTS=~/work/v2_harmonized/$V setsid nohup python3 run_finemap.py --arms CARMA --L $L --jobs 4 \
      --backup-prefix $B/$V &> ~/work/carma_lsweep_${V}_L${L}.log &
done; done
setsid nohup bash overnight_lsweep_carma.sh &> ~/work/overnight_lsweep_carma.log &     # CARMA cluster: waits, exports the sweep
setsid nohup bash overnight_lsweep.sh &> ~/work/overnight_lsweep.log &     # main cluster
```
- `overnight_lsweep.sh` waits for the chr3 CARMA fits, exports the primary run (below), runs `run_finemap.py --arms <arms>
  --L <n> --jobs 7 --backup-prefix $B/$V` for the six other arms, and exports them.
- As run, both scripts ended by committing their exports to the working repository and stopping their cluster; those lines
  are commented out in this release. Both clusters once stopped early; everything was relaunched unchanged and resumed from disk.
- **Check:** 250/250 loci per (arm, L, variant); argvs carry the sweep L and seed 42.

**Coverage 0.95, merged-LD arms**. SuSiE-merged L = 5 and FINEMAP-merged k = 5, `noambig`, for the
comparison with published 95% sets; `run_params.COVERAGE` unchanged; backup `<root>/fine_mapping_runs/cov095_2026-10/`.
```
HARMONIZED_INPUTS=~/work/v2_harmonized/noambig python3 run_finemap.py --arms SuSiE-merged --L 5 --coverage 0.95   # (inferred)
HARMONIZED_INPUTS=~/work/v2_harmonized/noambig python3 run_finemap.py --arms FINEMAP-merged --coverage 0.95       # (inferred)
```
**Check:** recorded argvs end `... 0.95 0.5 5` and `--n-causal-snps 5 ... --prob-cred-set 0.95`; SuSiE PIPs bit-identical
to the L = 5 sweep; FINEMAP k* identical at 250/250.

### 7 (export). PIP and credible-set tables → `results/finemap/`

**Purpose.** Read each arm's raw output, attach master alleles and PGC3 rsIDs, compute `common_purity` through `ld_io`.
```
HARMONIZED_INPUTS=~/work/v2_harmonized/noambig LD_ROOT=<ld_matrices> python export_results.py   # -> results/finemap/main/
HARMONIZED_INPUTS=~/work/v2_harmonized/primary LD_ROOT=<ld_matrices> python export_results.py   # -> results/finemap/with_ambiguous/
HARMONIZED_INPUTS=~/work/v2_harmonized/$V LD_ROOT=<local ld_matrices> python3 export_results.py --L $L --arms CARMA
```
- Sweep → `results/finemap/L_sweep/` (`L_sweep_with_ambiguous/` for `primary`, not shipped); coverage 0.95 (`--L 5 --coverage 0.95`, inferred) → `results/finemap/coverage_95/`.
- **Out:** `primary_<arm>[-noambig][_L{n}]_{pip,cs_purity}.csv.gz`, `primary_status*.csv.gz`.
- **Check:** `check_exports.py` (laptop; it also needs the `primary_status*.csv.gz` files, which are not shipped): 250/250 loci per arm, rows = logged domains, 0 rsIDs outside the master, 0 NaN PIPs,
  0 unconverged SuSiE-type fits; seeded vs unseeded CARMA in its `gate7_carma_seeded_vs_unseeded.csv`.

## Supplementary runs read by the paper (cluster; no new LD)

All with `LD_ROOT=~/work/ld_matrices SS_ROOT=~/work/PGC_SCZ_sumstats`.

- **Three EUR panels**: EUR-SuSiE-RSS, primary parameters, under GCTB-10k, AoU EUR and PolyFun on a common
  non-palindromic set, 249 loci, PolyFun as built. Out `results/eur_panels/`. Check: 747/747 converged.
  ```
  python eur_panels.py build  --loci-file ../inputs/loci_250.txt --out ~/work/eur_panels --jobs 6
  python eur_panels.py fit    --loci-file ../inputs/loci_250.txt --out ~/work/eur_panels --jobs 14
  python eur_panels.py export --loci-file ../inputs/loci_250.txt --out ~/work/eur_panels
  ```
- **SuSiEx with AoU EUR LD**: SuSiEx refitted with AoU EUR in place of GCTB, EAS unchanged, `noambig`, scratch only.
  Out `results/susiex_eur_panel/`. Check: the GCTB side reproduces the primary SuSiEx export exactly.
  ```
  python harmonize.py --arms susiex --eur-ld eur_aou --drop-ambiguous --loci-file ../inputs/loci_250.txt --out ~/work/v2_extra/noambig_euraou   # (paths inferred)
  HARMONIZED_INPUTS=~/work/v2_extra/noambig_euraou PRIMARY_WORK=~/work/extra_susiex_euraou_noambig python run_finemap.py --arms SuSiEx --no-backup
  python susiex_eur_panel.py GCTB_WORK GCTB_TREE AOU_WORK AOU_TREE OUT_DIR       # (inferred: usage line)
  ```
- **λ on the merged LD**: SuSiE-merged's recorded `estimate_s_rss` λ per locus →
  `results/lambda_merged/`. `python lambda_merged.py` (no arguments; inferred). This release writes the `lambda` column only.
- **susieR kriging at SuSiEx's false-positive loci**: 6 watched loci, 5 controls, exact SuSiEx inputs → `results/kriging/`.
  ```
  Rscript kriging.R <locus> ~/work/v2_harmonized/noambig/union/<locus>_susiex ../results/kriging 58749 14425
  ```
- **Fig. 2 locus data**: signed z, and r² under GCTB, AoU EAS and the merged LD → `results/fig2_loci/`.
  The `*_genes.csv` tracks were fetched on the laptop (`figures/fetch_gene_track.py`).
  ```
  python3 fig2_locus_data.py --tree ~/work/v2_harmonized/noambig \
      chr12_110129434_111489434 chr18_52716708_53630258:rs138834542 chr12_29898169_29948269:12:29935396_C_A
  python3 fig2_locus_data.py --tree ~/work/v2_harmonized/noambig --out ../results/fig2_loci \
      chr3_36834083_36998583:rs75968099 chr1_243512441_244077592:rs61833239
  ```
- **Seeded CARMA outliers**: `Outliers$Index` of the final seeded `noambig` run, mapped through the `.z`
  files CARMA read → `results/carma_outliers/`. Check: all 250 `.rds` postdate the seeded relaunch.
  ```
  Rscript carma_outliers.R ~/work/primary_v2_noambig/CARMA ~/work/v2_harmonized/noambig/merged ../results/carma_outliers/carma_outliers_seeded_noambig.tsv
  ```
- **Per-ancestry s, susieR `estimate_s_rss`**, via `pipeline/susier_s.py`, a port of susieR 0.14.2
  that matches R to 6e-9 → `results/s_statistic/`. Check: the `s_legacy` column reproduces the earlier s (max |Δ| 8e-14).
  ```
  python s_per_ancestry.py --loci-file ../inputs/loci_250.txt --out ../results/s_statistic
  ```
- **s for the three EUR panels, susieR** → `results/s_statistic/`. Check: the `s_legacy` column equals the three-panel run's own s to ~1e-12.
  ```
  python s_eur_panels.py --loci-file ../inputs/loci_250.txt --out ../results/s_statistic --jobs 6
  ```

## 8. Figures and tables (laptop)

**Purpose.** Rebuild every figure panel, composite and supplementary table from the committed results; no refits.
```
cd figures
bash make_figures.sh                     # writes figures/output/
OUT=/tmp/fig_check bash make_figures.sh  # writes elsewhere
```
- **In:** `results/`, `inputs/loci_*.txt`, `inputs/external/` (Yuan et al.'s SuSiEx output is not distributed; see
  `inputs/external/README.md`). Hand-made, no script: `figures/hand_made/` (`fig1A.png`, `supp fig1A.png`). Legends are in the manuscript.
- **Check:** every output was pixel-identical (PNG) or exactly equal (CSV, XLSX) to the figures and tables in the paper.

## Commands not found in the run records

1. `build_master_table.py` and `build_ld_eur_gctb.py`: exact argv (usage lines only); the 6-way split of the GCTB fetch.
2. `build_ld_aou.py`: full EAS command lines (only `--skip-a` and the checkpoint path recorded); `--out` for EUR.
3. How `all_targets.tsv` was made; the sync to `~/work/`. (`loci_250_pgc3scz.tsv` is PGC3-SCZ's locus table;
   the coordinates the fit needs are in `inputs/loci_250_coords.tsv`.)
4. `run_finemap.py` / `export_results.py` with `--coverage 0.95` (only the resulting argvs recorded).
5. `harmonize.py` `--loci-file`/`--out` and `susiex_eur_panel.py` arguments for the AoU-EUR SuSiEx run.
6. `lambda_merged.py`, `check_exports.py` and `figures/fetch_gene_track.py` invocations.
7. The `results/benchmark/` inputs read by step 8 (`export_published_pips.py`, `build_proxy_map.py`), which predate v2.
