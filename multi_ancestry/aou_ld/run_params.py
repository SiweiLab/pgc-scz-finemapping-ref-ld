#!/usr/bin/env python3
"""The one place the published analysis parameters live.

Standing rule added 2026-09-09 after the coverage-0.95 regression: no analysis
parameter may be a hardcoded literal in a sweep or rerun script. Import from here,
so a rerun cannot silently disagree with the original run.

Every value below is the value the ORIGINAL published run used, cross-checked
against the code in ../ORIGINAL_RUN_SPEC.md. Changing one changes what the
manuscript's numbers mean -- do not edit to make a rerun match; fix the rerun.
"""

# Credible-set coverage. The original ran 0.99 throughout (batch_perancestry_qcexcl.py:33).
# Note SuSiEx's own default is 0.95 and MultiSuSiE's wrapper default is 0.95, so this MUST
# be passed explicitly to both -- omitting it is what caused the 2026-09-08 regression.
COVERAGE = "0.99"

# Purity floor. SuSiEx calls it --min_purity, the SuSiE family min_abs_corr. 0.5 in all
# arms. SuSiEx's own default happens to be 0.5, but pass it explicitly anyway: a default
# that silently matches is indistinguishable from one that silently does not.
MIN_PURITY = "0.5"

# Maximum causal signals, AS RUN per arm -- each tool at its own published default.
# These are deliberately not uniform; see ORIGINAL_RUN_SPEC.md sec 3 item 1, which
# flags that the Methods text wrongly describes all three per-ancestry arms as L=10.
L_MULTISUSIE = 10
L_MESUSIE = 10
L_SUSIE_RSS = 10
N_SIG_SUSIEX = 5          # SuSiEx's --n_sig, the same quantity as L
KMAX_FINEMAP = 5          # FINEMAP's --n-causal-snps

# MultiSuSiE's cross-ancestry prior effect correlation. A fixed scalar, never estimated.
# Undeclared in the Methods; every published number used this value.
RHO_MULTISUSIE = 0.75

# SuSiEx marginal p-value gate. 1 disables it, leaving purity as the only CS filter.
PVAL_THRESH_SUSIEX = "1"

# Kai Yuan's published run passed study-wide N_eff fixed across all 250 loci (his
# submission script, shared 2026-09-11). Our primary run computes N per locus from each
# locus's median Weight -- an undocumented second difference from the published analysis.
N_GWAS_PUBLISHED_EUR = 58749
N_GWAS_PUBLISHED_EAS = 14425

# Merged-LD arms, fixed (Qin, 2026-09-23). Previously each merged locus took
# `n_samples = round(median(Weight))` of its OWN core-sumstats slice
# (finemap/build_finemap_input_gctb.py:67), so N varied locus to locus while the
# per-ancestry arms used a fixed value -- an inconsistency, not a design choice.
#
# Verified on the laptop 2026-09-23 against the public PGC3 files, whole-file median NEFF:
#     european  58,749.1   <- reproduces N_GWAS_PUBLISHED_EUR exactly
#     asian     14,424.8   <- reproduces N_GWAS_PUBLISHED_EAS exactly
#     core      73,173.9
# So Kai's two constants are simply each file's median NEFF, and the meta's is additive:
# 58,749.1 + 14,424.8 = 73,173.9. The merged arms therefore use the same convention with
# no special-casing. This is NOT the total 161,405 (= cases + controls, N_GWAS_MAWSON_TOTAL),
# which is a different quantity.
N_GWAS_MERGED = 73174

# FINEMAP credible-set probability. The merged arm passed this; the single-ancestry
# diagnostic arms did not (ORIGINAL_RUN_SPEC.md sec 3 item 4).
PROB_CRED_SET_FINEMAP = "0.99"

# The single-ancestry FINEMAP diagnostic arm (EUR-FINEMAP, EAS-FINEMAP) never passed
# --prob-cred-set at all, as run -- confirmed against `finemap --help`, whose own
# default is 0.95, not the 0.99 the Methods claims applies "throughout". Recorded here
# as an explicit constant (Task B8, 2026-09-11) so the omission is a documented fact
# about the original run, not an absent flag a rerun could silently "fix". Do not pass
# --prob-cred-set for this arm; this is not a param to set, it's a record that none was set.
PROB_CRED_SET_FINEMAP_SINGLEANC_OMITTED = True

# Per-subprocess BLAS/OpenMP pinning, required whenever loci run concurrently.
THREAD_LIMIT_ENV = {
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
}


# ---------------------------------------------------------------------------
# Sweep grids (Qin, 2026-09-09). These describe RERUNS, not the published run.
# ---------------------------------------------------------------------------

# L / n_sig / k_max values to sweep. 1 and 3 probe the low end, 5 is FINEMAP's
# published k_max and SuSiEx's mult-step starting point, 10 is the SuSiE family's
# published L -- so the grid contains both conventions' published values.
L_SWEEP_VALUES = [1, 3, 5, 10]

# SuSiEx ignores --n_sig entirely when mult-step is on (src/model.cpp: the mult-step
# branch calls init(5) and never reads par.n_sig), so the sweep must disable it or it
# sweeps nothing. This is a DEPARTURE from the published run, which used True --
# label these arms distinctly and keep the published-settings arm as a reference.
MULT_STEP_SUSIEX_SWEEP = "False"
MULT_STEP_SUSIEX_PUBLISHED = "True"

# Task B10 (2026-09-12): a purity x coverage diagnostic grid at published settings
# otherwise. --min_purity=0 makes SuSiEx emit every component with its own CS_PURITY
# (not just purity-survivors), so the purity axis is applied post hoc from these two
# coverage runs rather than needing a 4-cell grid on the command line. NOT the
# published MIN_PURITY (0.5) -- these are diagnostic arms, own --out_dir/label, never
# used in place of MIN_PURITY for a primary arm.
MIN_PURITY_ZERO = "0"
COVERAGE_B10_GRID = ["0.99", "0.95"]

# RSparsePro (Task B8 Task 2, 2026-09-11). K deliberately equals L_SUSIE_RSS/L_MULTISUSIE/
# L_MESUSIE=10 -- RSparsePro's K is the same "max number of single effects" quantity as
# SuSiE's L (see rsparsepro_ld.py's RSparsePro.__init__: self.k mirrors susieR's L exactly),
# so this records that equivalence rather than letting it look like an independent choice.
# cthres/eincre are RSparsePro's own tool defaults, passed explicitly for the same reason
# MIN_PURITY is passed explicitly to SuSiEx/MultiSuSiE above: a default that silently
# matches is indistinguishable from one that silently does not.
K_RSPARSEPRO = 10
CTHRES_RSPARSEPRO = 0.95
EINCRE_RSPARSEPRO = 1.5

# rho is MultiSuSiE-only. SuSiEx has no cross-ancestry effect-correlation parameter at
# all (its full option list carries none), and MESuSiE models sharing as a categorical
# per-credible-set configuration (EUR / EAS / shared) rather than a correlation, so
# there is nothing to sweep in either. rho=1 is not usable: with two ancestries it makes
# the prior covariance rank-1 and MultiSuSiE raises "rho must be positive definite".
RHO_SWEEP_VALUES = [0, 0.75, 0.99]
RHO_SWEEP_ARMS = ["MultiSuSiE"]


# ---------------------------------------------------------------------------
# core_sumstats provenance (Task B15, 2026-09-13). The original merged-LD run
# read `--core-sumstats-dir {FULL250}/core_sumstats` (batch_merged_ld_qcexcl.py:58),
# a file never backed up. The rebuild substituted core_sumstats_aligned_raw instead
# (batch_build_merged_ld_b1.py:20) -- confirmed pre-intersected, not raw
# (task_b15_out/phase1_report.md) -- which is what dropped 5/250 loci by one
# variant each (ORIGINAL_RUN_SPEC.md sec 3 item 7). This is the canonical
# replacement: never hardcode either path literal in a driver script again.
# ---------------------------------------------------------------------------

# The PGC3-SCZ trans-ancestry ("core") meta-analysis sumstats -- the actual source
# core_sumstats was always built from (confirmed against archive/eur_eas_merge_notes.md,
# which already used this exact file with prep_pgc_locus.py for this exact purpose).
# md5 1d8a9885cbd21b4a2f3dc50e7f9ed693, 7,585,077 data rows. Backed up at
# gs://dataproc-staging-wb-glittery-spinach-6060/head_node_backup_2026-08-07/work/pgc3_scz_raw/
# (md5-verified identical) before this task, so not re-backed-up by B15.
PGC3_CORE_META_SUMSTATS = "~/work/pgc3_scz_raw/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv.gz"

# core_sumstats, rebuilt from PGC3_CORE_META_SUMSTATS via data_prep/prep_pgc_locus.py
# (Task B15). This is the canonical --core-sumstats-dir for merge_eur_eas_ld.py going
# forward -- NOT core_sumstats_aligned_raw (pre-intersected leftover, kept only as
# Phase 1 evidence, never delete) and NOT the original {FULL250}/core_sumstats (lost).
# Backed up: gs://dataproc-staging-wb-glittery-spinach-6060/task_b15_backup_2026-09-13/core_sumstats_rebuilt/
CORE_SUMSTATS_REBUILT_DIR = "~/work/full250_qcexcl_derived/core_sumstats_rebuilt"

# finemap_in built from CORE_SUMSTATS_REBUILT_DIR (Task B15's merge_eur_eas_ld.py rerun),
# 215,341 variants total -- the corrected count, vs 215,336 in the pre-B15 finemap_in
# (STALE_finemap_in, built from the pre-intersected core_sumstats_aligned_raw). This is
# the canonical --in-files/FINEMAP_IN source for SuSiE-merged and FINEMAP-merged going
# forward (Task B18, 2026-09-14) -- the pre-B15 sweep (task_b1_out/l_sweep_*) used the
# stale tree and needs the rerun this task performs.
# Backed up: gs://dataproc-staging-wb-glittery-spinach-6060/task_b15_backup_2026-09-13/finemap_in_rebuilt/
FINEMAP_IN_REBUILT_DIR = "~/work/full250_qcexcl_derived/finemap_in_rebuilt"


# ---------------------------------------------------------------------------
# Mawson et al. 2026 matched-settings sensitivity arm (Task B20, 2026-09-15).
# Isolates the LD-source question their paper actually asks by matching everything
# else. See ESTABLISHED_FACTS.md sec 16.
# ---------------------------------------------------------------------------

# Mawson's own Methods: "the sample size was specified as n=161405 (the GWAS sample
# size)" -- the raw case+control total. NOT our canonical convention: our merged-LD
# n_samples (read by build_finemap_input_gctb.py's median-Weight default, and by both
# SuSiE-merged and FINEMAP-merged off the same master$n_samples) is 73,174, the SUM of
# N_GWAS_PUBLISHED_EUR + N_GWAS_PUBLISHED_EAS above -- Kai Yuan's own published SuSiEx
# --n_gwas arguments, an effective-N convention, not a mistake. This constant measures
# how much Mawson's total-N convention moves the answer; adopting it is not implied.
N_GWAS_MAWSON_TOTAL = 161405

# Mawson's coverage is STATED, not inferred (corrected 2026-09-16): their Results report
# "95% credible sets" as primary throughout, with 90/85/80% as deliberately relaxed
# variants. An earlier version of this comment called it unstated and reasoned from the
# susie_rss/FINEMAP defaults; the conclusion was right, the provenance was wrong.
COVERAGE_MAWSON = "0.95"

# Mawson's purity gate, confirmed 2026-09-16 from their Supplementary Note
# (Zotero KIHPS96Z/media-9.docx): "Credible sets that do not meet a certain purity
# threshold (default r=0.5) ... are consequently discarded." They never mention changing
# it, and the Note contains no min_abs_corr or coverage string at all. So their gate
# equals MIN_PURITY above -- purity is MATCHED to our canonical p05rebuilt run, not a
# difference. Note also susieR's own susie_get_cs default is 0.5; the 0.1 that caused
# the sec 3 item 6 defect was our run_susie_merged.R script default, never susieR's.
MIN_PURITY_MAWSON = "0.5"

# Mawson: "maximum number of causal variants was set to 5" for both susieR and FINEMAP.
# Matches KMAX_FINEMAP already; SuSiE-merged's own primary/canonical run uses L=10
# (archive/run_susie_merged.R), so L=5 is the one new value for that arm here.
L_MAWSON = 5


# ---------------------------------------------------------------------------
# Canonical per-ancestry inputs, one named home (Task B23, 2026-09-16). Before this,
# only the meta side had a constant (CORE_SUMSTATS_REBUILT_DIR, FINEMAP_IN_REBUILT_DIR
# above) -- EUR/EAS paths were picked by hand per driver, and the 250 loci were split
# across full250_qcexcl (11 loci, survived the 2026-09-08 cluster deletion) and
# full250_qcexcl_derived (that deletion's reconstruction). All five constants below were
# individually verified to cover all 250 loci with zero missing/zero extra, and (except
# LD_EUR_DIR, which has no independent pre-deletion copy to check against) verified
# byte-identical, decompressed, against the 11-locus survivor copies where one exists.
# Full verification: aou_ld/task_b23_out/{MANIFEST.csv, canonical_inputs_report.md}.
# Backed up: gs://dataproc-staging-wb-glittery-spinach-6060/task_b23_backup_2026-09-16/canonical_inputs_250_qcexcl/
# ---------------------------------------------------------------------------

# Root of the consolidated tree. Every constant below is a subdirectory of this one.
CANONICAL_INPUTS_DIR = "~/work/canonical_inputs_250_qcexcl"

# Per-ancestry GWAS sumstats (raw, not LD-intersected) -- read only for NEFF merge
# weights by the merged-LD build (ORIGINAL_RUN_SPEC.md sec 3 item 7), never for variant
# intersection. Source: full250_qcexcl_derived/{eur,eas}_sumstats, which already covers
# all 250 loci in one tree (not just the 239 "derived" ones -- confirmed byte-identical,
# decompressed, to full250_qcexcl's 11 independent survivor copies at every overlap).
# unified_sumstats/EUR and /EAS (read by merge_eur_eas_ld.py's actual invocation,
# task_b15_out/batch_build_merged_ld_rebuilt_b15.py) are themselves plain symlinks to
# these same two directories -- confirmed via `file`, not just by name.
EUR_SUMSTATS_DIR = f"{CANONICAL_INPUTS_DIR}/eur_sumstats"
EAS_SUMSTATS_DIR = f"{CANONICAL_INPUTS_DIR}/eas_sumstats"

# Trans-ancestry meta sumstats (raw, not LD-intersected). Content-identical to
# CORE_SUMSTATS_REBUILT_DIR above (same source file, copied into the new tree for the
# one-stop convenience the other four constants here provide) -- NOT
# core_sumstats_aligned or core_sumstats_aligned_raw (both STALE_-prefixed, both
# confirmed content-identical to EACH OTHER and both the flawed pre-B15 file -- Task
# B23's Q2). CORE_SUMSTATS_REBUILT_DIR is left unchanged; nothing that already imports
# it needs to change.
CORE_SUMSTATS_DIR = f"{CANONICAL_INPUTS_DIR}/core_sumstats"

# Per-ancestry LD, symlinked (not copied -- these are large). EUR:
# full250_qcexcl_derived/eur_ukb_ld_named, the only EUR LD tree on this cluster; no
# independent 11-locus survivor copy exists to check it against (the only surviving
# 11-locus EUR-side artifact, full250_qcexcl/merged_ld, is a LATER pipeline stage --
# already EUR+EAS merged, allele-flip columns keyed to EAS -- comparing it directly gave
# a spurious mismatch from comparing the wrong pair, not a real one; see the brief).
LD_EUR_DIR = f"{CANONICAL_INPUTS_DIR}/ld_eur"

# EAS: full250_qcexcl_derived/eas_npz_out_merged -- already the 11+239 consolidation,
# confirmed byte-identical (R matrix, np.array_equal) to full250_qcexcl/eas_npz_out's
# independent 11-locus survivor subdirectories at every one of the 11.
LD_EAS_DIR = f"{CANONICAL_INPUTS_DIR}/ld_eas"


# ---------------------------------------------------------------------------
# Version pins and previously-unrecorded defaults (Task B29, 2026-09-23).
# Every value below was read from the tool itself -- the binary's own output,
# formals() on the installed package, or the source at the pinned commit -- never
# from documentation, because B12, B28 and B29 all found docs disagreeing with code.
# These RECORD what the original runs used. Standing rule unchanged: never edit one
# to make a rerun agree; pass it explicitly in the rerun instead.
# Full table and evidence: ../VERSIONS.md, aou_ld/task_b29_out/.
# ---------------------------------------------------------------------------

# FINEMAP. Every run in this project used the official v1.4.2 x86_64 build (md5
# b5741e4374396696f511f2126e72efba, byte-identical to christianbenner.com's tarball).
# v1.4.2 is also the latest release, so the rerun uses the same binary.
# Two traps, both verified in B29:
#   1. v1.4.2 writes "Welcome to FINEMAP v1.4.1" into its .log_sss FILE (a stale string
#      in the binary) while printing v1.4.2 on stdout. B28 read the .log_sss banners and
#      wrongly concluded the runs used v1.4.1. Trust stdout or `--help`, never .log_sss.
#   2. --help in BOTH v1.4.1 and v1.4.2 says "--n-conv-sss ... Default is 1000", but
#      both binaries actually run with 100 when the flag is omitted (logged SETTINGS block,
#      FINEMAP's own example data). B28's "1000, 10x stricter" was the help text, not the
#      code. Every run in this project therefore used n-conv-sss = 100.
# FINEMAP has no seed option and is non-deterministic: two identical runs on its 55-SNP
# example give different .snp output in both versions (B18's 0.170 floor is inherent).
FINEMAP_VERSION = "1.4.2"
FINEMAP_BINARY_MD5 = "b5741e4374396696f511f2126e72efba"
# Effective defaults as run (never passed by any driver). Pass them explicitly from now on,
# so a future binary whose defaults differ cannot silently change the analysis.
N_CONV_SSS_FINEMAP = 100
PROB_CONV_SSS_TOL_FINEMAP = 0.001
N_ITER_FINEMAP = 100000
CORR_CONFIG_FINEMAP = 0.95
PRIOR_STD_FINEMAP = 0.05
N_CONFIGS_TOP_FINEMAP = 50000

# susieR. CRAN 0.14.2 -- the latest RELEASE, and what every run used (B21 packageVersion).
# GitHub master is 0.16.6 (8e56a8e, 821 commits past v0.14.2) with a REWRITTEN susie_rss:
# prior_variance=50 -> scaled_prior_variance=0.2, z_ld_weight and check_prior removed (our
# drivers' check_prior=FALSE fallback would stop working), new R_mismatch /
# unmappable_effects / greedy-L. Not a drop-in upgrade; the rerun stays on 0.14.2.
SUSIER_VERSION = "0.14.2"
# 0.14.2 susie_rss(z, R, n, ...) defaults that no driver passes. With n given, susie_rss
# routes through susie_suff_stat, so its defaults apply.
SUSIER_RSS_DEFAULTS_AS_RUN = dict(
    z_ld_weight=0, estimate_residual_variance=False, check_prior=True,
    scaled_prior_variance=0.2, estimate_prior_variance=True, estimate_prior_method="optim",
    max_iter=100, tol=1e-3, null_weight=0, standardize=True, refine=False, n_purity=100,
)
# The native single-ancestry SuSiE-RSS arms (susie/run_susie_gctb_aou.R) HARDCODE
# coverage=0.95, min_abs_corr=0.1, L=10 -- not COVERAGE / MIN_PURITY above. Recorded as run.
SINGLEANC_SUSIE_COVERAGE_AS_RUN = "0.95"
SINGLEANC_SUSIE_MIN_ABS_CORR_AS_RUN = "0.1"

# MultiSuSiE. pip install git+https://github.com/jordanero/MultiSuSiE, python 3.11 env.
# No tags; last upstream commit 10351eb (2026-02-15) predates every install, so every run
# used 10351eb and it is still the latest.
MULTISUSIE_COMMIT = "10351eb"
# multisusie_rss() library defaults. NOTE min_abs_corr's LIBRARY default is 0 and coverage
# 0.95; the 0.5 / 0.99 the runs used come only from our wrapper (run_multisusie.py) and
# batch_perancestry_qcexcl.py. Unpassed defaults as run:
MULTISUSIE_RSS_DEFAULTS_AS_RUN = dict(
    rho=0.75, scaled_prior_variance=0.2, standardize=False, pop_spec_standardization=True,
    estimate_residual_variance=True, estimate_prior_variance=True,
    estimate_prior_method="early_EM", pop_spec_effect_priors=True,
    iter_before_zeroing_effects=5, prior_tol=1e-9, max_iter=100, tol=1e-3,
    float_type="float32", low_memory_mode=False, single_population_mac_thresh=20,
    multi_population_maf_thresh=0,
)

# MESuSiE. remotes::install_github("borangao/MESuSiE"). No tags; last upstream commit
# 5a3dd32 (2025-08-02) predates every install -> every run used it, still latest.
# Upstream STILL hardcodes coverage=0.95 in meSuSie_get_cs at 5a3dd32, so
# run_mesusie.R's get_cs_fixed remains necessary.
MESUSIE_COMMIT = "5a3dd32"
# meSuSie_core(R_mat_list, summary_stat_list, L=10) -- only L is passed. Unpassed:
MESUSIE_CORE_DEFAULTS_AS_RUN = dict(
    residual_variance=None, prior_weights=None, ancestry_weight=None, optim_method="optim",
    estimate_residual_variance=False, max_iter=100, cor_method="min.abs.corr",
    cor_threshold=0.5,
)

# SuSiEx. git clone getian107/SuSiEx, built from src. v1.1.2-3-g1db2f58 (ORIGINAL_RUN_SPEC
# sec 4a); 1db2f58 is still upstream HEAD. Compiled-in defaults (src/data.hpp) for the
# flags batch_perancestry_qcexcl.py does NOT pass (it passes n_sig, level, pval_thresh,
# mult-step, keep-ambig):
SUSIEX_COMMIT = "1db2f58"
SUSIEX_DEFAULTS_AS_RUN = dict(maf=0.005, min_purity=0.5, tol=1e-4, max_iter=100, threads=1)

# RSparsePro. git clone zhwm/RSparsePro_LD (2026-09-08), no install step. No tags; HEAD
# 5990896 (2026-06-23) predates the clone -> every run used it, still latest. Its only
# 2026 change to src/rsparsepro_ld.py (34bef1b) altered no argparse default.
RSPARSEPRO_COMMIT = "5990896"
RSPARSEPRO_DEFAULTS_AS_RUN = dict(maxite=100, eps=1e-5, ubound=100000, minldthres=0.7,
                                  maxldthres=0.2, varemax=100.0, varemin=1e-3)

# CARMA (Task B29 Job 1 pilot only, 2 loci). Built via a user-writable conda env,
# ~/envs/carma (R 4.5.3, RcppGSL 0.3.14, GSL 2.7) -- NOT Dataproc's R 4.3.3. HEAD 856f0fb.
CARMA_COMMIT = "856f0fb"
CARMA_R_VERSION = "4.5.3"
# lambda.list is REQUIRED despite formals() defaulting it to NULL (the body calls
# log(lambda.list[[i]]) unconditionally). 1/sqrt(p) is the package example's value.
CARMA_LAMBDA_RULE = "1/sqrt(p)"
# CARMA's own formals as installed (agree with the man page on every value B29's brief lists):
CARMA_DEFAULTS_AS_RUN = dict(
    num_causal=10, rho_index=0.99, BF_index=10, effect_size_prior="Spike-slab",
    outlier_switch=True, outlier_BF_index=1 / 3.2, Max_Model_Dim=2e5, all_iter=3,
    all_inner_iter=10, input_alpha=0, epsilon_threshold=1e-5, y_var=1, tau=0.04,
    EM_dist="Logistic", prior_prob_computation="Logistic",
)



# ---------------------------------------------------------------------------
# CARMA (Task B29, 2026-09-23). CARMA's own installed defaults, passed explicitly.
# `rho.index` IS coverage -- the paper defines a credible set as "the smallest subset of
# correlated variants (with correlation > r) that has a probability rho or greater of
# containing at least one causal variant", i.e. SuSiE's sense. Its default already equals
# COVERAGE, but a default that silently matches is indistinguishable from one that silently
# does not, so it is passed. CARMA exposes NO purity/correlation argument.
# `lambda.list` is not a constant: run_carma_pilot_b29.R computes 1/sqrt(p) per locus, the
# value CARMA's own example uses, and it is required despite formals() defaulting it to NULL.
# ---------------------------------------------------------------------------
NUM_CAUSAL_CARMA = 10
BF_INDEX_CARMA = 10
OUTLIER_SWITCH_CARMA = True
# RNG seed (Qin, 2026-09-30, B55). CARMA's model search is stochastic, and run_carma_pilot_b29.R never set a seed, so
# every earlier CARMA run (v1 B30/B34, and v2 B55's first noambig pass) varies from run to run: two identical unseeded
# runs at chr5_86721950 differed by max |dPIP| 0.023. Passed as the runner's argument 7.
SEED_CARMA = 42



# ---------------------------------------------------------------------------
# 2026-09 primary run: values that DEPART from the as-run constants above (Task B30).
# Signed off in PRIMARY_RUN_PLAN.md sec 4 (Qin, 2026-09-23/24). The as-run constants stay
# untouched as the record of the original run; the primary-run driver imports these.
# ---------------------------------------------------------------------------

# RSparsePro --cthres is its coverage (PRIMARY_RUN_PLAN sec 4, verified from
# rsparsepro_ld.py's get_effect). Aligned with COVERAGE; the original ran its 0.95 default.
CTHRES_RSPARSEPRO_PRIMARY = 0.99

# FINEMAP --n-conv-sss: 1000, passed explicitly (PRIMARY_RUN_PLAN sec 4, footnote on the
# departure). Every original run used 100 (N_CONV_SSS_FINEMAP, B29).
N_CONV_SSS_FINEMAP_PRIMARY = 1000

# SuSiEx --n_gwas: the FIXED published N (Qin, 2026-09-24) -- N_GWAS_PUBLISHED_EUR/EAS above,
# as in Kai Yuan's published run. The original primary run passed each locus's median
# `Weight` instead (batch_perancestry_qcexcl.py:87), which is 57,602-58,749 for EUR.
SUSIEX_N_GWAS_RULE_PRIMARY = "fixed_published"
