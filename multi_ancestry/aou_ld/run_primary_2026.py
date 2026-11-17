#!/usr/bin/env python3
"""Single driver for the 2026-09 primary run: all ten arms, one parameter source.

**STATUS: verified in B30** (dry run checked against every runner, guards fired, ten-arm
smoke test at ATP2A2 -- aou_ld/PROGRESS_LOG.md). Drafted on the laptop; B30 fixed the
plumbing defects the dry run exposed (RSparsePro input format, FINEMAP masters pointing at
an earlier cluster, EUR-SuSiE-RSS's shared output file, wrong output
targets). The EUR arms' FINEMAP-format .z/.ld are rebuilt per locus from ld_per_anc_eur +
eur_sumstats rather than stored as a text copy of the LD and moved every tool into one conda env (envs/ilae4_finemap.yml).

Why this exists. The primary run was previously spread across `batch_perancestry_qcexcl.py`,
`batch_merged_ld_qcexcl.py`, `batch_native_singleanc_250_qcexcl.py` and the sweep drivers,
each passing analysis parameters independently. So: **every analysis parameter here is
imported from `run_params`. There are no literals.** `assert_no_literals()` fails the run if
that stops being true.

Design rules, all enforced below:
  1. One arm = one entry in `ARMS`: a `prep` step that builds any intermediate input the
     runner needs, a `cmd` builder, and the file whose existence means "done".
  2. Idempotent: an arm whose output file exists is skipped, so a killed run resumes.
  3. A run only counts as ok if its runner exits 0 AND the output file then exists.
  4. Per-locus work is a subprocess with BLAS/OpenMP pinned to 1 thread
     (`run_params.THREAD_LIMIT_ENV`); concurrency is across loci, never inside a tool.
  5. Failures are recorded per (arm, locus) and never abort the batch. Every run's
     stdout+stderr goes to {WORK}/{arm}/_logs/{locus}.log.
  6. Each arm is rsynced to gs:// (default: the manuscript bucket's
     fine_mapping_runs/primary_2026-09/<arm>/, `_inputs/` excluded) the moment its last
     locus finishes. --no-backup turns this off for tests.

Usage:
  python3 run_primary_2026.py --stage-inputs --loci chr12_110129434_111489434   # fetch inputs
  python3 run_primary_2026.py --dry-run --loci chr12_110129434_111489434        # review argv
  python3 run_primary_2026.py --loci chr12_110129434_111489434 --no-backup      # smoke test
  python3 run_primary_2026.py --jobs 14     # full run; backs up to fine_mapping_runs/primary_2026-09/

L sweep (B31): `--L n` overrides the arm's L-equivalent (table SWEEP_L below) for the selected
arms, writes to WORK/<arm>_L{n} (never WORK/<arm>), keeps state in WORK/state_L{n}.json, and
runs SuSiEx at MULT_STEP_SUSIEX_SWEEP (mult-step ignores --n_sig). Only SWEEP_L arms are
allowed, and a point the primary run already holds (PRIMARY_L) is refused, not recomputed.
assert_no_literals() checks the override actually reached every selected arm's argv.
  python3 run_primary_2026.py --L 3 --arms MultiSuSiE,CARMA --dry-run --loci chr12_110129434_111489434
"""
import argparse
import gzip
import concurrent.futures as cf
import csv
import json
import os
WORKSPACE_BUCKET = os.environ.get("WORKSPACE_BUCKET", "gs://YOUR-WORKSPACE-BUCKET")  # AoU Workbench sets this
import subprocess
import sys
import threading
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_params as P

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MA_DIR = os.path.dirname(SCRIPT_DIR)
FM_DIR = os.path.dirname(MA_DIR)          # Post-hoc_analyses/fine-mapping (holds susie/)
HOME = os.path.expanduser("~")
# Pipeline v2 (step 7): PRIMARY_V2_INPUTS = one harmonize_v2.py tree (~/work/v2_harmonized/{primary,noambig}).
# It supplies every arm's inputs and replaces the v1 trees below, which are deleted after B55.
V2 = os.path.expanduser(os.environ["PRIMARY_V2_INPUTS"]) if os.environ.get("PRIMARY_V2_INPUTS") else None
V2_VARIANT = os.path.basename(os.path.normpath(V2)) if V2 else None
WORK = os.path.expanduser(os.environ.get(
    "PRIMARY_WORK", f"~/work/primary_v2_{V2_VARIANT}" if V2 else "~/work/primary_2026"))

# --- Inputs: a local copy of the new bucket's ld_matrices/ layout, filled by --stage-inputs.
# These are inputs, not parameters.
INPUTS = os.path.expanduser(os.environ.get("PRIMARY_INPUTS", "~/work/primary_2026_inputs"))
# PRIMARY_UNION: B33 points this at ld_eur_eas_union_noambig/ (palindromic variants dropped)
UNION = os.path.expanduser(os.environ.get("PRIMARY_UNION", f"{INPUTS}/ld_eur_eas_union"))
# PRIMARY_INTERSECT / PRIMARY_MERGED: B34 points these at the *_noambig/ subsets
INTERSECT = os.path.expanduser(os.environ.get("PRIMARY_INTERSECT",
                                              f"{INPUTS}/ld_eur_eas_intersect"))  # MESuSiE
MERGED = os.path.expanduser(os.environ.get("PRIMARY_MERGED",
                                           f"{INPUTS}/ld_eur_eas_merged"))  # merged arms
# PRIMARY_DROP_AMBIGUOUS=1 (B34 §2b): the EUR-only arms build their .z/.ld with the builder's
# --drop-ambiguous, and EUR-SuSiE-RSS reads the filtered snps table the builder writes next to
# them instead of {EUR_LD}/{l}_snps.txt.gz. The inputs above are selected by path, not by this.
DROP_AMBIGUOUS = os.environ.get("PRIMARY_DROP_AMBIGUOUS", "") == "1"
EUR_LD = f"{INPUTS}/ld_per_anc_eur"           # {l}_R.npz, {l}_snps.txt.gz (EUR arms)
EUR_SS = f"{INPUTS}/eur_sumstats"             # {l}_sumstats.tsv.gz        (EUR arms)
# build_finemap_input_gctb.py wants {l}_ukb_R.npz / _ukb_snps.txt.gz; ld_per_anc_eur holds the
# same bytes as {l}_R.npz / _snps.txt.gz (B28), so the prep step symlinks them under these names.
EUR_LD_UKBNAMES = f"{INPUTS}/ld_per_anc_eur_ukbnames"
# MultiSuSiE reads the npz tree named here. v2: censored on the both-LD domain (2026-09-30,
# B54): a variant missing z in one ancestry gets z = 0 and a zero R row/col there, per the
# library's z-path contract, instead of v1's fitted-as-observed z = 0.
MULTISUSIE_NPZ = UNION
EUR_V2 = None
if V2:
    UNION, INTERSECT, MERGED = f"{V2}/union", f"{V2}/intersect", f"{V2}/merged"
    MULTISUSIE_NPZ, EUR_V2 = f"{V2}/union_censored_bothld", f"{V2}/eur"

BUCKET = f"{WORKSPACE_BUCKET}/manuscript_final_2026-09"
BUCKET_LD = f"{BUCKET}/ld_matrices"
BUCKET_EUR_SS = f"{BUCKET}/PGC_SCZ_sumstats/eur_sumstats"
FINEMAP_INPUT_BUILDER = os.path.join(FM_DIR, "finemap", "build_finemap_input_gctb.py")

# --- Tools: one conda env + three non-conda tools, built by envs/setup_ilae4_env.sh.
ENV = os.path.expanduser(os.environ.get("ILAE4_ENV", "~/envs/ilae4_finemap"))
RSCRIPT = f"{ENV}/bin/Rscript"
PYTHON = f"{ENV}/bin/python"
TOOLS = os.path.expanduser(os.environ.get("TOOLS_DIR", "~/tools"))
SUSIEX_BIN = f"{TOOLS}/SuSiEx/bin/SuSiEx"
PLINK_BIN = f"{TOOLS}/SuSiEx/utilities/plink"
RSPARSEPRO_PY = f"{TOOLS}/RSparsePro_LD/src/rsparsepro_ld.py"
FINEMAP_BIN = f"{TOOLS}/finemap_v1.4.2/finemap_v1.4.2_x86_64"


# --- L sweep (B31). None = primary run. Set once in main() from --L.
L_OVERRIDE = None
# --- Coverage override (B69). None = run_params.COVERAGE (0.99, never edited). Set from --coverage.
# Only the arms in COVERAGE_ARMS take it; outputs go to WORK/<arm>[_L{n}]_cov<NNN>.
COV_OVERRIDE = None
COVERAGE_ARMS = ("SuSiE-merged", "FINEMAP-merged")


def cov():
    return COV_OVERRIDE if COV_OVERRIDE is not None else P.COVERAGE


def cov_tag():
    return "" if COV_OVERRIDE is None else "_cov" + COV_OVERRIDE.replace(".", "").ljust(3, "0")

# arm -> the run_params constant holding its L-equivalent (max number of single effects).
SWEEP_L = {"SuSiEx": "N_SIG_SUSIEX", "MultiSuSiE": "L_MULTISUSIE", "MESuSiE": "L_MESUSIE",
           "SuSiE-merged": "L_SUSIE_RSS", "FINEMAP-merged": "KMAX_FINEMAP",
           "RSparsePro": "K_RSPARSEPRO", "CARMA": "NUM_CAUSAL_CARMA"}


def L_of(arm):
    return L_OVERRIDE if L_OVERRIDE is not None else getattr(P, SWEEP_L[arm])


def primary_L(arm):
    """The L the primary run used for this arm, or None when no primary output is the same
    fit at any L (SuSiEx: its sweep runs mult-step off, the primary ran it on)."""
    return None if arm == "SuSiEx" else getattr(P, SWEEP_L[arm])


def out_dir(arm):
    return os.path.join(WORK, (arm if L_OVERRIDE is None else f"{arm}_L{L_OVERRIDE}") + cov_tag())


def prep_dir(arm):
    return os.path.join(out_dir(arm), "_inputs")


def _link(src, dst):
    if not os.path.lexists(dst):
        os.symlink(src, dst)


def _write_master(path, z, ld, out_prefix, n):
    """FINEMAP master. FINEMAP writes .snp/.config/.cred/.log wherever the master says, so
    each FINEMAP arm gets its own master pointing outputs into that arm's out dir. The
    backed-up masters hard-code old-cluster paths and are never reused."""
    with open(path, "w") as f:
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(f"{z};{ld};{out_prefix}.snp;{out_prefix}.config;{out_prefix}.cred;"
                f"{out_prefix}.log;{n}\n")


# ---------------------------------------------------------------------------
# Arm registry. prep(locus) builds intermediate inputs; cmd(locus, info) returns argv;
# out is the file whose existence means "done". No numeric literals.
# ---------------------------------------------------------------------------

def _susiex(locus, info):
    """Matches batch_perancestry_qcexcl.py:94-107 except two signed-off changes:
    --min_purity 0 (gated post hoc, PRIMARY_RUN_PLAN 4.3) and --n_gwas = fixed published N
    (SUSIEX_N_GWAS_RULE_PRIMARY; the original used each locus's median Weight).
    Flag spellings: --key=value; --mult-step with a HYPHEN; --min_purity/--pval_thresh
    with underscores."""
    d = f"{UNION}/{locus}_susiex"
    c, a, b = info["chrom"], info["start"], info["stop"]
    return [SUSIEX_BIN,
            f"--sst_file={d}/EUR.sst.txt,{d}/EAS.sst.txt",
            f"--n_gwas={P.N_GWAS_PUBLISHED_EUR},{P.N_GWAS_PUBLISHED_EAS}",
            f"--ld_file={d}/EUR,{d}/EAS",
            f"--out_dir={out_dir('SuSiEx')}", f"--out_name={locus}_susiex",
            f"--chr={c}", f"--bp={a},{b}",
            "--chr_col=1,1", "--snp_col=2,2", "--bp_col=3,3", "--a1_col=4,4",
            "--a2_col=5,5", "--eff_col=6,6", "--se_col=7,7", "--pval_col=9,9",
            f"--plink={PLINK_BIN}",
            f"--n_sig={L_of('SuSiEx')}",
            f"--level={P.COVERAGE}",
            f"--min_purity={P.MIN_PURITY_ZERO}",
            f"--pval_thresh={P.PVAL_THRESH_SUSIEX}",
            f"--mult-step={P.MULT_STEP_SUSIEX_PUBLISHED if L_OVERRIDE is None else P.MULT_STEP_SUSIEX_SWEEP}",
            "--keep-ambig=False"]


def _multisusie(locus, info):
    """Every analysis parameter passed (run_multisusie.py has no defaults any more),
    including --population-sizes: the fixed published N, as for SuSiEx (B30). The original
    run used the npz's per-locus median-Weight `n` (EUR 57,602-58,749)."""
    return [PYTHON, os.path.join(MA_DIR, "multisusie", "run_multisusie.py"), locus,
            "--harmonized-dir", MULTISUSIE_NPZ, "--out", out_dir("MultiSuSiE"),
            "--use-z", "--L", str(L_of("MultiSuSiE")), "--coverage", P.COVERAGE,
            "--min-abs-corr", P.MIN_PURITY, "--rho", str(P.RHO_MULTISUSIE),
            "--population-sizes",
            f"EUR={P.N_GWAS_PUBLISHED_EUR},EAS={P.N_GWAS_PUBLISHED_EAS}"]


def _mesusie(locus, info):
    """run_mesusie.R: locus, ancestries, harm_dir, out_dir, coverage, L. Intersection tree."""
    return [RSCRIPT, os.path.join(MA_DIR, "mesusie", "run_mesusie.R"),
            locus, "EUR,EAS", INTERSECT, out_dir("MESuSiE"),
            P.COVERAGE, str(L_of("MESuSiE"))]


def _merged_prep(arm):
    """SuSiE-merged reads {dir}/{l}.z/.ld/.master and takes N from the master's n_samples,
    so give it a per-arm dir of symlinked .z/.ld plus a master with N_GWAS_MERGED."""
    def prep(locus):
        d = prep_dir(arm)
        _link(f"{MERGED}/{locus}.z", f"{d}/{locus}.z")
        _link(f"{MERGED}/{locus}.ld", f"{d}/{locus}.ld")
        _write_master(f"{d}/{locus}.master", f"{d}/{locus}.z", f"{d}/{locus}.ld",
                      f"{d}/{locus}", P.N_GWAS_MERGED)
    return prep


def _susie_merged(locus, info):
    """archive/run_susie_merged.R: locus, finemap_dir, out_dir, coverage, min_abs_corr, L.
    The fifth argument is what the original run omitted (-> min_abs_corr 0.1); the sixth was
    a hardcoded L = 10 until B31."""
    return [RSCRIPT, os.path.join(MA_DIR, "archive", "run_susie_merged.R"),
            locus, prep_dir("SuSiE-merged"), out_dir("SuSiE-merged"),
            cov(), P.MIN_PURITY, str(L_of("SuSiE-merged"))]


def _finemap_prep(arm, zld_dir, n):
    def prep(locus):
        o = out_dir(arm)
        _write_master(f"{o}/{locus}.master", f"{zld_dir}/{locus}.z", f"{zld_dir}/{locus}.ld",
                      f"{o}/{locus}", n)
    return prep


def _build_eur_zld(arm, locus):
    """EUR-only arms read FINEMAP-format {l}.z/{l}.ld. They are rebuilt here from
    ld_per_anc_eur + eur_sumstats with the same builder B8 used, instead of storing a text
    copy of the EUR LD in the bucket. B30 checked all 500 rebuilt files byte-identical (md5)
    to B8's. Each arm builds into its own _inputs/ so two arms never write one file."""
    d = prep_dir(arm)
    if os.path.exists(f"{d}/{locus}.z") and os.path.exists(f"{d}/{locus}.ld"):
        return
    if EUR_V2:  # v2: harmonize_v2.py's eur/ tree, joined by vkey and signed by v2_io
        _link(f"{EUR_V2}/{locus}.z", f"{d}/{locus}.z")
        _link(f"{EUR_V2}/{locus}.ld", f"{d}/{locus}.ld")
        _write_master(f"{d}/{locus}.master", f"{d}/{locus}.z", f"{d}/{locus}.ld", f"{d}/{locus}",
                      P.N_GWAS_PUBLISHED_EUR)
        return
    os.makedirs(EUR_LD_UKBNAMES, exist_ok=True)
    for src, dst in ((f"{EUR_LD}/{locus}_R.npz", f"{EUR_LD_UKBNAMES}/{locus}_ukb_R.npz"),
                     (f"{EUR_LD}/{locus}_snps.txt.gz", f"{EUR_LD_UKBNAMES}/{locus}_ukb_snps.txt.gz")):
        try:
            _link(src, dst)
        except FileExistsError:
            pass                                   # the other EUR arm linked it first
    subprocess.run([PYTHON, FINEMAP_INPUT_BUILDER, locus, "--ld-dir", EUR_LD_UKBNAMES,
                    "--sumstats-dir", EUR_SS, "--out", d, "--n", str(P.N_GWAS_PUBLISHED_EUR)]
                   + (["--drop-ambiguous"] if DROP_AMBIGUOUS else []),
                   check=True, capture_output=True, text=True,
                   env={**os.environ, **P.THREAD_LIMIT_ENV})


def _eur_finemap_prep(locus):
    _build_eur_zld("EUR-FINEMAP", locus)
    _finemap_prep("EUR-FINEMAP", prep_dir("EUR-FINEMAP"), P.N_GWAS_PUBLISHED_EUR)(locus)


def _finemap(arm):
    """FINEMAP v1.4.2, --n-conv-sss 1000 passed explicitly (N_CONV_SSS_FINEMAP_PRIMARY; every
    original run used 100). No --seed exists."""
    def cmd(locus, info):
        return [FINEMAP_BIN, "--sss", "--in-files", f"{out_dir(arm)}/{locus}.master",
                "--n-causal-snps", str(L_of(arm) if arm in SWEEP_L else P.KMAX_FINEMAP),
                "--n-conv-sss", str(P.N_CONV_SSS_FINEMAP_PRIMARY),
                "--prob-cred-set", cov() if arm in COVERAGE_ARMS else P.COVERAGE, "--log"]
    return cmd


def _rsparsepro_prep(locus):
    """rsparsepro_ld.py reads --z as TAB-separated with columns RSID, Z. The FINEMAP .z is
    space-separated (rsid ... beta se, se = 1 so beta is Z) -- convert, as B8/B16 did."""
    dst = f"{prep_dir('RSparsePro')}/{locus}_rsparsepro_input.txt"
    if os.path.exists(dst):
        return
    with open(f"{MERGED}/{locus}.z") as fin, open(dst + ".tmp", "w") as fout:
        r = csv.DictReader(fin, delimiter=" ", skipinitialspace=True)
        w = csv.writer(fout, delimiter="\t")
        w.writerow(["RSID", "Z"])
        for row in r:
            w.writerow([row["rsid"], row["beta"]])
    os.replace(dst + ".tmp", dst)


def _rsparsepro(locus, info):
    return [PYTHON, RSPARSEPRO_PY,
            "--z", f"{prep_dir('RSparsePro')}/{locus}_rsparsepro_input.txt",
            "--ld", f"{MERGED}/{locus}.ld",
            "--K", str(L_of("RSparsePro")),
            "--cthres", str(P.CTHRES_RSPARSEPRO_PRIMARY),
            "--eincre", str(P.EINCRE_RSPARSEPRO),
            "--save", f"{out_dir('RSparsePro')}/{locus}"]


def _carma(locus, info):
    """run_carma_pilot_b29.R: locus, in_dir, out_dir, rho.index, num.causal, BF.index, seed.
    Reads only {in_dir}/{l}.z and .ld. The seed (B55) makes CARMA's stochastic search reproducible."""
    return [RSCRIPT, os.path.join(SCRIPT_DIR, "run_carma_pilot_b29.R"),
            locus, MERGED, out_dir("CARMA"),
            P.COVERAGE, str(L_of("CARMA")), str(P.BF_INDEX_CARMA), str(P.SEED_CARMA)]


def _eur_susie_prep(locus):
    # The builder also writes {l}.master with n_samples = N_GWAS_PUBLISHED_EUR, which is
    # where run_susie_gctb_aou.R reads N from.
    _build_eur_zld("EUR-SuSiE-RSS", locus)
    if DROP_AMBIGUOUS:
        d = prep_dir("EUR-SuSiE-RSS")
        with open(f"{d}/{locus}.z") as f:
            nz = sum(1 for _ in f) - 1
        with gzip.open(f"{d}/{locus}_snps.txt.gz", "rt") as f:
            ns = sum(1 for _ in f) - 1
        assert nz == ns, f"{locus}: filtered .z has {nz} rows, snps table {ns}"
    os.makedirs(f"{out_dir('EUR-SuSiE-RSS')}/{locus}", exist_ok=True)


def _eur_snps_file(locus):
    if EUR_V2:
        return f"{EUR_V2}/{locus}_snps.txt.gz"
    return (f"{prep_dir('EUR-SuSiE-RSS')}/{locus}_snps.txt.gz" if DROP_AMBIGUOUS
            else f"{EUR_LD}/{locus}_snps.txt.gz")


def _eur_susie(locus, info):
    """susie/run_susie_gctb_aou.R: analysis, locus, finemap_in_dir, snps_file, out_dir,
    L, coverage, min_abs_corr. Its output filenames are NOT locus-specific, so each locus
    gets its own out dir (as batch_native_singleanc_250_qcexcl.py did), and the snps file is
    the per-locus EUR LD snps table."""
    return [RSCRIPT, "--vanilla", os.path.join(FM_DIR, "susie", "run_susie_gctb_aou.R"),
            "eur", locus, prep_dir("EUR-SuSiE-RSS"), _eur_snps_file(locus),
            f"{out_dir('EUR-SuSiE-RSS')}/{locus}",
            str(P.L_SUSIE_RSS), P.COVERAGE, P.MIN_PURITY]


def _noprep(locus):
    return None


ARMS = {
    # name            group           prep                        cmd              done-file
    "SuSiEx":         ("per-ancestry", _noprep,        _susiex,         "{o}/{l}_susiex.cs"),
    "MultiSuSiE":     ("per-ancestry", _noprep,        _multisusie,     "{o}/{l}_multisusie_pip.tsv"),
    "MESuSiE":        ("per-ancestry", _noprep,        _mesusie,        "{o}/{l}_mesusie_pip.tsv"),
    "SuSiE-merged":   ("merged",  _merged_prep("SuSiE-merged"), _susie_merged, "{o}/{l}_susie_pip.tsv"),
    "FINEMAP-merged": ("merged",  _finemap_prep("FINEMAP-merged", MERGED, P.N_GWAS_MERGED),
                                                       _finemap("FINEMAP-merged"), "{o}/{l}.snp"),
    "RSparsePro":     ("merged",       _rsparsepro_prep, _rsparsepro,   "{o}/{l}.rsparsepro.txt"),
    "CARMA":          ("merged",       _noprep,        _carma,          "{o}/{l}_carma_pip.tsv"),
    "EUR-SuSiE-RSS":  ("EUR-only",     _eur_susie_prep, _eur_susie,
                       "{o}/{l}/eur_susie_gctb_aou_pip_all.tsv"),
    # The original never passed --prob-cred-set here (FINEMAP's 0.95 default,
    # ORIGINAL_RUN_SPEC 3.4). Passing COVERAGE aligns the arm; deliberate departure.
    "EUR-FINEMAP":    ("EUR-only", _eur_finemap_prep,
                                                       _finemap("EUR-FINEMAP"), "{o}/{l}.snp"),
}


def assert_no_literals():
    """Every aligned parameter must come from run_params, not from this file."""
    required = ["COVERAGE", "MIN_PURITY", "MIN_PURITY_ZERO", "L_MULTISUSIE", "L_MESUSIE",
                "L_SUSIE_RSS", "N_SIG_SUSIEX", "KMAX_FINEMAP", "RHO_MULTISUSIE",
                "PVAL_THRESH_SUSIEX", "MULT_STEP_SUSIEX_PUBLISHED", "K_RSPARSEPRO",
                "CTHRES_RSPARSEPRO_PRIMARY", "EINCRE_RSPARSEPRO", "N_GWAS_PUBLISHED_EUR",
                "N_GWAS_PUBLISHED_EAS", "N_GWAS_MERGED", "N_CONV_SSS_FINEMAP_PRIMARY",
                "NUM_CAUSAL_CARMA", "BF_INDEX_CARMA", "SEED_CARMA", "THREAD_LIMIT_ENV"]
    missing = [k for k in required if not hasattr(P, k)]
    if missing:
        sys.exit(f"run_params.py is missing required constants: {missing}")
    if P.COVERAGE != "0.99":
        sys.exit(f"COVERAGE is {P.COVERAGE!r}; the 2026-09 run is signed off at 0.99.")
    if str(P.CTHRES_RSPARSEPRO_PRIMARY) != P.COVERAGE:
        sys.exit("RSparsePro --cthres is its coverage and must equal COVERAGE.")
    if L_OVERRIDE is not None:
        for k in ("L_SWEEP_VALUES", "MULT_STEP_SUSIEX_SWEEP"):
            if not hasattr(P, k):
                sys.exit(f"run_params.py is missing required constant {k}")


def _argv_L(arm, argv):
    """Read the L-equivalent back out of a built argv, by the runner's own contract."""
    def after(flag):
        return argv[argv.index(flag) + 1] if flag in argv else None
    if arm == "SuSiEx":
        v = [a.split("=", 1)[1] for a in argv if a.startswith("--n_sig=")]
        return v[0] if len(v) == 1 else None
    return {"MultiSuSiE": lambda: after("--L"), "MESuSiE": lambda: argv[-1],
            "SuSiE-merged": lambda: argv[-1], "FINEMAP-merged": lambda: after("--n-causal-snps"),
            "RSparsePro": lambda: after("--K"), "CARMA": lambda: argv[-3]}[arm]()


def _argv_cov(arm, argv):
    """Read the coverage back out of a built argv (SuSiE-merged: 3 before the end; FINEMAP)."""
    if arm == "SuSiE-merged":
        return argv[-3]
    if arm == "FINEMAP-merged":
        return argv[argv.index("--prob-cred-set") + 1]
    return None


def assert_override_reached(arms, locus):
    """With --L: every selected arm's argv must carry L, write under WORK/<arm>_L{n}, and
    (SuSiEx) run mult-step off. A sweep that silently runs at one L produces a flat curve."""
    if L_OVERRIDE is None and COV_OVERRIDE is None:
        return
    bad = []
    for a in arms:
        argv = ARMS[a][2](locus, LOCI_INFO[locus])
        if L_OVERRIDE is not None:
            got = _argv_L(a, argv)
            if got != str(L_OVERRIDE):
                bad.append(f"{a}: argv carries L={got!r}, not {L_OVERRIDE}")
        if COV_OVERRIDE is not None:
            got = _argv_cov(a, argv)
            if got != COV_OVERRIDE:
                bad.append(f"{a}: argv carries coverage={got!r}, not {COV_OVERRIDE}")
        tag = (a if L_OVERRIDE is None else f"{a}_L{L_OVERRIDE}") + cov_tag()
        if not any(tag in x for x in argv) or out_dir(a) != os.path.join(WORK, tag):
            bad.append(f"{a}: output not under {os.path.join(WORK, tag)}")
        # never into a primary or any L-sweep dir (cov dirs carry the _cov suffix)
        toks = [t for x in argv for y in x.split("=") for t in y.split(",")]
        prot = [os.path.join(WORK, a)] + [os.path.join(WORK, f"{a}_L{n}") for n in P.L_SWEEP_VALUES]
        for prim in prot:
            if any(t == prim or t.startswith(prim + os.sep) for t in toks):
                bad.append(f"{a}: argv points into the protected output dir {prim}")
        if COV_OVERRIDE is not None and not out_dir(a).endswith(cov_tag()):
            bad.append(f"{a}: output dir {out_dir(a)} lacks the {cov_tag()} suffix")
        if a == "SuSiEx" and L_OVERRIDE is not None:
            # checked against the behaviour, not the constant: mult-step ON ignores --n_sig,
            # so a wrong MULT_STEP_SUSIEX_SWEEP must not be able to vouch for itself
            ms = [x.split("=", 1)[1] for x in argv if x.startswith("--mult-step=")]
            if len(ms) != 1 or ms[0].lower() != "false":
                bad.append(f"SuSiEx: --mult-step={ms} is not off; --n_sig would be ignored")
    if bad:
        sys.exit("L override did not reach the argv:\n  " + "\n  ".join(bad))


def assert_tools(arms):
    need = {"SuSiEx": [SUSIEX_BIN, PLINK_BIN], "RSparsePro": [RSPARSEPRO_PY, PYTHON],
            "FINEMAP-merged": [FINEMAP_BIN], "EUR-FINEMAP": [FINEMAP_BIN],
            "MultiSuSiE": [PYTHON]}
    miss = sorted({p for a in arms for p in need.get(a, [RSCRIPT]) if not os.path.exists(p)})
    if miss:
        sys.exit(f"missing tools (run envs/setup_ilae4_env.sh): {miss}")


def stage_inputs(loci, all_loci):
    """Copy the inputs these loci need from the buckets into INPUTS. Whole prefixes when
    running all 250, per-locus files otherwise. Copy only; nothing at source is touched."""
    if V2:
        sys.exit("v2 inputs are built on the cluster by harmonize_v2.py, not staged from the bucket")
    for d in (UNION, INTERSECT, MERGED, EUR_LD, EUR_SS):
        os.makedirs(d, exist_ok=True)
    if all_loci:
        jobs = [(f"{BUCKET_LD}/ld_eur_eas_union/", UNION, None),
                (f"{BUCKET_LD}/ld_eur_eas_intersect/", INTERSECT, None),
                (f"{BUCKET_LD}/ld_eur_eas_merged/", MERGED, None),
                (f"{BUCKET_LD}/ld_per_anc_eur/", EUR_LD, None),
                (f"{BUCKET_EUR_SS}/", EUR_SS, None)]
        for src, dst, excl in jobs:
            cmd = ["gsutil", "-m", "-q", "rsync", "-r"] + (["-x", excl] if excl else []) + [src, dst]
            subprocess.run(cmd, check=True)
        return
    for l in loci:
        srcs = [(f"{BUCKET_LD}/ld_eur_eas_union/{l}_harmonized.npz", UNION),
                (f"{BUCKET_LD}/ld_eur_eas_intersect/{l}.*", INTERSECT),
                (f"{BUCKET_LD}/ld_eur_eas_merged/{l}.z", MERGED),
                (f"{BUCKET_LD}/ld_eur_eas_merged/{l}.ld", MERGED),
                (f"{BUCKET_LD}/ld_per_anc_eur/{l}_snps.txt.gz", EUR_LD),
                (f"{BUCKET_LD}/ld_per_anc_eur/{l}_R.npz", EUR_LD),
                (f"{BUCKET_EUR_SS}/{l}_sumstats.tsv.gz", EUR_SS)]
        for src, dst in srcs:
            subprocess.run(["gsutil", "-m", "-q", "cp", src, dst + "/"], check=True)
        subprocess.run(["gsutil", "-m", "-q", "cp", "-r",
                        f"{BUCKET_LD}/ld_eur_eas_union/{l}_susiex", UNION + "/"], check=True)


LOCI_INFO = {}


def run_one(arm, locus):
    o = out_dir(arm)
    os.makedirs(o, exist_ok=True)
    os.makedirs(prep_dir(arm), exist_ok=True)
    os.makedirs(f"{o}/_logs", exist_ok=True)
    _, prep, builder, out_tpl = ARMS[arm]
    target = out_tpl.format(o=o, l=locus)
    if os.path.exists(target):
        return arm, locus, "skip", ""
    try:
        prep(locus)
    except Exception as e:
        return arm, locus, "FAIL", f"prep: {e!r}"[-400:]
    cmd = builder(locus, LOCI_INFO[locus])
    env = {**os.environ, **P.THREAD_LIMIT_ENV}
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    dt = time.time() - t0
    with open(f"{o}/_logs/{locus}.log", "w") as f:
        f.write("argv: " + " ".join(cmd) + f"\nexit: {p.returncode}  wall_s: {dt:.1f}\n\n")
        f.write(p.stdout + "\n--- stderr ---\n" + p.stderr)
    if p.returncode != 0:
        return arm, locus, "FAIL", (p.stderr or p.stdout)[-400:]
    if not os.path.exists(target):
        return arm, locus, "FAIL", f"exit 0 but no output {target}"
    return arm, locus, "ok", f"{dt:.1f}s"


# Outputs go under the manuscript bucket's own layout, next to the earlier runs in
# fine_mapping_runs/. `_inputs/` is excluded: it holds symlinks to ld_matrices/ (rsync would
# follow them and duplicate the LD) and derived dumps that the prep step rebuilds.
BACKUP_DEFAULT = (f"{WORKSPACE_BUCKET}/manuscript_final_2026-09/"
                  "fine_mapping_runs/primary_2026-09")
BACKUP_EXCLUDE = r"(^|.*/)_inputs/.*"


BACKUP_V2_DEFAULT = (f"{WORKSPACE_BUCKET}/manuscript_final_2026-09/"
                     "fine_mapping_runs/primary_v2_2026-09")

BACKUP_SWEEP_DEFAULT = (f"{WORKSPACE_BUCKET}/manuscript_final_2026-09/"
                        "fine_mapping_runs/l_sweep_2026-09")


BACKUP_COV_DEFAULT = (f"{WORKSPACE_BUCKET}/manuscript_final_2026-09/"
                      "fine_mapping_runs/cov095_2026-10")


def backup_arm(arm, prefix):
    dst = f"{prefix}/{os.path.basename(out_dir(arm))}/"
    r = subprocess.run(["gsutil", "-m", "-q", "rsync", "-r", "-c", "-x", BACKUP_EXCLUDE,
                        out_dir(arm) + "/", dst], capture_output=True, text=True)
    print(f"  [backup] {arm} -> {dst} exit={r.returncode} {r.stderr[-200:]}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arms", default=",".join(ARMS), help="comma-separated subset")
    ap.add_argument("--loci", default=None, help="comma-separated subset, else all 250")
    ap.add_argument("--loci-file", default=os.path.expanduser(
        "~/work/ref_inputs/loci_250_pgc3scz.tsv"))
    ap.add_argument("--jobs", type=int,
                    default=int(os.environ.get("N_CONCURRENT", max(1, os.cpu_count() - 2))))
    ap.add_argument("--dry-run", action="store_true",
                    help="print the exact argv per arm for the first locus, run nothing")
    ap.add_argument("--stage-inputs", action="store_true",
                    help="copy the needed inputs from gs:// into PRIMARY_INPUTS, then exit")
    ap.add_argument("--backup-prefix", default=None,
                    help="gs:// prefix; each arm is rsynced there when its last locus ends "
                         "(default: primary_2026-09/, or l_sweep_2026-09/ with --L)")
    ap.add_argument("--L", type=int, default=None,
                    help="L sweep: override the L-equivalent of the selected SWEEP_L arms "
                         "and write to WORK/<arm>_L{n}")
    ap.add_argument("--coverage", default=None,
                    help="B69: override the credible-set coverage of SuSiE-merged / FINEMAP-merged "
                         "(e.g. 0.95); writes WORK/<arm>[_L{n}]_cov095, run_params.COVERAGE untouched")
    ap.add_argument("--no-backup", action="store_true", help="skip the gs:// backup (tests)")
    ap.add_argument("--prep-only", action="store_true",
                    help="run each arm's prep step for every locus, fit nothing (B34 §3 "
                         "checks the built inputs before any fit)")
    args = ap.parse_args()

    global L_OVERRIDE, COV_OVERRIDE
    L_OVERRIDE = args.L
    COV_OVERRIDE = args.coverage
    if COV_OVERRIDE is not None:
        if COV_OVERRIDE == P.COVERAGE:
            sys.exit(f"--coverage {COV_OVERRIDE} equals run_params.COVERAGE: that is the primary run")
        off = [a for a in args.arms.split(",") if a.strip() not in COVERAGE_ARMS]
        if off:
            sys.exit(f"--coverage applies only to {list(COVERAGE_ARMS)}; got {off} (pass --arms)")
    state_path = os.path.join(WORK, ("state.json" if args.L is None else f"state_L{args.L}.json")
                              .replace(".json", cov_tag() + ".json"))
    if COV_OVERRIDE is not None:
        if args.backup_prefix is None:
            args.backup_prefix = f"{BACKUP_COV_DEFAULT}/{V2_VARIANT or 'v1'}"
        elif "primary_" in args.backup_prefix or "l_sweep" in args.backup_prefix:
            sys.exit(f"--coverage must not back up under a primary/L-sweep prefix: {args.backup_prefix}")
    elif args.backup_prefix is None:
        args.backup_prefix = (f"{BACKUP_V2_DEFAULT}/{V2_VARIANT}" if V2 and args.L is None else
                              BACKUP_DEFAULT if args.L is None else BACKUP_SWEEP_DEFAULT)
    if V2 and DROP_AMBIGUOUS:
        sys.exit("v2 selects no-palindromic by PRIMARY_V2_INPUTS=.../noambig; unset PRIMARY_DROP_AMBIGUOUS")
    assert_no_literals()
    if DROP_AMBIGUOUS and os.path.normpath(WORK) == os.path.normpath(
            os.path.expanduser("~/work/primary_2026")):
        sys.exit("PRIMARY_DROP_AMBIGUOUS=1 needs its own PRIMARY_WORK, not ~/work/primary_2026")
    arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    bad = [a for a in arms if a not in ARMS]
    if bad:
        sys.exit(f"unknown arm(s): {bad}; known: {list(ARMS)}")
    if args.L is not None:
        if args.L not in P.L_SWEEP_VALUES:
            sys.exit(f"--L {args.L} is not in run_params.L_SWEEP_VALUES {P.L_SWEEP_VALUES}")
        if args.arms == ",".join(ARMS):
            arms = list(SWEEP_L)
        off = [a for a in arms if a not in SWEEP_L]
        if off:
            sys.exit(f"--L applies only to the sweep arms {list(SWEEP_L)}; got {off}")
        reuse = [a for a in arms if primary_L(a) == args.L]
        if reuse:
            sys.exit(f"the primary run already is the L={args.L} point for {reuse}: reuse "
                     f"{[os.path.join(WORK, a) for a in reuse]}, do not recompute")

    # SuSiEx needs --chr/--bp, so every arm gets the locus's coordinates
    import pandas as pd
    t = pd.read_csv(args.loci_file, sep="\t")
    for r in t.itertuples():
        LOCI_INFO[f"chr{r.chr}_{r.start}_{r.stop}"] = dict(chrom=r.chr, start=r.start,
                                                           stop=r.stop)
    if args.loci:
        loci = [l.strip() for l in args.loci.split(",") if l.strip()]
        miss = [l for l in loci if l not in LOCI_INFO]
        if miss:
            sys.exit(f"locus not in {args.loci_file}: {miss}")
    else:
        loci = list(LOCI_INFO)
    assert_override_reached(arms, loci[0])

    if args.stage_inputs:
        stage_inputs(loci, all_loci=args.loci is None)
        print(f"staged inputs for {len(loci)} loci into {INPUTS}")
        return

    if args.dry_run:
        print(f"WORK={WORK}\nV2={V2}\nMULTISUSIE_NPZ={MULTISUSIE_NPZ}\nINPUTS={INPUTS}\nENV={ENV}\nL override={L_OVERRIDE}\n"
              f"coverage={cov()}  purity={P.MIN_PURITY}  "
              f"(SuSiEx runs at {P.MIN_PURITY_ZERO}, gated post hoc)\n")
        for a in arms:
            grp, _, builder, out_tpl = ARMS[a]
            print(f"--- {a}  [{grp}]")
            print("    " + " ".join(builder(loci[0], LOCI_INFO[loci[0]])))
            print(f"    -> {out_tpl.format(o=out_dir(a), l=loci[0])}\n")
        return

    if args.prep_only:
        def prep_one(a, l):
            os.makedirs(prep_dir(a), exist_ok=True)
            try:
                ARMS[a][1](l)
                return a, l, ""
            except Exception as e:  # noqa: BLE001
                return a, l, repr(e)[-300:]
        with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
            res = list(ex.map(lambda al: prep_one(*al), [(a, l) for a in arms for l in loci]))
        bad = [r for r in res if r[2]]
        for a, l, e in bad:
            print(f"  PREP FAIL {a} {l}: {e}", flush=True)
        print(f"prep-only: {len(res) - len(bad)} ok, {len(bad)} failed")
        return

    assert_tools(arms)
    os.makedirs(WORK, exist_ok=True)
    # Largest loci first within each arm, so the long jobs do not trail at the end.
    size = {l: os.path.getsize(f"{MERGED}/{l}.ld") if os.path.exists(f"{MERGED}/{l}.ld") else 0
            for l in loci}
    tasks = [(a, l) for a in arms for l in sorted(loci, key=lambda x: -size[x])]
    remaining = {a: len(loci) for a in arms}
    lock = threading.Lock()
    print(f"{len(arms)} arms x {len(loci)} loci = {len(tasks)} jobs, --jobs {args.jobs}",
          flush=True)
    results, t0 = [], time.time()
    backups = []
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(run_one, a, l) for a, l in tasks]
        for fu in cf.as_completed(futs):
            arm, locus, status, msg = fu.result()
            results.append(dict(arm=arm, locus=locus, status=status, msg=msg))
            if status == "FAIL":
                print(f"  FAIL {arm} {locus}: {msg}", flush=True)
            with lock:
                remaining[arm] -= 1
                done = remaining[arm] == 0
            if done:
                print(f"  [arm complete] {arm}", flush=True)
                if args.backup_prefix and not args.no_backup:
                    th = threading.Thread(target=backup_arm, args=(arm, args.backup_prefix))
                    th.start()
                    backups.append(th)
            with open(state_path, "w") as f:
                json.dump(results, f, indent=1)
    for th in backups:
        th.join()

    print(f"\ndone in {time.time()-t0:.0f}s -> {state_path}")
    for a in arms:
        r = [x for x in results if x["arm"] == a]
        n_ok = sum(x["status"] == "ok" for x in r)
        n_sk = sum(x["status"] == "skip" for x in r)
        n_f = sum(x["status"] == "FAIL" for x in r)
        print(f"  {a:<16} ok={n_ok:<4} skipped={n_sk:<4} FAILED={n_f}")


if __name__ == "__main__":
    main()
