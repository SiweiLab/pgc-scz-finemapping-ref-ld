#!/usr/bin/env python3
"""v2 step 3 for all loci in one pass: ancestry LD from the AoU ACAF splitMT.

Stage A (Hail, one session): read the splitMT once, filter to one point interval per target
  position, semi-join rows to the targets' (locus, alleles) in both orderings, semi-join columns
  to the ancestry sample set (asserted), keep GT, EAS-subset variant_qc, flag fail_cr / fail_mono,
  checkpoint (before the row filter, so the per-locus filter counts come from the checkpoint
  without a second scan).
Stage B: per chromosome (or per locus when a chromosome's int8 matrix would exceed --collect-gb),
  collect ALT dosage as int8 (missing = -1) with its row table.
Stage C: per locus, replicate hl.row_correlation in numpy: mean-impute missing per variant,
  center, scale to unit L2 norm, R = X X^T in float64, stored float32.
Stage D: write the extractor-format files ({locus}_eas_R.npz + {locus}_eas_snps.tsv + {locus}.z)
  and run align_ld_aou.py (3c) into OUT/ld_v2/eas_aou/.

Acceptance before anything is uploaded: the --pilot loci must equal a reference tree (an earlier
staged R) to max |delta| <= 1e-6, with identical variant lists, order and sign.

Usage:
  python build_ld_aou.py --out ~/work/ld_aou --ckpt gs://.../eas_targets.mt
      [--ancestry-code eas] [--expected-n 12374] [--exclude-samples gs://...] [--pilot L1,..]
      [--pilot-ref DIR] [--point-intervals | --locus-intervals] [--collect-gb 2]
No sample ID is written to OUT (the checkpoint's column keys stay in the bucket).
"""
import argparse
import glob
import json
import os
WORKSPACE_BUCKET = os.environ.get("WORKSPACE_BUCKET", "gs://YOUR-WORKSPACE-BUCKET")  # AoU Workbench sets this
import subprocess
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ACAF = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/splitMT/hail.mt"
ANC = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
EXCL = (f"{WORKSPACE_BUCKET}/manuscript_final_2026-09/aou_samples/"
        "combined_exclude_samples.txt")
TG_GS = f"{WORKSPACE_BUCKET}/manuscript_final_2026-09/ld_v2_extract/all_targets.tsv"
PILOT = ("chr12_110129434_111489434,chr2_97500526_98599526,chr3_48182866_53553354,"
         "chr8_144792601_144940101,chr22_40612819_42707289")

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--ckpt", required=True)
ap.add_argument("--ancestry-code", default="eas")
ap.add_argument("--expected-n", type=int, default=12374)
ap.add_argument("--exclude-samples", default=EXCL)
ap.add_argument("--min-call-rate", type=float, default=0.95)
ap.add_argument("--pilot", default=PILOT)
ap.add_argument("--pilot-ref", default=os.path.expanduser("~/work/ld_aou_pilot_ref/ld_v2/eas_aou"))
ap.add_argument("--locus-intervals", action="store_true", help="fallback: 250 locus windows, not point intervals")
ap.add_argument("--collect-gb", type=float, default=2.0)
ap.add_argument("--n-partitions", type=int, default=1000,
                help="naive_coalesce before the Stage A checkpoint (no shuffle; no-op if already fewer)")
ap.add_argument("--pilot-check", choices=["reference", "eur_gctb"], default="reference",
                help="reference: EAS against --pilot-ref; eur_gctb: AoU EUR against the GCTB EUR LD")
ap.add_argument("--gctb-dir", default=os.path.expanduser("~/work/ld_eur_gctb/ld_v2/eur_gctb"))
ap.add_argument("--bm-tmp", default="/tmp/b48_bm", help="default-FS (HDFS) dir for per-unit BlockMatrix writes")
ap.add_argument("--pilot-only", action="store_true", help="stop after the pilot check")
ap.add_argument("--resume", action="store_true", help="skip loci whose extract npz already exists")
ap.add_argument("--skip-a", action="store_true", help="reuse an existing checkpoint")
ap.add_argument("--fixed-ckpt", default=None, help="unfiltered, recomputed-QC checkpoint (default: --ckpt with _fixed)")
a = ap.parse_args()
OUT = os.path.expanduser(a.out)
a.fixed_ckpt = a.fixed_ckpt or a.ckpt.replace(".mt", "_fixed.mt")
AC = a.ancestry_code.lower()      # file names, the ld_v2/{AC}_aou tree and the {AC}_af column
for d in ("extract", "logs"):
    os.makedirs(os.path.join(OUT, d), exist_ok=True)
LOG = open(os.path.join(OUT, "logs", "onepass.log"), "a")
TIMING = {}


def say(*x):
    s = time.strftime("%H:%M:%S ") + " ".join(str(v) for v in x)
    print(s, flush=True)
    LOG.write(s + "\n")
    LOG.flush()


def timed(name, t0):
    TIMING[name] = round(time.time() - t0, 1)
    say(f"[time] {name}: {TIMING[name]} s")
    json.dump(TIMING, open(os.path.join(OUT, "logs", "timing.json"), "w"), indent=1)


def cluster_shape():
    r = subprocess.run("yarn node -list 2>/dev/null | head -1", shell=True, capture_output=True, text=True)
    return r.stdout.strip()


import hail as hl  # noqa: E402
from hail_ld import load_eas_ids, load_exclude_ids  # noqa: E402

hl.init(default_reference="GRCh38", tmp_dir=f"/tmp/b45_hail_{a.ancestry_code}", log=os.path.join(OUT, "logs", "hail.log"),
        gcs_requester_pays_configuration=os.environ["GOOGLE_PROJECT"])

# ---- targets -------------------------------------------------------------------------------
tg = pd.read_csv(os.path.join(HERE, "..", "inputs", "hg38_targets", "intervals.tsv"), sep="\t")
loci = tg["locus"].tolist()
ht = hl.import_table(TG_GS, types={"position": hl.tint32})
ht = ht.rename({"locus": "locus_id"})
ht = ht.annotate(locus=hl.locus("chr" + ht.chromosome, ht.position),
                 allele1=ht.allele1.upper(), allele2=ht.allele2.upper())
exact = (ht.annotate(alleles=[ht.allele1, ht.allele2])
         .union(ht.annotate(alleles=[ht.allele2, ht.allele1]))
         .key_by("locus", "alleles").select("rsid", "locus_id"))
tloc = pd.concat(pd.read_csv(f, sep="\t", dtype=str).assign(locus_id=os.path.basename(f).replace("_hg38_targets.tsv", ""))
                 for f in sorted(glob.glob(os.path.join(HERE, "..", "inputs", "hg38_targets", "*_hg38_targets.tsv"))))
n_targets = tloc.groupby("locus_id").size()

# ---- Stage A --------------------------------------------------------------------------------
if not a.skip_a:
    say(f"Stage A start; cluster: {cluster_shape()}")
    ids = load_eas_ids(ANC, a.ancestry_code) - load_exclude_ids(a.exclude_samples)
    say(f"[samples] {a.ancestry_code} after exclusion: {len(ids)}")
    assert len(ids) == a.expected_n, f"{len(ids)} samples, expected {a.expected_n}"
    samp = hl.Table.parallelize([hl.struct(s=x) for x in ids], hl.tstruct(s=hl.tstr)).key_by("s")

    t0 = time.time()
    mt = hl.read_matrix_table(ACAF).drop("variant_qc", "info")
    if a.locus_intervals:
        ivs = [hl.parse_locus_interval(x) for x in tg["interval"]]
        say(f"[A] {len(ivs)} locus intervals")
    else:
        pos = tloc[["chromosome", "position"]].drop_duplicates()
        ivs = [hl.Interval(hl.Locus("chr" + c, int(p), "GRCh38"), hl.Locus("chr" + c, int(p), "GRCh38"),
                           includes_start=True, includes_end=True) for c, p in zip(pos["chromosome"], pos["position"])]
        say(f"[A] {len(ivs):,} point intervals")
    mt = hl.filter_intervals(mt, ivs)
    mt = mt.semi_join_rows(exact)
    mt = mt.semi_join_cols(samp)
    mt = mt.select_entries("GT")
    mt = mt.annotate_rows(vkey=exact[mt.row_key].rsid, locus_id=exact[mt.row_key].locus_id)
    # The ACAF MT carries filtered entries (FT). variant_qc already counts them as not
    # called, but no other GT expression sees them (an earlier Stage B crash), so count them, then
    # unfilter them into missing GT before variant_qc; the checkpoint then needs no fix pass.
    mt = mt.annotate_rows(n_filtered_entries=hl.int64(a.expected_n) - hl.agg.count())
    mt = mt.unfilter_entries()
    mt = hl.variant_qc(mt)
    mt = mt.annotate_rows(fail_cr=mt.variant_qc.call_rate < a.min_call_rate,
                          fail_mono=(mt.variant_qc.AF[1] == 0) | (mt.variant_qc.AF[1] == 1))
    n_part = mt.n_partitions()
    mt = mt.naive_coalesce(a.n_partitions)
    say(f"[A] partitions before naive_coalesce {n_part:,}, after {mt.n_partitions():,}")
    timed("A_plan_python", t0)
    t0 = time.time()
    mt = mt.checkpoint(a.ckpt, overwrite=True)
    timed("A_scan_and_checkpoint", t0)
    say(f"Stage A done; cluster: {cluster_shape()}")

mt = hl.read_matrix_table(a.ckpt)
n_cols = mt.count_cols()
say(f"[A] checkpoint: {mt.count_rows():,} rows x {n_cols:,} samples")
assert n_cols == a.expected_n

# The ACAF MT carries filtered entries (FT). Hail's aggregators skip them, so
# Stage A's variant_qc call rates were too high, and or_else() never runs on them (the Stage B
# crash). Count them per row, unfilter them into missing GT, and recompute variant_qc and both
# row filters on the unfiltered MT. AF is unchanged (filtered entries were never counted as
# called). The lenient checkpoint holds every variant the corrected filter keeps.
# (variant_qc's call rate was already right -- n_filtered is in its denominator -- so the
# recompute changes nothing; only the unfilter matters.) A checkpoint written after this fix is already unfiltered.
LEGACY = "n_filtered_entries" not in mt.row.dtype
if LEGACY:
    t0 = time.time()
    mt = mt.annotate_rows(n_filtered_entries=n_cols - hl.agg.count(), fail_cr_stage_a=mt.fail_cr)
    mt = mt.unfilter_entries().drop("variant_qc", "fail_cr", "fail_mono")
    mt = hl.variant_qc(mt)
    mt = mt.annotate_rows(fail_cr=mt.variant_qc.call_rate < a.min_call_rate,
                          fail_mono=(mt.variant_qc.AF[1] == 0) | (mt.variant_qc.AF[1] == 1))
    mt = mt.checkpoint(a.fixed_ckpt, overwrite=True)
    timed("fix_unfilter_recompute_checkpoint", t0)
else:
    mt = mt.annotate_rows(fail_cr_stage_a=mt.fail_cr)
n_filt = mt.aggregate_rows(hl.agg.sum(mt.n_filtered_entries))
say(f"[A] filtered entries (now missing GT) {n_filt:,} of {mt.count_rows() * n_cols:,}")

# Filter counts per locus, from the checkpoint.
t0 = time.time()
r = mt.rows()
kept = ~r.fail_cr & ~r.fail_mono
fc = r.group_by(r.locus_id).aggregate(
    matched_in_acaf=hl.agg.count(),
    n_filtered_entries=hl.agg.sum(r.n_filtered_entries),
    fail_callrate_recomputed=hl.agg.count_where(r.fail_cr & ~r.fail_cr_stage_a),
    kept_stage_a=hl.agg.count_where(~r.fail_cr_stage_a & ~r.fail_mono),
    fail_callrate_only=hl.agg.count_where(r.fail_cr & ~r.fail_mono),
    fail_monomorphic_only=hl.agg.count_where(~r.fail_cr & r.fail_mono),
    fail_both=hl.agg.count_where(r.fail_cr & r.fail_mono),
    kept=hl.agg.count_where(kept),
    callrate_min_kept=hl.agg.filter(kept, hl.agg.min(r.variant_qc.call_rate)),
    af_min_kept=hl.agg.filter(kept, hl.agg.min(r.variant_qc.AF[1])),
    af_max_kept=hl.agg.filter(kept, hl.agg.max(r.variant_qc.AF[1]))).to_pandas().set_index("locus_id")
fc = fc.rename(columns={"af_min_kept": f"{AC}_af_min_kept", "af_max_kept": f"{AC}_af_max_kept"})
fc.insert(0, "targets", n_targets.reindex(fc.index))
fc = fc.reindex(loci)
fc["targets"] = n_targets.reindex(loci).values
fc = fc.fillna({"matched_in_acaf": 0, "n_filtered_entries": 0, "fail_callrate_recomputed": 0, "kept_stage_a": 0,
               "fail_callrate_only": 0, "fail_monomorphic_only": 0, "fail_both": 0, "kept": 0})
fc.index.name = "locus"
fc.to_csv(os.path.join(OUT, "filter_counts_per_locus.csv"))
tot = fc[["targets", "matched_in_acaf", "n_filtered_entries", "fail_callrate_recomputed", "kept_stage_a",
          "fail_callrate_only", "fail_monomorphic_only", "fail_both", "kept"]].sum()
say("[counts] totals: " + ", ".join(f"{k} {int(v):,}" for k, v in tot.items()))
timed("A_filter_counts", t0)

mt = mt.filter_rows(~mt.fail_cr & ~mt.fail_mono)
mt = mt.annotate_rows(af=mt.variant_qc.AF[1], call_rate=mt.variant_qc.call_rate)


# ---- Stages B + C + D per collect unit ------------------------------------------------------
def corr(G, chunk_bytes=2e9):
    """hl.row_correlation in numpy: per variant, mean-impute missing (-1) with the mean over
    called samples, center, scale to unit L2 norm; R = X X^T in float64, returned as float32.
    Summed over sample chunks, so peak memory is ~k x chunk x 8 bytes, not k x n (EUR's
    chr3 locus is 10k x 267k). Same arithmetic as the earlier one-shot version up to summation order."""
    k, n = G.shape
    step = max(256, int(chunk_bytes / 8 / max(k, 1)))
    s = np.zeros(k)
    c = np.zeros(k)
    for j in range(0, n, step):
        g = G[:, j:j + step]
        ok = g >= 0
        s += np.where(ok, g, 0).sum(1, dtype=np.float64)
        c += ok.sum(1)
    mean = s / c
    P = np.zeros((k, k))
    ss = np.zeros(k)
    for j in range(0, n, step):
        g = G[:, j:j + step]
        X = np.where(g >= 0, g, 0).astype(np.float64)
        X -= mean[:, None]
        X[g < 0] = 0.0          # imputed to the mean -> 0 after centering
        ss += (X * X).sum(1)
        P += X @ X.T
    d = 1.0 / np.sqrt(ss)
    return (P * d[:, None] * d[None, :]).astype(np.float32)


def write_locus(L, rows, R):
    d = os.path.join(OUT, "extract", L)
    os.makedirs(d, exist_ok=True)
    z = pd.DataFrame({"rsid": rows["vkey"].values, "chromosome": rows["contig"].str.replace("chr", "", regex=False).values,
                      "position": rows["position"].values, "allele1": rows["alt"].values, "allele2": rows["ref"].values,
                      f"{AC}_af": rows["af"].values, "call_rate": rows["call_rate"].values})
    z.to_csv(os.path.join(d, f"{L}_{AC}_snps.tsv"), sep=" ", index=False)
    z.to_csv(os.path.join(d, f"{L}.z"), sep=" ", index=False)
    np.savez_compressed(os.path.join(d, f"{L}_{AC}_R.npz"), R=R)   # last: its presence marks the locus done


TC_PATH = os.path.join(OUT, "logs", "per_locus_time.csv")
FAILED = []


def collect(sub, n):
    """ALT dosage as int8 (missing -> -1), rows in key order. The BlockMatrix is written to an
    explicit default-FS path, read back 1024 rows at a time (each read < 2^31 entries, so
    to_numpy stays on its local-file path), then deleted."""
    path = f"{a.bm_tmp}/{int(time.time() * 1000)}.bm"
    hl.linalg.BlockMatrix.write_from_entry_expr(hl.or_else(sub.GT.n_alt_alleles(), -1), path, block_size=1024)
    bm = hl.linalg.BlockMatrix.read(path)
    G = np.empty((bm.n_rows, bm.n_cols), dtype=np.int8)
    for i in range(0, bm.n_rows, 1024):
        G[i:i + 1024] = bm[i:min(i + 1024, bm.n_rows), :].to_numpy().astype(np.int8)
    hl.current_backend().fs.rmtree(path)
    assert G.shape[0] == n
    return G


def run_units(units, tag):
    tc = {}
    for chrom, ls in units:
        if a.resume:
            ls = [L for L in ls if not os.path.exists(os.path.join(OUT, "extract", L, f"{L}_{AC}_R.npz"))]
            if not ls:
                continue
        try:
            tc.update(run_unit(chrom, ls, tag))
        except Exception as e:  # noqa: BLE001 -- unattended run: log, carry on with the next unit
            FAILED.extend(ls)
            say(f"[FAIL {tag}] {chrom or ls[0]} ({len(ls)} loci): {type(e).__name__}: {str(e)[:300]}")
    return tc


def run_unit(chrom, ls, tag):
    tc = {}
    t0 = time.time()
    if chrom is None:  # per-locus collect
        sub = hl.filter_intervals(mt, [hl.parse_locus_interval(tg.set_index("locus").at[ls[0], "interval"])])
        sub = sub.filter_rows(sub.locus_id == ls[0])
    else:
        sub = hl.filter_intervals(mt, [hl.parse_locus_interval(chrom)])
        sub = sub.filter_rows(hl.literal(set(ls)).contains(sub.locus_id))
    sub = sub.annotate_rows(contig=sub.locus.contig, position=sub.locus.position,
                            ref=sub.alleles[0], alt=sub.alleles[1])
    rows = sub.rows().key_by().select("vkey", "locus_id", "contig", "position", "ref", "alt",
                                      "af", "call_rate").to_pandas()
    G = collect(sub, len(rows))
    assert G.shape == (len(rows), n_cols)
    tb = time.time() - t0
    say(f"[B {tag}] {chrom or ls[0]}: collected {G.shape[0]:,} x {G.shape[1]:,} int8 in {tb:.1f} s")
    for L in ls:
        t1 = time.time()
        m = (rows["locus_id"] == L).to_numpy()
        R = corr(G[m])
        write_locus(L, rows[m].reset_index(drop=True), R)
        tc[L] = dict(collect_unit=chrom or L, collect_s=round(tb, 1), corr_s=round(time.time() - t1, 2), n=int(m.sum()))
        say(f"[C {tag}] {L}: {int(m.sum()):,} variants, corr {time.time() - t1:.2f} s")
        pd.DataFrame([dict(locus=L, **tc[L])]).to_csv(TC_PATH, mode="a", index=False,
                                                      header=not os.path.exists(TC_PATH))
    del G
    return tc


def units_for(ls):
    kept = fc.loc[ls, "kept"]
    by_chr = {}
    for L in ls:
        by_chr.setdefault(L.split("_")[0], []).append(L)
    units = []
    for c, g in by_chr.items():
        gb = kept[g].sum() * n_cols / 1e9
        units += [(c, g)] if gb <= a.collect_gb else [(None, [L]) for L in g]
    return units


def align(ls):
    ls = [L for L in ls if os.path.exists(os.path.join(OUT, "extract", L, f"{L}_{AC}_R.npz"))]
    if ls:
        subprocess.run([sys.executable, os.path.join(HERE, "align_ld_aou.py"), "--loci", ",".join(ls),
                        "--extract-dir", os.path.join(OUT, "extract"), "--out-dir", OUT, "--anc", AC], check=True)


pilot = a.pilot.split(",")
t0 = time.time()
tc = run_units(units_for(pilot), "pilot")
align(pilot)
timed("pilot_BCD", t0)

def row_corr(A, B):
    """Per-row Pearson correlation of A and B, off-diagonal entries only."""
    k = A.shape[0]
    off = ~np.eye(k, dtype=bool)
    a_, b_ = np.where(off, A, np.nan), np.where(off, B, np.nan)
    am = a_ - np.nanmean(a_, 1, keepdims=True)
    bm_ = b_ - np.nanmean(b_, 1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.nansum(am * bm_, 1) / np.sqrt(np.nansum(am ** 2, 1) * np.nansum(bm_ ** 2, 1))


def load_signed(root):
    sn = pd.read_csv(f"{root}/snps.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    v = sn["sign"].astype(float).to_numpy()
    return sn, np.load(f"{root}/R.npz")["R"].astype(np.float64) * np.outer(v, v)


def eur_z(L):
    f = os.path.expanduser(f"~/work/PGC_SCZ_sumstats/eur_sumstats/{L}_sumstats.tsv.gz")
    if not os.path.exists(f):
        return pd.Series(dtype=float)
    e = pd.read_csv(f, sep="\t", dtype=str)
    c = L.split("_")[0][3:]
    e["vkey"] = [f"{c}:{p}:{min(x.upper(), y.upper())}:{max(x.upper(), y.upper())}"
                 for p, x, y in zip(e["POS"], e["Allele1"], e["Allele2"])]
    return e.drop_duplicates("vkey").set_index("vkey")["Zscore"].astype(float)


def pilot_check_eur_gctb(pilot):
    """EUR pilot gate: PSD, and AoU EUR vs the GCTB eur_gctb LD (UKB) over shared vkeys, both signed
    to the master's ref_a1 via their own sign columns."""
    acc, flagged, all_rc = [], [], []
    for L in pilot:
        sa, Ra = load_signed(os.path.join(OUT, "ld_v2", f"{AC}_aou", L))
        min_eig = float(np.linalg.eigvalsh(Ra).min())
        r = dict(locus=L, n_aou=len(sa), min_eig=min_eig)
        gd = os.path.join(a.gctb_dir, L)
        if os.path.exists(f"{gd}/R.npz"):
            sg, Rg = load_signed(gd)
            pg = pd.Series(np.arange(len(sg)), index=sg["vkey"])
            ia = np.array([i for i, v in enumerate(sa["vkey"]) if v in pg.index], dtype=int)
            ig = pg[sa["vkey"].iloc[ia]].to_numpy()
            A, B = Ra[np.ix_(ia, ia)], Rg[np.ix_(ig, ig)]
            iu = np.triu_indices(len(ia), 1)
            rc = row_corr(A, B)
            all_rc.append(rc)
            z = eur_z(L)
            r.update(n_gctb=len(sg), n_shared=len(ia), median_abs_dR=float(np.median(np.abs(A - B)[iu])),
                     median_rowcorr=float(np.nanmedian(rc)), n_rowcorr_lt_05=int((rc < 0.5).sum()),
                     n_rowcorr_neg=int((rc < 0).sum()))
            for k in np.where(rc < 0.5)[0]:
                v = sa["vkey"].iloc[ia[k]]
                flagged.append(dict(locus=L, vkey=v, rsid=sa["rsid"].iloc[ia[k]], row_corr=float(rc[k]),
                                    eur_abs_z=abs(z.get(v, np.nan)), eur_af=sa[f"{AC}_af"].iloc[ia[k]]))
        acc.append(r)
        say("[accept] " + ", ".join(f"{k} {v:.4g}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    acc = pd.DataFrame(acc)
    acc.to_csv(os.path.join(OUT, "pilot_acceptance.csv"), index=False)
    pd.DataFrame(flagged).to_csv(os.path.join(OUT, "pilot_rowcorr_lt05.csv"), index=False)
    pooled = float(np.nanmedian(np.concatenate(all_rc))) if all_rc else np.nan
    psd_ok = bool((acc["min_eig"] >= -1e-4).all())
    say(f"[accept] PSD ok {psd_ok}; median row corr vs eur_gctb: pooled {pooled:.4f}, per locus "
        f"{', '.join(f'{x:.3f}' for x in acc.get('median_rowcorr', pd.Series(dtype=float)))}; "
        f"rows < 0.5: {len(flagged)}")
    return psd_ok and pooled >= 0.8


def pilot_check_reference(pilot):
    # Acceptance 1 against the reference tree. The reference tree used the same lenient call rate, so the
    # addendum's revised check: this run's variants are a subset of the reference's, in the same relative order,
    # with identical sign, and R equal over the shared variants to <= 1e-6.
    ok = True
    acc = []
    for L in pilot:
        rn = os.path.join(OUT, "ld_v2", f"{AC}_aou", L)
        ro = os.path.join(a.pilot_ref, L)
        sn = pd.read_csv(f"{rn}/snps.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
        so = pd.read_csv(f"{ro}/snps.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
        pos_o = {v: i for i, v in enumerate(so["vkey"])}
        extra = [v for v in sn["vkey"] if v not in pos_o]
        io = np.array([pos_o[v] for v in sn["vkey"] if v in pos_o], dtype=int)
        subset = not extra
        order = subset and bool(np.all(np.diff(io) > 0))
        same_sign = subset and sn["sign"].tolist() == so["sign"].iloc[io].tolist()
        d = (np.abs(np.load(f"{rn}/R.npz")["R"].astype(np.float64)
                    - np.load(f"{ro}/R.npz")["R"].astype(np.float64)[np.ix_(io, io)]).max() if subset else np.inf)
        acc.append(dict(locus=L, n_b43=len(so), n_b45=len(sn), dropped_vs_b43=len(so) - len(io), extra_vs_b43=len(extra),
                        subset=subset, same_order=order, same_sign=same_sign, max_abs_delta=d))
        say(f"[accept] {L}: new {len(sn)} vs reference {len(so)} (dropped {len(so) - len(io)}, extra {len(extra)}), "
            f"subset {subset}, order {order}, same sign {same_sign}, max|delta| {d:.3g}")
        ok &= subset and order and same_sign and d <= 1e-6
    pd.DataFrame(acc).to_csv(os.path.join(OUT, "pilot_acceptance.csv"), index=False)
    return ok


ok = pilot_check_reference(pilot) if a.pilot_check == "reference" else pilot_check_eur_gctb(pilot)
if not ok:
    say("[accept] FAIL: stopping before the other loci")
    sys.exit(2)
say("[accept] PASS at all pilot loci")

if a.pilot_only:
    say("pilot only: done")
    sys.exit(0)
rest = [L for L in loci if L not in set(pilot)]
t0 = time.time()
tc.update(run_units(units_for(rest), "rest"))
timed("rest_BC", t0)
t0 = time.time()
align(rest)
timed("rest_D", t0)
if FAILED:
    say(f"[FAIL] {len(FAILED)} loci not built: {','.join(FAILED)}")
say(f"done; cluster: {cluster_shape()}")
