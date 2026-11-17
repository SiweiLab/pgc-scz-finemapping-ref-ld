#!/usr/bin/env python3
"""Sensitivity / FDR of fine-mapping runs against a fixed published benchmark.

Design settled 2026-09-21. The three decisions that shape it:

  * **No noise floor.** FDR(t) = |S_C(t) \\ P_B| / |S_C(t)| -- every call outside the
    benchmark's positive set is a false discovery. This is exactly 1 - precision, and it
    is what removes the undecided middle band that made the previous FDR curve
    non-monotone for reasons that had nothing to do with the methods.
  * **Sensitivity is capped by `max_sensitivity`, and must be read against it.** Because s
    is fixed and t sweeps, at t > s the benchmark's own positive set includes variants the
    benchmark does not call at t either. At t = 0.9 only 13 of the 40 positives have
    benchmark PIP > 0.9, so 0.325 is a perfect score, not a poor one.

  * **s is FIXED, only t sweeps.** P_B is computed once. Letting s track t re-labels the
    reference at every point, which makes sensitivity rise at high t simply because the
    benchmark's surviving positives are its easiest calls -- the curve then is not a
    trajectory.
  * **Loci where the benchmark produced no credible set are excluded.** 57 of 250 have
    every PIP at 0 because nothing survived --min_purity=0.5, not because SuSiEx judged
    the locus empty. Keeping them makes every comparator call there a false discovery by
    construction.

  * **Error bars bootstrap LOCI, not variants.** |P_B| is only ~40 variants, so the point
    estimates are noisy and a bar is not optional. Variants inside a locus are in LD and
    nowhere near independent, so a binomial interval on the variant count would be far too
    narrow; we resample the 193 loci with replacement (2000 draws) and take percentiles.
    The benchmark's positive set is resampled along with them, which is deliberate.

Joining uses the metric layer's `variant_key` -- (bp, unordered allele pair) with an rsID
fallback -- because the per-ancestry arms carry dbSNP rsIDs while the merged-LD arms carry
`chr:pos:a1:a2` tags, and a join on rsid alone silently returns nothing.

Cluster level groups variants by 95% credible-set membership from EITHER side, so a call
that lands on a different variant in the same set still counts. Cluster PIP is
1 - prod(1 - PIP) within the group; the swap rate is the fraction of exact-variant misses
that are cluster-level hits.

A library: the figure scripts import BENCH, BENCH_RAW, bench_cs_loci and curves from it.
"""
import glob
import os
import re

import numpy as np
import pandas as pd

from sensitivity_fdr import load_pips, variant_key

BENCH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "benchmark",
                     "pip_vectors_susiex_published_full.csv.gz")
# Yuan et al. raw SuSiEx output; local only (git-ignored), available from the authors on request
BENCH_RAW = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "inputs", "external", "202.EUR_ASN.SuSiEx_spec")
MIN_CALLS = 20          # below this the ratio is unstable; flagged, not dropped

# (label, file, method-in-file, L, rho, ld_source)
def bench_cs_loci(src):
    """Loci where the benchmark produced at least one credible set, and its CS membership."""
    keep, mem = set(), []
    for snp in sorted(glob.glob(os.path.join(os.path.expanduser(src), "*.cs99.snp"))):
        stem = snp[:-4]
        first = open(stem + ".summary").readline().strip()
        m = re.match(r"#\s*chr(\w+):(\d+)-(\d+)", first)
        lid = f"chr{m.group(1)}_{m.group(2)}_{m.group(3)}"
        lid = {"chr15_34646330_34663925": "chr15_34645987_34664593"}.get(lid, lid)
        cs = stem + ".cs"
        head = open(cs).readline()
        if head.startswith("NULL"):
            continue                      # no set survived purity -- exclude the locus
        keep.add(lid)
        d = pd.read_csv(cs, sep="\t")
        mem.append(pd.DataFrame({"locus": lid, "bp": d.BP,
                                 "cs": lid + "|B" + d.CS_ID.astype(str)}))
    return keep, (pd.concat(mem, ignore_index=True) if mem else pd.DataFrame())


def _per_locus(mat, codes_order, bounds):
    """Sum a [n x T] (or [n]) array within each locus. reduceat, not add.at: this runs
    once per array per run and add.at on 140k x 86 is minutes, not seconds."""
    return np.add.reduceat(mat[codes_order], bounds, axis=0)


def _boot_ci(num, den, W, ci):
    """Percentile CI of a ratio-of-sums statistic under a bootstrap over LOCI.

    `num`/`den` are per-locus [nL x T] (or [nL]) totals; W is the [B x nL] matrix of
    resample weights. Resampling whole loci -- not variants -- is the point: variants
    inside a locus are in LD and nowhere near independent, so a binomial/Wald interval
    on the variant count would be far too narrow. It also means the benchmark's positive
    set is itself resampled, so |P_B| varies across replicates; that uncertainty is real
    and belongs in the bar, given |P_B| is only ~40.
    """
    n = W @ num
    d = W @ den
    if d.ndim == 1:
        d = d[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(d > 0, n / np.where(d > 0, d, 1), np.nan)
    lo = np.nanpercentile(r, 100 * (1 - ci) / 2, axis=0)
    hi = np.nanpercentile(r, 100 * (1 + ci) / 2, axis=0)
    se = np.nanstd(r, axis=0, ddof=1)
    return lo, hi, se


def curves(u, s, thresholds, group=None, n_boot=2000, ci=0.95, rng=None):
    """Hard and expected sensitivity/FDR. `group` -> cluster level, else variant level.

    Adds a locus-level bootstrap CI on sensitivity and FDR, and `predicted_fdr` --
    1 - mean(own PIP) over the calls, which is MultiSuSiE's benchmark-free expected FDR
    (sec 30). `expected_fdr` is the same shape computed from the BENCHMARK's PIPs; the two
    answer different questions and both are kept.
    """
    d = u[["locus", "pb", "pc"]].copy()
    if group is not None:
        d = (d.assign(g=group)
               .groupby("g", as_index=False)
               .agg(locus=("locus", "first"),
                    pb=("pb", lambda v: 1 - np.prod(1 - v)),
                    pc=("pc", lambda v: 1 - np.prod(1 - v))))
    pb, pc = d.pb.to_numpy(), d.pc.to_numpy()
    P = pb > s
    thr = np.asarray(thresholds, dtype=float)
    S = pc[:, None] > thr[None, :]                       # [n x T]

    codes, loci = pd.factorize(d.locus, sort=True)
    order = np.argsort(codes, kind="stable")
    bounds = np.searchsorted(codes[order], np.arange(len(loci)))
    NS = _per_locus(S.astype(np.int64), order, bounds)                   # [nL x T]
    TP = _per_locus((S & P[:, None]).astype(np.int64), order, bounds)
    NP = _per_locus(P.astype(np.int64), order, bounds)                   # [nL]
    SPC = _per_locus(np.where(S, pc[:, None], 0.0), order, bounds)
    SPB = _per_locus(np.where(S, pb[:, None], 0.0), order, bounds)
    MPB = _per_locus(pb, order, bounds)

    nS, tp = NS.sum(0), TP.sum(0)
    n_P, mass_U = int(NP.sum()), float(MPB.sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        sens = tp / n_P if n_P else np.full_like(thr, np.nan)
        fdr = np.where(nS > 0, (nS - tp) / np.where(nS > 0, nS, 1), np.nan)
        pred = np.where(nS > 0, 1 - SPC.sum(0) / np.where(nS > 0, nS, 1), np.nan)
        exp_s = SPB.sum(0) / mass_U if mass_U else np.full_like(thr, np.nan)
        exp_f = np.where(nS > 0, 1 - SPB.sum(0) / np.where(nS > 0, nS, 1), np.nan)
        # Ceiling on sensitivity: the fraction of P_B the BENCHMARK ITSELF puts above t.
        # s is fixed at 0.5 while t sweeps (sec 28), so at t > s the reference set contains
        # positives the benchmark does not itself call at t -- a comparator that reproduced
        # the benchmark exactly would still score only this. Without it a bar at t = 0.9
        # reads as "found 30% of the signal" when it is in fact at the ceiling.
        max_s = (pb[P][:, None] > thr[None, :]).mean(axis=0) if n_P else np.full_like(thr, np.nan)

    if n_boot and len(loci) > 1:
        rng = rng or np.random.default_rng(0)
        W = rng.multinomial(len(loci), np.full(len(loci), 1 / len(loci)),
                            size=n_boot).astype(np.float64)
        s_lo, s_hi, s_se = _boot_ci(TP, NP, W, ci)
        f_lo, f_hi, f_se = _boot_ci(NS - TP, NS, W, ci)
    else:
        nan = np.full_like(thr, np.nan)
        s_lo = s_hi = s_se = f_lo = f_hi = f_se = nan

    out = pd.DataFrame(dict(
        threshold=np.round(thr, 4), n_U=len(d), n_loci=len(loci), n_PB=n_P,
        n_SC=nS, n_true_pos=tp,
        sensitivity=sens, sensitivity_lo=s_lo, sensitivity_hi=s_hi, sensitivity_se=s_se,
        fdr=fdr, fdr_lo=f_lo, fdr_hi=f_hi, fdr_se=f_se,
        predicted_fdr=pred,
        expected_sensitivity=exp_s, expected_fdr=exp_f, max_sensitivity=max_s,
        n_boot=int(n_boot), ci=ci,
        unstable=nS < MIN_CALLS,
    ))
    return out, pd.Series(P, index=d.index)
