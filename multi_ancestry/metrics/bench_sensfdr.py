#!/usr/bin/env python3
"""Sensitivity / FDR of fine-mapping runs against a fixed published benchmark.

Design settled with Qin 2026-09-21. The three decisions that shape it:

  * **No noise floor.** FDR(t) = |S_C(t) \\ P_B| / |S_C(t)| -- every call outside the
    benchmark's positive set is a false discovery. This is exactly 1 - precision, and it
    is what removes the undecided middle band that made the previous FDR curve
    non-monotone for reasons that had nothing to do with the methods (§28).
  * **Sensitivity is capped by `max_sensitivity`, and must be read against it.** Because s
    is fixed and t sweeps, at t > s the benchmark's own positive set includes variants the
    benchmark does not call at t either. At t = 0.9 only 13 of the 40 positives have
    benchmark PIP > 0.9, so 0.325 is a perfect score, not a poor one.

  * **s is FIXED, only t sweeps.** P_B is computed once. Letting s track t re-labels the
    reference at every point, which makes sensitivity rise at high t simply because the
    benchmark's surviving positives are its easiest calls -- the curve then is not a
    trajectory (§28, and the matched-sweep tables in metrics/README.md).
  * **Loci where the benchmark produced no credible set are excluded.** 57 of 250 have
    every PIP at 0 because nothing survived --min_purity=0.5, not because SuSiEx judged
    the locus empty. Keeping them makes every comparator call there a false discovery by
    construction (§29).

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

Usage:
  python3 bench_sensfdr.py --smoke-locus chr12_110129434_111489434     # inspect first
  python3 bench_sensfdr.py --out-prefix bench            # full run
  python3 bench_sensfdr.py --out-prefix bench --n-boot 0  # point estimates only, faster
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

from sensitivity_fdr import load_pips, variant_key

BENCH = "pip_vectors_susiex_published_full.csv.gz"
# Yuan et al. raw SuSiEx output; local only (git-ignored), available from the authors on request
KAI = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "inputs", "202.EUR_ASN.SuSiEx_spec")
MIN_CALLS = 20          # below this the ratio is unstable; flagged, not dropped

# (label, file, method-in-file, L, rho, ld_source)
def run_catalog():
    R = []
    for L in (1, 3, 5, 10):
        for m in ("SuSiEx-ours", "MultiSuSiE", "MESuSiE", "SuSiE-merged", "FINEMAP-merged"):
            name = f"{m}_L{L}_nomultistep" if m == "SuSiEx-ours" else f"{m}_L{L}"
            R.append((m, "pip_vectors_L_sweep_v2.csv.gz", name, L, np.nan, "merged/per-anc"))
    for rho in ("0", "0.75", "0.99"):
        R.append(("MultiSuSiE", "../aou_ld/task_b1_out/rho_sweep_pip_vectors.csv.gz",
                  f"MultiSuSiE_rho{rho}", 10, float(rho), "merged/per-anc"))
    for tag, f in (("gctb10k", "eur_susie_gctb10k_pip_vectors.csv.gz"),
                   ("polyfun337k", "eur_susie_polyfun337k_pip_vectors.csv.gz")):
        R.append(("EUR-SuSiE-RSS", f"../aou_ld/task_b22_out/{f}",
                  f"EUR-SuSiE-RSS_{tag}", 5, np.nan, tag))
    # RSparsePro (Task B27 Job 1). bp/a1/a2 backfilled from the tool's own rsid tag by
    # aou_ld/parse_rsparsepro_keys_b27.py -- the raw exports carry no bp/a1/a2 at all, which
    # made variant_key's rsid fallback compare a "chr10:105040690:C:T"-style tag against every
    # other arm's "105040690:C:T" allele key and match nothing (0/231,272 verified before the
    # fix). RSparsePro-merged and RSparsePro-EUR are NOT on matched inputs -- B16 rebuilt only
    # the merged arm on the corrected finemap_in_rebuilt; RSparsePro-EUR still reads the
    # original native single-ancestry build. Both the original merged export and the rebuilt
    # one are included below, unreduced, so this does not silently pick one for Qin.
    R.append(("RSparsePro", "../aou_ld/task_b27_out/rsparsepro_pip_vectors_keyed.csv.gz",
              "RSparsePro-merged", 10, np.nan, "merged"))
    R.append(("RSparsePro", "../aou_ld/task_b27_out/rsparsepro_merged_rebuilt_pip_vectors_keyed.csv.gz",
              "RSparsePro-merged-rebuilt", 10, np.nan, "merged"))
    R.append(("RSparsePro", "../aou_ld/task_b27_out/rsparsepro_pip_vectors_keyed.csv.gz",
              "RSparsePro-EUR", 10, np.nan, "EUR-only"))
    return R


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


def comparator_cs():
    """Credible-set membership for our arms, where B12 exported it."""
    out = []
    for f in glob.glob("../aou_ld/task_b12_out/*_cs_purity.csv.gz"):
        d = pd.read_csv(f, low_memory=True,
                        usecols=["method", "locus", "bp", "a1", "a2", "rsid", "cs_ids_all"])
        d = d[d.cs_ids_all.notna()]
        if len(d):
            d["cs"] = d.locus + "|C" + d.cs_ids_all.astype(str)
            out.append(d[["method", "locus", "bp", "a1", "a2", "rsid", "cs"]])
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s", type=float, default=0.5, help="FIXED benchmark threshold")
    ap.add_argument("--step", type=float, default=0.01)
    ap.add_argument("--tmin", type=float, default=0.05)
    ap.add_argument("--tmax", type=float, default=0.95)
    ap.add_argument("--kai", default=KAI)
    ap.add_argument("--n-boot", type=int, default=2000,
                    help="bootstrap replicates over LOCI for the CIs; 0 disables")
    ap.add_argument("--ci", type=float, default=0.95)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--loci", default=None,
                    help="file of locus ids, one per line, to restrict to (e.g. "
                         "loci_single175.txt). Intersected with the benchmark's own "
                         "has-a-credible-set set, never a superset of it")
    ap.add_argument("--smoke-locus", default=None)
    ap.add_argument("--out-prefix", default="bench")
    args = ap.parse_args()

    thresholds = np.round(np.arange(args.tmin, args.tmax + 1e-9, args.step), 4)
    rng = np.random.default_rng(args.seed)
    keep_loci, bcs = bench_cs_loci(args.kai)
    print(f"[info] benchmark loci with >=1 credible set: {len(keep_loci)} "
          f"(excluded {250 - len(keep_loci)} purity-empty)")
    if args.loci:
        want = set(open(args.loci).read().split())
        keep_loci = keep_loci & want
        print(f"[info] restricted to {args.loci}: {len(want)} requested, "
              f"{len(keep_loci)} kept after intersecting with the above")

    b = load_pips(BENCH)
    b = b[b.locus.isin(keep_loci)][["locus", "key", "bp", "pip"]].rename(columns={"pip": "pb"})
    b = b.groupby(["locus", "key"], as_index=False).agg(pb=("pb", "max"), bp=("bp", "first"))

    ccs = comparator_cs()
    if len(ccs):
        ccs["key"] = variant_key(ccs)
    bcs = bcs.merge(b[["locus", "bp", "key"]], on=["locus", "bp"], how="inner")
    bench_map = bcs[["locus", "key", "cs"]].drop_duplicates(["locus", "key"])
    comp_map = ({m: g[["locus", "key", "cs"]].drop_duplicates(["locus", "key"])
                 for m, g in ccs.groupby("method")} if len(ccs) else {})

    cache, var_rows, clu_rows, smoke = {}, [], [], None
    for meth, f, name, L, rho, ld in run_catalog():
        if f not in cache:
            cache[f] = load_pips(f)
        c = cache[f]
        c = c[c.method == name]
        if c.empty:
            print(f"[warn] {name} not in {f} -- skipped")
            continue
        c = c[c.locus.isin(keep_loci)][["locus", "key", "pip"]].rename(columns={"pip": "pc"})
        c = c.groupby(["locus", "key"], as_index=False).pc.max()
        shared = set(b.locus) & set(c.locus)
        u = b[b.locus.isin(shared)].merge(c[c.locus.isin(shared)], on=["locus", "key"],
                                          how="inner")
        if args.smoke_locus:
            u = u[u.locus == args.smoke_locus]
            if u.empty:
                continue

        cv, P = curves(u, args.s, thresholds, n_boot=args.n_boot,
                       ci=args.ci, rng=rng)
        cv.insert(0, "ld_source", ld); cv.insert(0, "rho", rho)
        cv.insert(0, "L", L); cv.insert(0, "run", name); cv.insert(0, "method", meth)
        var_rows.append(cv)

        # ---- cluster level: one group per credible set from EITHER side ----
        # merge, never iterrows: this runs once per run and the frames are ~200k rows
        gb = bench_map                       # benchmark CS, run-independent
        gc = comp_map.get(meth)              # this arm's CS, if B12 exported it
        g = u[["locus", "key"]].merge(gb, on=["locus", "key"], how="left")
        if gc is not None:
            g = g.merge(gc, on=["locus", "key"], how="left", suffixes=("", "_c"))
            g["cs"] = g.cs.fillna(g.cs_c)
        solo = pd.Series("solo|" + u.locus.values + "|" + u.key.values, index=g.index)
        g = g.cs.fillna(solo).values
        cl, Pc = curves(u, args.s, thresholds, group=g,
                        n_boot=args.n_boot, ci=args.ci, rng=rng)
        # swap rate: exact-variant misses that are cluster-level hits, at t = s
        S = u.pc > args.s
        miss = P & ~S
        gm = pd.Series(g, index=u.index)
        hit_groups = set(gm[S & P]) | set(gm[u.pc > args.s])
        swap = float(gm[miss].isin(hit_groups).mean()) if int(miss.sum()) else np.nan
        cl.insert(0, "swap_rate_at_s", swap)
        cl.insert(0, "ld_source", ld); cl.insert(0, "rho", rho)
        cl.insert(0, "L", L); cl.insert(0, "run", name); cl.insert(0, "method", meth)
        clu_rows.append(cl)

        if smoke is None:
            smoke = (name, u, P, cv)

    var = pd.concat(var_rows, ignore_index=True)
    clu = pd.concat(clu_rows, ignore_index=True)

    if args.smoke_locus:
        name, u, P, cv = smoke
        print(f"\n===== SMOKE TEST — locus {args.smoke_locus}, run {name} =====")
        print(f"|U| = {len(u)}   |P_B| (benchmark PIP > {args.s}) = {int(P.sum())}")
        print("\nP_B members:")
        print(u[P][["locus", "key", "bp", "pb", "pc"]].to_string(index=False))
        print(f"\noutput rows at t = 0.10 / 0.50 / 0.90:")
        print(cv[cv.threshold.isin([0.10, 0.50, 0.90])].to_string(index=False))
        return

    var.to_csv(f"{args.out_prefix}_variant_level.csv", index=False)
    clu.to_csv(f"{args.out_prefix}_cluster_level.csv", index=False)
    print(f"[done] {args.out_prefix}_variant_level.csv  ({len(var):,} rows, "
          f"{var.run.nunique()} runs)")
    print(f"[done] {args.out_prefix}_cluster_level.csv  ({len(clu):,} rows)")


if __name__ == "__main__":
    main()
