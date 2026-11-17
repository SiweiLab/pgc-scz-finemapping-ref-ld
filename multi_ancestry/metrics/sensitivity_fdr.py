#!/usr/bin/env python3
"""Credible-set-free sensitivity/FDR curves against a published benchmark.

Implements the metric requested in co-author review:

  sensitivity(t) = of variants the BENCHMARK puts above t, the fraction the
                   comparator also puts above t
  FDR(t)         = of variants the COMPARATOR puts above t, the fraction the
                   benchmark calls noise (below --noise-threshold)

Both are variant-exact and credible-set-free, swept over t in small steps, so a
locus where a method reports no surviving credible set still contributes its
PIPs (chr4/rs35734242 is exactly that case). This is the point of the metric:
the credible-set-based concordance it replaces cannot see those loci at all.

Identifier handling (do not simplify this away -- see PROGRESS_LOG.md 2026-09-08):
the three per-ancestry arms key variants by dbSNP rsID while the two merged-LD
arms key them by a `chr:pos:a1:a2` tag, so a join on `rsid` alone silently
fails across arms. We join on (bp, {a1,a2}) as an unordered allele pair, with
an rsID fallback for the 4.8% of SuSiEx rows carrying no allele info.

Benchmark absence is NOT missingness: the published tables report only variants
above their own PIP floor, so a variant absent from the benchmark is below that
floor. It is scored as benchmark PIP 0, which is why --noise-threshold must be
>= the floor (--benchmark-floor) to be computable from public data.

Usage:
  python3 sensitivity_fdr.py \
      --pip-vectors ../aou_ld/task_a_out/pip_vectors_full250_qcexcl.csv.gz \
      --benchmark   ../aou_ld/task_a_out/published_susiex_pip_vectors.csv.gz \
      --loci        <optional file, one locus id per line, to restrict stratum> \
      --out         curves.csv

  # validate the machinery with no benchmark file, by holding out one arm:
  python3 sensitivity_fdr.py --pip-vectors ... --self-test SuSiEx-ours
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

PER_ANCESTRY = ["MultiSuSiE", "MESuSiE", "SuSiEx-ours"]
MERGED = ["SuSiE-merged", "FINEMAP-merged"]


_MASTER_ALLELES = None


def _backfill_alleles(df):
    """Fill blank a1/a2 from the v2 master by (locus, rsid), PGC3's own IDs (2026-10-01).

    Our copy of Yuan's posterior (reconstruct_published_full_posterior.py) has no alleles of its own: they were backfilled
    from v1's PIP vectors, so every variant v1 lacked (37,687 rows, 28,361 at the 193 benchmark loci, incl. 2 positives
    and all 293 multi-allelic ones) stayed blank and fell back to an rsid key that never joins an allele-keyed export.
    v1's own exports have no blank alleles, so their keys are unchanged."""
    global _MASTER_ALLELES
    if not {"locus", "rsid"} <= set(df.columns):
        return df
    blank = (df["a1"].isna() | df["a2"].isna()).to_numpy()
    if not blank.any():
        return df
    if _MASTER_ALLELES is None:
        mp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "aou_ld", "v2_inputs",
                          "master_variants_250.tsv.gz")
        m = pd.read_csv(mp, sep="\t", dtype=str, keep_default_na=False, usecols=["locus", "rsid", "ref_a1", "ref_a2"])
        _MASTER_ALLELES = {(l, r): (x, y) for l, r, x, y in zip(m.locus, m.rsid, m.ref_a1, m.ref_a2)}
    df = df.copy()
    df["a1"] = df["a1"].astype(object)
    df["a2"] = df["a2"].astype(object)
    for i in np.flatnonzero(blank):
        hit = _MASTER_ALLELES.get((df["locus"].iat[i], str(df["rsid"].iat[i])))
        if hit:
            df.iat[i, df.columns.get_loc("a1")], df.iat[i, df.columns.get_loc("a2")] = hit
    left = int((df["a1"].isna() | df["a2"].isna()).sum())
    if left:
        print(f"[warn] variant_key: {left:,} rows still have no alleles and are keyed by bare rsid", file=sys.stderr)
    return df


def variant_key(df):
    """(bp, unordered allele pair) where alleles exist, else the raw rsid.

    Blank alleles are first filled from the master by (locus, rsid) (_backfill_alleles), so an rsid key is now the
    exception and is reported.

    Unordered because a1/a2 orientation is not guaranteed to agree across arms;
    verified on chr12 that ordered and unordered joins give the same overlap, but
    unordered is the safe default for the other 249 loci, which were not checked.
    """
    df = _backfill_alleles(df)
    a1 = df["a1"].astype("string").str.upper()
    a2 = df["a2"].astype("string").str.upper()
    have = (a1.notna() & a2.notna()).to_numpy(dtype=bool)
    # fill before comparing: pandas NA in a <= b raises rather than returning False
    a1f = a1.fillna("").to_numpy(dtype=object)
    a2f = a2.fillna("").to_numpy(dtype=object)
    swap = a1f <= a2f
    lo = np.where(swap, a1f, a2f)
    hi = np.where(swap, a2f, a1f)
    allele_key = (df["bp"].astype("Int64").astype(str).to_numpy(dtype=object)
                  + ":" + lo + ":" + hi)
    return pd.Series(np.where(have, allele_key, df["rsid"].astype(str).to_numpy()),
                     index=df.index)


def load_pips(path):
    d = pd.read_csv(path)
    missing = {"method", "locus", "rsid", "bp", "pip"} - set(d.columns)
    if missing:
        sys.exit(f"[error] {path} missing columns: {sorted(missing)}")
    for c in ("a1", "a2"):
        if c not in d.columns:
            d[c] = pd.NA
    d["key"] = variant_key(d)
    return d





def top_variant_only(d):
    """Keep just the highest-PIP variant per (method, locus).

    Makes the sweep rank-sensitive. The default metric asks "did you put mass above t on
    this variant", which a variant sitting at rank 2 satisfies as well as one at rank 1 --
    so it cannot see the difference the manuscript's argmax comparison reports. Restricted
    to each source's own top variant, a locus counts only if the method's *winner* is the
    benchmark's winner, which is the manuscript's question asked at every threshold.

    Ties at the maximum are broken by first occurrence; exact ties are rare and a tie means
    the method did not resolve a winner anyway.
    """
    idx = d.groupby(["method", "locus"])["pip"].idxmax()
    return d.loc[idx].reset_index(drop=True)


def wilson(k, n, z=1.96):
    """Wilson score interval for a proportion.

    Used rather than the normal approximation because these counts are small
    (tens of variants) and the proportions sit near 0 and 1, where the normal
    interval is badly calibrated and can run outside [0, 1].
    """
    if not n:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z / d * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, c - h), min(1.0, c + h))


def load_proxy_map(path):
    """(locus, key) -> set of keys in near-perfect LD with it.

    Expected columns: locus, rsid_a, bp_a, a1_a, a2_a, rsid_b, bp_b, a1_b, a2_b, r2
    -- keyed the same way as the PIP tables so the two join. Pairs are symmetric;
    the file need only carry each once.
    """
    p = pd.read_csv(path)
    need = {"locus", "rsid_a", "bp_a", "rsid_b", "bp_b", "r2"}
    missing = need - set(p.columns)
    if missing:
        sys.exit(f"[error] {path} missing columns: {sorted(missing)}")
    for side in ("a", "b"):
        sub = p[[f"rsid_{side}", f"bp_{side}", f"a1_{side}", f"a2_{side}"]].copy()
        sub.columns = ["rsid", "bp", "a1", "a2"]
        p[f"key_{side}"] = variant_key(sub).values
    m = {}
    for locus, ka, kb in zip(p["locus"], p["key_a"], p["key_b"]):
        m.setdefault((locus, ka), set()).add(kb)
        m.setdefault((locus, kb), set()).add(ka)
    print(f"[info] proxy map: {len(p)} pairs, {len(m)} variants with >=1 proxy")
    return m


def _best_over_proxies(j, value_col, proxy_map):
    """Max of value_col over the variant itself and its near-perfect proxies.

    A benchmark variant counts as recovered if the comparator commits to it OR to
    anything in near-perfect LD with it; symmetrically a comparator hit is not a
    false discovery if it sits in near-perfect LD with a real benchmark signal.
    Without this, two variants at r2=0.995 score as a total miss AND a false
    discovery -- which is exactly the chr1 rs11587347/rs72769124 case in the
    manuscript, and why exact matching understates every method.
    """
    lookup = dict(zip(zip(j["locus"], j["key"]), j[value_col]))
    out = []
    for locus, key, own in zip(j["locus"], j["key"], j[value_col]):
        best = own
        for pk in proxy_map.get((locus, key), ()):
            v = lookup.get((locus, pk))
            if v is not None and v > best:
                best = v
        out.append(best)
    return np.asarray(out)


def curve(bench, comp, thresholds, noise_threshold, bench_threshold=None, proxy_map=None,
          restrict_to_shared=True):
    """One method's sensitivity/FDR curve.

    bench, comp: DataFrames with columns (locus, key, pip), already restricted to
    the locus stratum. Joined on (locus, key) -- never on rsid alone.

    Two modes, and the choice matters (see PROGRESS_LOG / plan sec 1.3):

    bench_threshold=None -- "matched sweep", the requested spec: the SAME t gates
      both sides. The benchmark's positive set therefore SHRINKS as t rises, so
      the sensitivity denominator moves with the numerator and the parametric
      curve is NOT monotone -- it doubles back on itself and is hard to read.

    bench_threshold=B -- "fixed benchmark": the positive set is pinned at
      benchmark PIP > B once, then only the comparator's t is swept. Both
      sensitivity and FDR then fall monotonically with t, giving the clean
      precision-recall shape the figure wants. This is also the shape the MACHINE
      preprint uses for its replication-failure rate (0.9 to identify, 0.1 to
      call a miss).
    """
    b = bench.rename(columns={"pip": "pip_bench"})[["locus", "key", "pip_bench"]]
    c = comp.rename(columns={"pip": "pip_comp"})[["locus", "key", "pip_comp"]]
    j = b.merge(c, on=["locus", "key"], how="outer", indicator=True)
    # presence on the benchmark side, captured BEFORE fillna(0) below erases the
    # difference between "the benchmark scored this at 0" and "the benchmark never saw it"
    j["fitted"] = j["_merge"].isin(["both", "left_only"])

    # A join with no overlap is always a key-space bug, never a result, but every
    # downstream quantity is still well-defined on it -- sensitivity 0, FDR 1, at every
    # threshold -- so it produces a plottable curve instead of an error. That happened
    # 2026-09-12: an export written without a1/a2 fell back to raw rsIDs (see
    # variant_key) while the benchmark keyed on bp:allele-pair, and the figure was built
    # before anyone noticed. Fail loudly instead.
    n_both = int((j["_merge"] == "both").sum())
    if n_both == 0:
        bk = b["key"].head(3).tolist()
        ck = c["key"].head(3).tolist()
        sys.exit(
            f"[error] benchmark and comparator share ZERO variants after the join.\n"
            f"        benchmark keys look like: {bk}\n"
            f"        comparator keys look like: {ck}\n"
            f"        variant_key uses (bp, unordered allele pair) where a1/a2 are present "
            f"and falls back to the raw rsid where they are not, so a file missing alleles "
            f"cannot join to one carrying them. Add a1/a2 to whichever side lacks them.")
    smaller = min(len(b), len(c))
    if n_both < 0.05 * smaller:
        print(f"[warn] only {n_both:,} variants join, against {smaller:,} on the smaller "
              f"side ({100*n_both/smaller:.1f}%). That is low enough to suggest a partial "
              f"key-space mismatch rather than a genuine variant-set difference -- check "
              f"a1/a2 coverage on both sides before trusting these curves.", file=sys.stderr)
    j = j.drop(columns="_merge")

    # absence on either side means "below that source's reporting floor", not unknown
    j["pip_bench"] = j["pip_bench"].fillna(0.0)
    j["pip_comp"] = j["pip_comp"].fillna(0.0)

    if proxy_map:
        # comparator evidence reachable from each benchmark variant, and vice versa
        j["pip_comp_nbr"] = _best_over_proxies(j, "pip_comp", proxy_map)
        j["pip_bench_nbr"] = _best_over_proxies(j, "pip_bench", proxy_map)
    else:
        j["pip_comp_nbr"] = j["pip_comp"]
        j["pip_bench_nbr"] = j["pip_bench"]

    rows = []
    # did the benchmark fit this variant at all? (recorded before the fillna above makes
    # "absent" and "zero" indistinguishable -- see the false-discovery note below)
    fixed_pos = j["pip_bench"] > bench_threshold if bench_threshold is not None else None
    for t in thresholds:
        bench_pos = fixed_pos if fixed_pos is not None else (j["pip_bench"] > t)
        comp_pos = j["pip_comp"] > t
        n_bench = int(bench_pos.sum())
        n_comp = int(comp_pos.sum())
        # recovered: comparator commits at the variant OR at a near-perfect proxy
        recovered = int((bench_pos & (j["pip_comp_nbr"] > t)).sum())
        # false discovery: no real benchmark signal at the variant OR any proxy.
        #
        # A variant the benchmark never FITTED is not a false discovery -- there is no
        # ground truth at it to disagree with, so it lies outside the comparison domain
        # rather than inside it as an error. Counting those (the default before
        # 2026-09-14) inflates FDR by 0.07-0.13 depending on the arm, and by different
        # amounts per arm, because the arms differ in how much of their own variant set
        # the benchmark shares. `fitted` marks presence in the benchmark's own vector.
        if restrict_to_shared:
            scoreable = comp_pos & j["fitted"]
            n_comp = int(scoreable.sum())
            false_disc = int((scoreable & (j["pip_bench_nbr"] < noise_threshold)).sum())
        else:
            false_disc = int((comp_pos & (j["pip_bench_nbr"] < noise_threshold)).sum())
        rows.append({
            "threshold": round(float(t), 4),
            "n_benchmark_positive": n_bench,
            "n_comparator_positive": n_comp,
            "n_recovered": recovered,
            "sensitivity": recovered / n_bench if n_bench else np.nan,
            "n_false_discoveries": false_disc,
            "fdr": false_disc / n_comp if n_comp else np.nan,
            "sensitivity_lo": wilson(recovered, n_bench)[0],
            "sensitivity_hi": wilson(recovered, n_bench)[1],
            "fdr_lo": wilson(false_disc, n_comp)[0],
            "fdr_hi": wilson(false_disc, n_comp)[1],
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pip-vectors", required=True)
    ap.add_argument("--benchmark", help="tidy published-benchmark PIP vectors "
                                        "(same schema); omit only with --self-test")
    ap.add_argument("--benchmark-method",
                    help="which source in --benchmark to score against, e.g. "
                         "SuSiEx-published. Required when the file holds more than one.")
    ap.add_argument("--self-test", metavar="ARM",
                    help="use this arm as a stand-in benchmark to exercise the code "
                         "before the published table is available")
    ap.add_argument("--loci", help="file of locus ids, one per line, to restrict the stratum")
    ap.add_argument("--noise-threshold", type=float, default=0.1,
                    help="benchmark PIP below which a comparator hit counts as a false "
                         "discovery (default 0.1)")
    ap.add_argument("--benchmark-floor", type=float, default=0.1,
                    help="the published table's own PIP reporting floor; "
                         "--noise-threshold below this is not computable (default 0.1)")
    ap.add_argument("--benchmark-threshold", type=float, default=None,
                    help="pin the benchmark's positive set at this PIP and sweep only the "
                         "comparator (recommended, gives monotone curves). Omit for the "
                         "literal matched sweep, where both sides move together.")
    ap.add_argument("--no-restrict-to-benchmark-loci", action="store_true",
                    help="score every locus, including those the benchmark never covers. "
                         "Almost always wrong: comparator hits at an uncovered locus become "
                         "false discoveries by construction.")
    ap.add_argument("--count-unfitted-as-false", action="store_true",
                    help="count comparator calls at variants the benchmark never fitted as "
                         "false discoveries. This was the behaviour before 2026-09-14 and is "
                         "almost certainly wrong -- there is no ground truth at such a variant "
                         "to disagree with. Kept only to reproduce earlier numbers.")
    ap.add_argument("--proxy-map", help="high-LD variant pairs (r2 above the study's "
                                        "0.95 cutoff); enables LD-aware scoring")
    ap.add_argument("--top-variant-only", action="store_true",
                    help="restrict both sides to each locus's single highest-PIP "
                         "variant, making the sweep rank-sensitive")
    ap.add_argument("--min-threshold", type=float, default=None,
                    help="lowest threshold to sweep (default: the benchmark floor, "
                         "below which sensitivity is undefined)")
    ap.add_argument("--step", type=float, default=0.05)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.noise_threshold < args.benchmark_floor:
        sys.exit(f"[error] --noise-threshold {args.noise_threshold} is below the benchmark's "
                 f"reporting floor {args.benchmark_floor}: variants between them are absent "
                 f"from the published table for unknown reasons, so FDR is not computable "
                 f"there. Raise --noise-threshold, or lower --benchmark-floor if the "
                 f"published table really does report that low.")

    d = load_pips(args.pip_vectors)
    proxy_map = load_proxy_map(args.proxy_map) if args.proxy_map else None
    if proxy_map is None:
        print("[info] EXACT-MATCH scoring (no --proxy-map): a variant in near-perfect LD "
              "with the true one counts as both a miss and a false discovery")

    if args.self_test:
        if args.self_test not in set(d["method"]):
            sys.exit(f"[error] --self-test arm {args.self_test!r} not in {sorted(set(d['method']))}")
        bench = d[d["method"] == args.self_test]
        comparators = [m for m in d["method"].unique() if m != args.self_test]
        bench_label = f"SELF-TEST: {args.self_test} (stand-in)"
        print(f"[self-test] benchmark = {args.self_test} (stand-in, NOT the published result)")
    elif args.benchmark:
        bench = load_pips(args.benchmark)
        avail = sorted(set(bench["method"]))
        if args.benchmark_method:
            if args.benchmark_method not in avail:
                sys.exit(f"[error] --benchmark-method {args.benchmark_method!r} not in {avail}")
            bench = bench[bench["method"] == args.benchmark_method]
        elif len(avail) > 1:
            sys.exit(f"[error] {args.benchmark} holds several benchmark sources {avail}. "
                     f"Pooling them would double-count variants and mix two different "
                     f"ground truths -- pass --benchmark-method to pick one.")
        bench_label = str(bench["method"].iloc[0])
        print(f"[info] benchmark = {bench_label} "
              f"({len(bench)} variants, {bench['locus'].nunique()} loci)")
        comparators = list(d["method"].unique())
    else:
        sys.exit("[error] give --benchmark, or --self-test ARM to exercise the code")

    if args.loci:
        keep = {l.strip() for l in open(args.loci) if l.strip()}
        d = d[d["locus"].isin(keep)]
        bench = bench[bench["locus"].isin(keep)]
        print(f"[info] restricted to {len(keep)} loci")

    # A locus the benchmark does not cover has no ground truth, only silence. Scoring
    # comparator hits there would count EVERY one as a false discovery by construction
    # -- the published SuSiEx table lists only loci that have a credible set at all, so
    # 114 of our 250 target loci are simply absent from it. Restrict to loci the
    # benchmark actually speaks to.
    if args.top_variant_only:
        n_before = (len(d), len(bench))
        d, bench = top_variant_only(d), top_variant_only(bench)
        print(f"[info] TOP-VARIANT-ONLY: comparator {n_before[0]:,} -> {len(d):,} rows, "
              f"benchmark {n_before[1]:,} -> {len(bench):,} (one per method x locus)")

    if not args.no_restrict_to_benchmark_loci:
        bench_loci = set(bench["locus"].unique())
        dropped = sorted(set(d["locus"].unique()) - bench_loci)
        d = d[d["locus"].isin(bench_loci)]
        print(f"[info] scoring on {len(bench_loci)} loci the benchmark covers; "
              f"dropped {len(dropped)} loci absent from it (no ground truth there)")

    # The sweep cannot start below the benchmark's own reporting floor. Sensitivity
    # at t asks "of variants the benchmark puts above t, how many did we keep?" --
    # below the floor we cannot see the benchmark's variants at all, so every
    # threshold under it silently reuses the floor's variant set and manufactures
    # resolution the published table does not have.
    lo = args.min_threshold if args.min_threshold is not None else args.benchmark_floor
    if lo < args.benchmark_floor:
        sys.exit(f"[error] --min-threshold {lo} is below the benchmark's reporting floor "
                 f"{args.benchmark_floor}; sensitivity is undefined there.")
    thresholds = np.round(np.arange(lo, 0.95 + 1e-9, args.step), 4)
    print(f"[info] sweeping {len(thresholds)} thresholds, {lo:g} to {thresholds[-1]:g}, "
          f"step {args.step:g}")
    out = []
    for m in comparators:
        cv = curve(bench, d[d["method"] == m], thresholds, args.noise_threshold,
                   bench_threshold=args.benchmark_threshold, proxy_map=proxy_map,
                   restrict_to_shared=not args.count_unfitted_as_false)
        cv.insert(0, "method", m)
        base = m.split("_")[0]
        cv.insert(1, "arm", "per-ancestry" if base in PER_ANCESTRY
                  else "merged-LD" if base in MERGED else "other")
        # provenance travels with the numbers so a figure can never misstate its
        # own ground truth (this file was briefly generated against a pooled
        # SuSiEx+PGC-FINEMAP benchmark before that was caught)
        cv["benchmark"] = bench_label
        cv["benchmark_loci"] = bench["locus"].nunique()
        cv["benchmark_threshold"] = args.benchmark_threshold
        cv["noise_threshold"] = args.noise_threshold
        cv["fdr_domain"] = ("comparator calls at benchmark-fitted variants only"
                            if not args.count_unfitted_as_false else "all comparator calls")
        cv["matching"] = ("LD-aware (r2>0.95)" if proxy_map else "exact variant match") \
            + (", top variant only" if args.top_variant_only else "")
        out.append(cv)
    res = pd.concat(out, ignore_index=True)
    res.to_csv(args.out, index=False)
    print(f"[done] wrote {args.out} ({len(res)} rows, {len(comparators)} methods, "
          f"{len(thresholds)} thresholds, noise<{args.noise_threshold})")

    show = res[np.isclose(res["threshold"], 0.5)]
    if len(show):
        print("\nAt t=0.5:")
        print(show[["method", "arm", "n_benchmark_positive", "sensitivity",
                    "n_comparator_positive", "fdr"]].to_string(index=False))


if __name__ == "__main__":
    main()
