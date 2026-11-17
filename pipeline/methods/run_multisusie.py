#!/usr/bin/env python3
"""
Run MultiSuSiE (multisusie_rss) on a harmonized locus.

Reads {locus}_harmonized.npz produced by ../harmonize.py and calls the
summary-statistics entry point. Effect sizes + SEs are passed directly (b_list /
s_list); pass --use-z to use z_list instead.

Install: see https://github.com/jordanero/MultiSuSiE
  conda env create -f environment.yml && conda activate MultiSuSiE && pip install .

Outputs (in --out):
  {locus}_multisusie_pip.tsv   rsID, chr, pos, alleles, pip
  {locus}_multisusie_cs.tsv    credible-set membership (cs_id, rsID, pip)

Usage:
  python run_multisusie.py chr1_5495783_8495783 \\
      --harmonized-dir /scratch/.../ma_harmonized \\
      --out            /scratch/.../multisusie_out \\
      --L 10
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd
import MultiSuSiE


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("locus")
    p.add_argument("--harmonized-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--L", type=int, required=True, help="Max number of causal effects")
    p.add_argument("--use-z", action="store_true",
                   help="Pass z_list instead of b_list/s_list")
    p.add_argument("--coverage", type=float, required=True)
    p.add_argument("--min-abs-corr", type=float, required=True,
                   help="Purity threshold for credible sets")
    p.add_argument("--rho", type=float, required=True,
                   help="Cross-population effect-size correlation (fixed scalar). "
                        "Library default is 0.75 if omitted -- every published run "
                        "of this wrapper used that default implicitly, since this "
                        "flag didn't exist before 2026-09-08.")
    p.add_argument("--population-sizes", required=True,
                   help="GWAS N per ancestry, keyed by name, e.g. EUR=58749,EAS=14425. "
                        "Required: before this the wrapper used the npz's per-locus "
                        "median-Weight `n` (EUR 57,602-58,749 across loci), which no caller "
                        "could override.")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)

    d = np.load(os.path.join(args.harmonized_dir, f"{args.locus}_harmonized.npz"),
                allow_pickle=True)
    K = len(d["ancestries"])
    R_list = [d["R"][k].astype(np.float64) for k in range(K)]
    # Population sizes come from the caller, matched to the npz's ancestry order BY NAME, so
    # a reordered argument cannot silently swap EUR and EAS.
    ancestries = [str(a) for a in d["ancestries"]]
    n_by_anc = dict(kv.split("=") for kv in args.population_sizes.split(","))
    if sorted(n_by_anc) != sorted(ancestries):
        sys.exit(f"--population-sizes names {sorted(n_by_anc)} do not match the npz "
                 f"ancestries {ancestries}")
    N_list = [int(n_by_anc[a]) for a in ancestries]
    varY_list = [float(x) for x in d["varY"]]

    kwargs = dict(R_list=R_list, varY_list=varY_list, population_sizes=N_list,
                  L=args.L, coverage=args.coverage, min_abs_corr=args.min_abs_corr)
    # rho is required now: the published runs used 0.75 but never passed it, so it silently
    # took whatever multisusie_rss defaults to. Passing it makes the value a recorded fact.
    kwargs["rho"] = args.rho
    if args.use_z:
        kwargs["z_list"] = [d["z"][k].astype(np.float64) for k in range(K)]
    else:
        kwargs["b_list"] = [d["beta"][k].astype(np.float64) for k in range(K)]
        kwargs["s_list"] = [d["se"][k].astype(np.float64) for k in range(K)]

    fit = MultiSuSiE.multisusie_rss(**kwargs)

    pip = np.asarray(fit.pip)
    out = pd.DataFrame({
        "rsid": d["rsids"], "chr": d["chr"], "pos": d["pos"],
        "effect_allele": d["effect_allele"], "other_allele": d["other_allele"],
        "pip": pip,
    })
    out.to_csv(os.path.join(args.out, f"{args.locus}_multisusie_pip.tsv"),
               sep="\t", index=False)

    # fit.sets == [cs, purity, claimed_coverage, include_mask]; cs is a length-L list
    # of variant-index arrays. include_mask[k] is MultiSuSiE's authoritative flag for
    # whether single-effect k forms a purity-passing (min_abs_corr) credible set.
    cs_rows = []
    sets = getattr(fit, "sets", None)
    cs_list = sets[0] if isinstance(sets, (list, tuple)) and len(sets) else []
    purity = sets[1] if isinstance(sets, (list, tuple)) and len(sets) > 1 else None
    cov = sets[2] if isinstance(sets, (list, tuple)) and len(sets) > 2 else None
    mask = sets[3] if isinstance(sets, (list, tuple)) and len(sets) > 3 else None

    def _scalar(arr, k):
        try:
            v = np.ravel(arr[k])
            return float(v[0]) if v.size else np.nan
        except Exception:
            return np.nan

    n_cs = 0
    n_dropped_nan_purity = 0
    for k, members in enumerate(cs_list):
        members = list(members)
        if len(members) == 0:
            continue
        if mask is not None and not bool(np.ravel(mask)[k]):   # purity-pass per MultiSuSiE
            continue
        # ── NaN-purity guard (added 2026-07-30) ──────────────
        # MultiSuSiE's own purity filter cannot reject a credible set whose purity is
        # NaN, so such sets reach us with include_mask=True even though their purity
        # was never actually established. Chain, from MultiSuSiE/susiepy_ss.py:
        #   1. our v2-union harmonization keeps SNPs missing from ONE ancestry's
        #      sumstats; for a few low-frequency EUR-only variants the SNP is also
        #      absent from the EAS LD panel, so its whole EAS LD row/col is NaN;
        #   2. multisusie_rss calls susie_get_cs(..., n_purity=np.inf), overriding
        #      that function's own default of 100, so purity is computed over EVERY
        #      member pair rather than a <=100-member subsample -- the bad SNP is
        #      included with certainty rather than ~3% of the time;
        #   3. get_purity_x reduces with functools.reduce(np.maximum, ...), and
        #      np.maximum PROPAGATES NaN (np.fmax would not), so one missing-ancestry
        #      entry makes the cross-ancestry max NaN even where EUR is complete;
        #      np.min of that is NaN;
        #   4. the filter is `include_mask[purity < min_abs_corr] = False`, and
        #      `NaN < 0.5` is False, so the mask is never cleared and the set is kept.
        # Net: any CS containing one such SNP bypasses the purity filter entirely,
        # which is how four loci emitted 22 near-duplicate ~3,500-SNP "credible sets"
        # (pairwise Jaccard 0.93-0.999, member PIPs summing to ~9 instead of ~1).
        # Treating an unestablished purity as a FAILURE is the conservative reading of
        # what min_abs_corr is for, and it makes the unfiltered and QC-excluded
        # pipelines agree exactly at all four affected loci.
        purity_k = _scalar(purity, k)
        if purity is not None and np.isnan(purity_k):
            n_dropped_nan_purity += 1
            continue
        n_cs += 1
        for j in members:
            cs_rows.append({
                "cs_id": f"L{k+1}", "rsid": d["rsids"][j], "pip": pip[j],
                "min_abs_corr": _scalar(purity, k),
                "claimed_coverage": _scalar(cov, k),
            })
    pd.DataFrame(cs_rows, columns=["cs_id", "rsid", "pip", "min_abs_corr",
                                   "claimed_coverage"]).to_csv(
        os.path.join(args.out, f"{args.locus}_multisusie_cs.tsv"),
        sep="\t", index=False)
    msg = (f"MultiSuSiE {args.locus}: {n_cs} credible sets, "
           f"{len(cs_rows)} CS members, {(pip>0.9).sum()} SNPs PIP>0.9")
    if n_dropped_nan_purity:
        msg += (f" [dropped {n_dropped_nan_purity} CS with NaN purity -- "
                f"purity filter could not evaluate them; see NaN-purity guard above]")
    print(msg)
    # fit.converged is set True on IBSS convergence break, False otherwise (never
    # unset -- see MultiSuSiE.susiepy_ss.susie_multi_ss). No equivalent "lambda"
    # (LD-consistency) diagnostic exists on this fit object -- that quantity is
    # susieR::estimate_s_rss, a property of (z, R, n), already harvested once
    # per locus/ancestry by harvest_s_lambda_b1.py, not of this fit.
    print(f"converged={bool(fit.converged)}")


if __name__ == "__main__":
    main()
