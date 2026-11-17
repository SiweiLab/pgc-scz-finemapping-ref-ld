#!/usr/bin/env python3
"""Per-ancestry z/LD s with the real susieR estimate_s_rss (susier_s.py, verified against R).

Per locus and ancestry, non-palindromic master variants with
that ancestry's sumstats AND its reference LD. EUR = EUR z vs GCTB-10k, EAS = EAS z vs AoU EAS. Read-only; no fits.

Columns written to OUT/s_susier_per_ancestry.csv:
  n_snp       variants in the set
  n_gwas      median per-variant GWAS N from the sumstats (Weight / NEFF), NaN if the file carries neither
  max_abs_z   largest |z| in the set
  s_legacy       ld_diagnostic.estimate_s_scaled_legacy, then named estimate_s_and_kriging (the estimate first reported; free scale, no n adjustment)
  s_susier_n  susieR estimate_s_rss(z, R, n = n_gwas)     <- the primary value
  s_susier    susieR estimate_s_rss(z, R) without n
  LD_ROOT=... SS_ROOT=... python s_per_ancestry.py --loci-file ../inputs/loci_250.txt --out ../results/s_statistic
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
import ld_io  # noqa: E402
from ld_diagnostic import estimate_s_scaled_legacy  # noqa: E402
from susier_s import estimate_s_rss  # noqa: E402

PANELS = (("eur", "eur_gctb"), ("eas", "eas_aou"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci-file", default=os.path.join(HERE, "..", "inputs", "loci_250.txt"))
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results", "s_statistic"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for i, L in enumerate([l.strip() for l in open(a.loci_file) if l.strip()], 1):
        m = ld_io.load_master(L)
        keep = set(m.index[~m["palindromic"]])
        for anc, src in PANELS:
            ss = ld_io.load_ss(anc, L)
            R, snps = ld_io.load_ld(src, L)
            common = [v for v in m.index if v in keep and v in ss.index and v in snps.index and np.isfinite(ss.at[v, "z"])]
            row = dict(locus=L, ancestry=anc.upper(), panel=src, n_snp=len(common), n_gwas=np.nan, max_abs_z=np.nan,
                       s_legacy=np.nan, s_susier_n=np.nan, s_susier=np.nan)
            if len(common) >= 2:
                Rc = ld_io.subset(R, snps, common)[0].astype(np.float64)
                z = ss.loc[common, "z"].to_numpy(float)
                n_gwas = float(np.nanmedian(ss.loc[common, "n"])) if "n" in ss.columns else np.nan
                row.update(n_gwas=n_gwas, max_abs_z=float(np.max(np.abs(z))),
                           s_legacy=float(estimate_s_scaled_legacy(z, Rc)[0]),
                           s_susier_n=estimate_s_rss(z, Rc, n_gwas) if np.isfinite(n_gwas) else np.nan,
                           s_susier=estimate_s_rss(z, Rc))
            rows.append(row)
        if i % 25 == 0:
            print(f"{i} loci", flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(a.out, "s_susier_per_ancestry.csv"), index=False)
    cols = ["s_legacy", "s_susier_n", "s_susier"]
    print(t.groupby("ancestry")[cols].describe(percentiles=[0.25, 0.5, 0.75]).T.to_string())
    print("Spearman with s_legacy:\n", t.groupby("ancestry")[cols].corr(method="spearman").to_string())
    print("Spearman with max_abs_z:\n", t.groupby("ancestry")[cols + ["max_abs_z"]].corr(method="spearman")["max_abs_z"].to_string())
    for c in cols:
        mx = t.groupby("locus")[c].max()
        print(f"{c}: loci with max-over-ancestries > 0.005 / 0.01 / 0.02:",
              [(mx > x).sum() for x in (0.005, 0.01, 0.02)], "| > 0.01:", sorted(mx[mx > 0.01].index))


if __name__ == "__main__":
    main()
