#!/usr/bin/env python3
"""B71: recompute B53's per-panel s (Supp Fig 10 a, b) with the real susieR estimate_s_rss (susier_s.py).

Same variant set as eur_panel_compare_b53.py build_one: EUR sumstats, non-palindromic, present in all three EUR panels
(GCTB-10k, AoU EUR, PolyFun UKB), in master order; 249 loci (chr3_48182866_53553354 has no PolyFun). Read-only:
writes no z/ld files and runs no fits.

OUT/s_stats_susier.csv columns: locus, panel, n (variants), n_gwas (median per-variant GWAS N, as B70),
s_b53 (ld_diagnostic.estimate_s_scaled_legacy, then named estimate_s_and_kriging, what B53 reported), s_susier_n (primary), s_susier (no n), min_eig.
  V2_LD_ROOT=... V2_SS_ROOT=... python s_stat_eur_panels_b71.py --loci-file v2_inputs/loci_250.txt --out task_b71_out
"""
import argparse
import concurrent.futures as cf
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import v2_io  # noqa: E402
from harmonize_v2 import ordered  # noqa: E402
from ld_diagnostic import estimate_s_scaled_legacy  # noqa: E402
from susier_s import estimate_s_rss  # noqa: E402

PANELS = ("eur_gctb", "eur_aou", "eur_polyfun")
NO_POLYFUN = "chr3_48182866_53553354"


def one(L):
    m = v2_io.load_master(L)
    ss = v2_io.load_ss("eur", L)
    ld = {p: v2_io.load_ld(p, L) for p in PANELS}
    base = set(ss.index) & set(m.index[~m["palindromic"]])
    common = ordered(m, base.intersection(*[set(ld[p][1].index) for p in PANELS]))
    z = ss.loc[common, "z"].to_numpy(float)
    n_gwas = float(np.nanmedian(ss.loc[common, "n"])) if "n" in ss.columns and len(common) else np.nan
    rows = []
    for p in PANELS:
        r = dict(locus=L, panel=p, n=len(common), n_gwas=n_gwas, s_b53=np.nan, s_susier_n=np.nan, s_susier=np.nan,
                 min_eig=np.nan)
        if len(common) > 1:
            R = v2_io.subset(*ld[p], common)[0].astype(np.float64)
            r.update(s_b53=float(estimate_s_scaled_legacy(z, R)[0]),
                     s_susier_n=estimate_s_rss(z, R, n_gwas) if np.isfinite(n_gwas) else np.nan,
                     s_susier=estimate_s_rss(z, R), min_eig=float(np.linalg.eigvalsh(R).min()))
        rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci-file", default=os.path.join(HERE, "v2_inputs", "loci_250.txt"))
    ap.add_argument("--out", default=os.path.join(HERE, "task_b71_out"))
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    loci = [l.strip() for l in open(a.loci_file) if l.strip() and l.strip() != NO_POLYFUN]
    rows = []
    with cf.ProcessPoolExecutor(a.jobs) as ex:
        for i, r in enumerate(ex.map(one, loci), 1):
            rows += r
            if i % 25 == 0:
                print(f"{i} loci", flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(a.out, "s_stats_susier.csv"), index=False)
    cols = ["s_b53", "s_susier_n", "s_susier"]
    print(t.groupby("panel")[cols].describe(percentiles=[0.25, 0.5, 0.75]).T.to_string())
    w = t.pivot(index="locus", columns="panel", values="s_susier_n")
    for p in ("eur_aou", "eur_polyfun"):
        print(f"s_susier_n: {p} > eur_gctb at {(w[p] > w['eur_gctb']).sum()}/{len(w)} loci")
    print("Spearman s_b53 vs s_susier_n per panel:",
          {p: round(g.s_b53.corr(g.s_susier_n, method="spearman"), 4) for p, g in t.groupby("panel")})


if __name__ == "__main__":
    main()
