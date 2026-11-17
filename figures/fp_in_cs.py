#!/usr/bin/env python3
"""False positives inside credible sets, per arm, v2 (same join and label as Fig. 1B in-CS).

FP = a CS member of the arm whose benchmark (Yuan et al.) PIP <= 0.1; TP = PIP > 0.1.
Loci: the 193 where the benchmark formed a CS; variants: the arm x benchmark inner join on
variant_key, so a member the benchmark never tested is counted separately (not_in_bench),
never as an FP. Reported for all members and for members with our PIP > 0.1 / 0.5 / 0.9.
  python3 fp_in_cs.py [--export ../results/finemap/main --suffix=-noambig]
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_sensfdr import bench_cs_loci, BENCH, BENCH_RAW
from sensitivity_fdr import load_pips, variant_key
from plot_cs_split import METHOD_ORDER, disp

HERE = os.path.dirname(os.path.abspath(__file__))
from paths import V2
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--export", default=os.path.join(HERE, "..", "results", "finemap", "main"))
ap.add_argument("--suffix", default="-noambig")
ap.add_argument("--label-pip", type=float, default=0.1)
ap.add_argument("--out", default=os.path.join(V2, "fp_in_cs.csv"))
a = ap.parse_args()

full, _ = bench_cs_loci(BENCH_RAW)
b = load_pips(BENCH)
b = b[b.locus.isin(full)].groupby(["locus", "key"], as_index=False).pip.max().rename(columns={"pip": "pb"})
rows = []
for arm in METHOD_ORDER:
    f = os.path.join(a.export, f"primary_{arm}{a.suffix}_cs_purity.csv.gz")
    if not os.path.exists(f):
        continue
    d = pd.read_csv(f, low_memory=False)
    d = d[d.locus.isin(full) & d.cs_ids_all.notna()].copy()
    d["key"] = variant_key(d)
    u = d.merge(b, on=["locus", "key"], how="left")
    inb = u.pb.notna()
    r = dict(arm=disp(arm), loci_with_cs=d.locus.nunique(), cs_members=len(u), not_in_bench=int((~inb).sum()))
    for s, tag in ((None, "all"), (0.1, "pip>0.1"), (0.5, "pip>0.5"), (0.9, "pip>0.9")):
        k = inb & (True if s is None else u.pip > s)
        fp = k & (u.pb <= a.label_pip); tp = k & (u.pb > a.label_pip)
        r[f"FP_{tag}"], r[f"TP_{tag}"] = int(fp.sum()), int(tp.sum())
        r[f"FDR_{tag}"] = round(fp.sum() / k.sum(), 3) if k.sum() else float("nan")
        if tag == "pip>0.5":
            r["loci_with_FP_pip>0.5"] = u.loc[fp, "locus"].nunique()
    rows.append(r)
t = pd.DataFrame(rows)
t.to_csv(a.out, index=False)
print(t.to_string(index=False))
