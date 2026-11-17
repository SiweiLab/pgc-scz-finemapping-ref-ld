#!/usr/bin/env python3
"""Supplement: SuSiEx with GCTB EUR + AoU EAS LD (primary) against AoU EUR + AoU EAS, no-palindromic, in-CS.

Reads the Extra run already on disk (results/susiex_eur_panel/): both sides read by
export_readers.susiex, each CS gated on SuSiEx's own purity >= 0.5 against its own EUR panel. No refits.
Against the benchmark (Yuan et al.), on the 193 loci where it formed a CS, joined on (locus, rsid):
  in-CS calls at our PIP > 0.5 / 0.9 split into TP (benchmark PIP > 0.1) and FP (<= 0.1); in-CS AUPRC as Fig. 1B;
  exact recovery as Fig. 1C scenario C (exactly one CS, its top member = the benchmark lead) on the 122/175/193 strata.
Writes output/supp_susiex_eur_panel.csv and prints the FP lists.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bench_sensfdr import bench_cs_loci, BENCH, BENCH_RAW
from paths import V2

X = os.path.join(HERE, "..", "results", "susiex_eur_panel")
full, _ = bench_cs_loci(BENCH_RAW)
b = pd.read_csv(BENCH, dtype={"rsid": str})
b = b[b.locus.isin(full)].groupby(["locus", "rsid"], as_index=False).pip.max().rename(columns={"pip": "pb"})
# the benchmark lead as Fig. 1c defines it: Yuan's published max-PIP variant (Supplementary Table 26), via the per-locus
# recovery table. An idxmax over the PIP vector breaks ties by file order and disagreed at chr19_19324669, where two variants tie
lead = (pd.read_csv(os.path.join(V2, "fig1c_193_C_per_locus.csv")).drop_duplicates("locus")
        .set_index("locus").pub_rsid)
strata = {"193": set(full),
          "175": {l.strip() for l in open(os.path.join(HERE, "..", "inputs", "loci_single175.txt")) if l.strip()} & set(full),
          "122": {l.strip() for l in open(os.path.join(HERE, "..", "inputs", "loci_confident122.txt")) if l.strip()} & set(full)}


def load(tag):
    p = pd.read_csv(f"{X}/susiex_{tag}_noambig_pip.csv.gz", dtype={"rsid": str})
    c = pd.read_csv(f"{X}/susiex_{tag}_noambig_cs.csv.gz")
    mem = c.assign(rsid=c.members.str.split(";")).explode("rsid")[["locus", "cs_id", "rsid"]]
    p = p[p.locus.isin(full)].merge(b, on=["locus", "rsid"], how="left")
    p["incs"] = p.set_index(["locus", "rsid"]).index.isin(mem.set_index(["locus", "rsid"]).index)
    return p, c[c.locus.isin(full)]


def pr_points(u):
    """(recall, precision) over our PIP threshold, in-CS calls, label benchmark PIP > 0.1 (as Fig. 1b); padded with
    precision 0 out to recall 1 beyond the run's reach."""
    pos = (u.pb > 0.1).to_numpy(); nP = pos.sum(); pts = [(0.0, 1.0)]
    for s in np.unique(np.concatenate([np.arange(0.01, 1.0, 0.01), [0.1, 0.5, 0.9]]))[::-1]:
        call = (u.pip.to_numpy() > s) & u.incs.to_numpy(); n = call.sum()
        if n:
            pts.append(((call & pos).sum() / nP, (call & pos).sum() / n))
    r, pr = zip(*sorted(pts)); r, pr = list(r), list(pr)
    if r[-1] < 1:
        r += [r[-1], 1.0]; pr += [0.0, 0.0]
    return r, pr


def auprc(u):
    r, pr = pr_points(u)
    return float(np.trapezoid(pr, r)) if hasattr(np, "trapezoid") else float(np.trapz(pr, r))


RUNS = (("gctb", "GCTB EUR + AoU EAS (primary)"), ("euraou", "AoU EUR + AoU EAS"))


def exact(p, c):
    """Exact recovery per stratum: exactly one CS and its top member is the benchmark lead (Fig. 1c, scenario C)."""
    ncs = c.groupby("locus").size()
    top = p[p.incs].sort_values("pip").groupby("locus").tail(1).set_index("locus").rsid
    return {st: sum(1 for l in L if ncs.get(l, 0) == 1 and top.get(l) == lead.get(l)) for st, L in strata.items()}


def main():
    rows, fps = [], {}
    for tag, name in RUNS:
        p, c = load(tag)
        u = p[p.pb.notna()]
        r = dict(LD=name, variants=len(p), CSs=len(c), loci_with_CS=c.locus.nunique(), AUPRC_inCS=round(auprc(u), 3))
        for s in (0.5, 0.9):
            k = u.incs & (u.pip > s)
            r[f"TP_pip>{s}"], r[f"FP_pip>{s}"] = int((k & (u.pb > 0.1)).sum()), int((k & (u.pb <= 0.1)).sum())
        fps[name] = u[u.incs & (u.pip > 0.5) & (u.pb <= 0.1)][["locus", "rsid", "pip", "pb"]]
        for st, k in exact(p, c).items():
            r[f"exact_{st}"] = f"{k}/{len(strata[st])} ({100 * k / len(strata[st]):.1f}%)"
        rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(V2, "supp_susiex_eur_panel.csv"), index=False)
    pd.set_option("display.width", 250)
    print(t.T.to_string(header=False))
    for k, v in fps.items():
        print(f"\nin-CS FP at PIP > 0.5, {k}:\n" + (v.to_string(index=False) if len(v) else "  none"))


if __name__ == "__main__":
    main()
