#!/usr/bin/env python3
"""Extra: SuSiEx with AoU EUR + AoU EAS LD vs the primary v2 SuSiEx (GCTB EUR + AoU EAS), noambig.
Same z, same argv (run_finemap.py's _susiex), same EAS inputs; only the EUR LD panel (and so the EUR domain) differs.
Both runs are read with export_readers.susiex (max PIP over CSs; CSs gated on SuSiEx's own purity >= MIN_PURITY).
Usage: python susiex_eur_panel.py GCTB_WORK GCTB_TREE AOU_WORK AOU_TREE OUT_DIR"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipeline"))
import export_readers as E  # noqa: E402
import run_finemap as D  # noqa: E402  (E.out_dir reads D.WORK, so D.WORK must be set per run)

def read(work, tree):
    D.WORK, D.L_OVERRIDE, E.UNION = work, None, os.path.join(tree, "union")
    assert E.out_dir("SuSiEx") == os.path.join(work, "SuSiEx"), E.out_dir("SuSiEx")
    loci = [l.strip() for l in open(os.path.join(HERE, "..", "inputs", "loci_250.txt")) if l.strip()]
    pip, cs = [], []
    for l in loci:
        rows, c, native, _, st = E.susiex(l, {})
        pip += [(l, r[0], float(r[4])) for r in rows]
        cs += [(l, k, len(v), native[k], ";".join(v)) for k, v in c.items()]
    return (pd.DataFrame(pip, columns=["locus", "rsid", "pip"]),
            pd.DataFrame(cs, columns=["locus", "cs_id", "size", "purity", "members"]))

gw, gt, aw, at, out = [os.path.expanduser(x) for x in sys.argv[1:6]]
os.makedirs(out, exist_ok=True)
(pg, cg), (pa, ca) = read(gw, gt), read(aw, at)
pg.to_csv(f"{out}/susiex_gctb_noambig_pip.csv.gz", index=False); pa.to_csv(f"{out}/susiex_euraou_noambig_pip.csv.gz", index=False)
cg.to_csv(f"{out}/susiex_gctb_noambig_cs.csv.gz", index=False); ca.to_csv(f"{out}/susiex_euraou_noambig_cs.csv.gz", index=False)
incs = lambda c: {(l, m) for l, ms in zip(c.locus, c.members) for m in ms.split(";")}
sg, sa = incs(cg), incs(ca)
summ = []
for name, p, c in (("GCTB EUR", pg, cg), ("AoU EUR", pa, ca)):
    summ.append(dict(panel=name, variants=len(p), gated_cs=len(c), loci_with_cs=c.locus.nunique(),
                     pip_gt_05=int((p.pip > .5).sum()), pip_gt_095=int((p.pip > .95).sum()),
                     median_cs_size=float(c["size"].median()) if len(c) else np.nan))
S = pd.DataFrame(summ); print(S.to_string(index=False))
j = pg.merge(pa, on=["locus", "rsid"], suffixes=("_gctb", "_aou"))
j["d"] = (j.pip_gctb - j.pip_aou).abs()
print(f"\nshared variants {len(j):,} (GCTB-only {len(pg) - len(j):,}, AoU-only {len(pa) - len(j):,})")
print(f"PIP Pearson r {np.corrcoef(j.pip_gctb, j.pip_aou)[0, 1]:.4f}; max |dPIP| {j.d.max():.3f}; "
      f"loci with max |dPIP| > 0.1: {j.loc[j.d > .1, 'locus'].nunique()}; > 0.5: {j.loc[j.d > .5, 'locus'].nunique()}")
print(f"0.5 crossings {int(((j.pip_gctb > .5) != (j.pip_aou > .5)).sum())}; 0.95 crossings {int(((j.pip_gctb > .95) != (j.pip_aou > .95)).sum())}")
tg = pg.loc[pg.groupby("locus").pip.idxmax()].set_index("locus"); ta = pa.loc[pa.groupby("locus").pip.idxmax()].set_index("locus")
top = tg.join(ta, lsuffix="_gctb", rsuffix="_aou"); top["same"] = top.rsid_gctb == top.rsid_aou
print(f"top variant same at {int(top.same.sum())}/{len(top)} loci")
nb = pd.DataFrame({"gctb": cg.groupby("locus").size(), "aou": ca.groupby("locus").size()}).fillna(0).astype(int)
nb = nb.reindex(sorted(set(pg.locus)), fill_value=0)
print(f"gated CS count per locus: same {int((nb.gctb == nb.aou).sum())}, AoU more {int((nb.aou > nb.gctb).sum())}, GCTB more {int((nb.gctb > nb.aou).sum())}")
print(f"CS membership (variant in a gated CS): both {len(sg & sa):,}, GCTB only {len(sg - sa):,}, AoU only {len(sa - sg):,}")
big = j.loc[j.groupby("locus").d.idxmax()].sort_values("d", ascending=False)
big = big[big.d > .3]
print("\nloci with max |dPIP| > 0.3 (variant, PIP GCTB -> AoU; gated CS counts):")
for r in big.itertuples():
    print(f"  {r.locus} {r.rsid}: {r.pip_gctb:.3f} -> {r.pip_aou:.3f}  (CS {nb.at[r.locus, 'gctb']} -> {nb.at[r.locus, 'aou']})")
top.to_csv(f"{out}/top_variant.csv"); nb.to_csv(f"{out}/gated_cs_per_locus.csv"); S.to_csv(f"{out}/summary.csv", index=False)
big.to_csv(f"{out}/loci_dpip_gt_0.3.csv", index=False)
# benchmark top-13 (Fig 2D set)
MET = os.path.join(HERE, "..", "figures")
k193 = set(pd.read_csv(os.path.join(MET, "figures", "fig1c_193_C_per_locus.csv")).locus)
b = pd.read_csv(os.path.join(MET, "pip_vectors_susiex_published_full.csv.gz"))
b = b[b.locus.isin(k193) & (b.pip > .9)]
print("\nbenchmark top-13: lead PIP GCTB / AoU, top variant GCTB / AoU")
for r in b.sort_values("locus").itertuples():
    g = pg[(pg.locus == r.locus) & (pg.rsid == r.rsid)].pip; a = pa[(pa.locus == r.locus) & (pa.rsid == r.rsid)].pip
    print(f"  {r.locus} {r.rsid}: {g.iloc[0] if len(g) else float('nan'):.3f} / {a.iloc[0] if len(a) else float('nan'):.3f}; "
          f"top {tg.at[r.locus, 'rsid']} / {ta.at[r.locus, 'rsid']}")
