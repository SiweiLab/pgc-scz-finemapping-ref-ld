#!/usr/bin/env python3
"""Gate 7 of pipeline v2: did the fits run cleanly? Laptop-side, from the committed exports only.

Per arm x variant (task_v2_out/{noambig,primary}/) and per sweep point (task_v2_out/lsweep/<V>/):
  loci with output (expect 250), rows = the domain B52/B54 logged, rsids not in the master (expect 0),
  status: converged = FALSE, n_cs, loci with a CS; CSs without a common_purity; PIP > 0.5 / 0.95 counts;
  beside v1's export (b59d67a: task_b30_out primary, task_b34_out noambig) for orientation.
Seeded vs unseeded CARMA: the interim exports (CARMA = the unseeded first pass) against the final seeded export.
Writes task_v2_out/gate7_tables.md (+ .csv per table); REPORT.md is written from these.
"""
import glob
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
T = os.path.join(HERE, "task_v2_out")
ARMS = ["SuSiEx-ours", "MultiSuSiE", "MESuSiE", "SuSiE-merged", "FINEMAP-merged", "RSparsePro", "CARMA",
        "EUR-SuSiE-RSS", "EUR-FINEMAP"]
DOMAIN = {  # B52 gate 6 (+ B54): rows each arm must have
    "noambig": {"SuSiEx-ours": 228397, "MultiSuSiE": 196077, "MESuSiE": 133647, "SuSiE-merged": 188125,
                "FINEMAP-merged": 188125, "RSparsePro": 188125, "CARMA": 188125, "EUR-SuSiE-RSS": 194135,
                "EUR-FINEMAP": 194135},
    "primary": {"SuSiEx-ours": 269466, "MultiSuSiE": 231220, "MESuSiE": 157998, "SuSiE-merged": 221830,
                "FINEMAP-merged": 221830, "RSparsePro": 221830, "CARMA": 221830, "EUR-SuSiE-RSS": 228853,
                "EUR-FINEMAP": 228853}}
V1 = {"noambig": ("task_b34_out", "-noambig"), "primary": ("task_b30_out", "")}
M = pd.read_csv(os.path.join(HERE, "v2_inputs", "master_variants_250.tsv.gz"), sep="\t", dtype=str, keep_default_na=False)
MS = set(zip(M.locus, M.rsid))


def summarise(pipf, csf):
    p = pd.read_csv(pipf, dtype={"rsid": str}, low_memory=False)
    c = pd.read_csv(csf, dtype={"rsid": str}, low_memory=False)
    inc = c[c.cs_ids_all.notna()]
    if len(inc):
        e = inc.assign(cs=inc.cs_ids_all.astype(str).str.split(";")).explode("cs")
        ncs = (e.locus + "|" + e.cs.astype(str)).nunique()
    else:
        ncs = 0
    return dict(loci=p.locus.nunique(), rows=len(p), not_in_master=int(sum((l, r) not in MS for l, r in zip(p.locus, p.rsid))),
                n_cs=int(ncs), loci_with_cs=int(inc.locus.nunique()),
                cs_no_common_purity=int(inc.drop_duplicates(["locus", "cs_ids_all"]).common_purity.isna().sum()),
                pip05=int((p.pip > 0.5).sum()), pip095=int((p.pip > 0.95).sum()), nan_pip=int(p.pip.isna().sum()))


def main():
    out, rows = [], []
    for V, suf in (("noambig", "-noambig"), ("primary", "")):
        st = pd.read_csv(os.path.join(T, V, "primary_status.csv.gz"))
        for a in ARMS:
            lab = a + suf
            s = summarise(os.path.join(T, V, f"primary_{lab}_pip.csv.gz"), os.path.join(T, V, f"primary_{lab}_cs_purity.csv.gz"))
            ss = st[st.method == lab]
            conv = ss.converged.astype(str).str.upper()
            r = dict(variant=V, arm=a, **s, domain=DOMAIN[V][a], rows_eq_domain=s["rows"] == DOMAIN[V][a],
                     status_rows=len(ss), converged_false=int((conv == "FALSE").sum()))
            vd, vs = V1[V]
            vf = os.path.join(HERE, vd, f"primary_{a}{vs}_pip.csv.gz")
            if a == "SuSiEx-ours" and V == "noambig":
                vf = os.path.join(HERE, "task_b33_out", "primary_SuSiEx-noambig_pip.csv.gz")
            if os.path.exists(vf):
                v1 = pd.read_csv(vf, usecols=["pip"])
                r.update(v1_rows=len(v1), v1_pip05=int((v1.pip > 0.5).sum()), v1_pip095=int((v1.pip > 0.95).sum()))
            rows.append(r)
    main_t = pd.DataFrame(rows)
    main_t.to_csv(os.path.join(T, "gate7_main.csv"), index=False)
    out.append("## Main run, 9 arms x 2 variants\n\n" + main_t.to_markdown(index=False))

    # seeded vs unseeded CARMA
    cr = []
    for V, suf in (("noambig", "-noambig"), ("primary", "")):
        s = pd.read_csv(os.path.join(T, V, f"primary_CARMA{suf}_cs_purity.csv.gz"), dtype={"rsid": str}, low_memory=False)
        u = pd.read_csv(os.path.join(T, f"{V}_interim", f"primary_CARMA{suf}_cs_purity.csv.gz"), dtype={"rsid": str}, low_memory=False)
        loci = sorted(set(u.locus) & set(s.locus))
        j = s[s.locus.isin(loci)].merge(u[u.locus.isin(loci)], on=["locus", "rsid"], suffixes=("_s", "_u"))
        d = (j.pip_s - j.pip_u).abs()
        top = lambda x: x.loc[x.groupby("locus").pip.idxmax(), ["locus", "rsid"]].set_index("locus").rsid
        ts, tu = top(s[s.locus.isin(loci)]), top(u[u.locus.isin(loci)])
        ncs = lambda x: x[x.cs_ids_all.notna()].drop_duplicates(["locus", "cs_ids_all"]).groupby("locus").size().reindex(loci, fill_value=0)
        cr.append(dict(variant=V, loci_compared=len(loci), max_abs_dpip=round(float(d.max()), 4),
                       loci_max_dpip_gt_0_05=int(j.loc[d > 0.05, "locus"].nunique()), loci_max_dpip_gt_0_2=int(j.loc[d > 0.2, "locus"].nunique()),
                       top_variant_changed=int((ts != tu.reindex(ts.index)).sum()), loci_cs_count_changed=int((ncs(s) != ncs(u)).sum()),
                       pip05_seeded=int((s.pip > 0.5).sum()), pip05_unseeded=int((u.pip > 0.5).sum())))
    ct = pd.DataFrame(cr)
    ct.to_csv(os.path.join(T, "gate7_carma_seeded_vs_unseeded.csv"), index=False)
    out.append("## CARMA, seeded (final) vs unseeded first pass (interim)\n\n" + ct.to_markdown(index=False))

    # L sweep
    sw = []
    for V, suf in (("noambig", "-noambig"), ("primary", "")):
        for f in sorted(glob.glob(os.path.join(T, "lsweep", V, "primary_*_pip.csv.gz"))):
            lab = os.path.basename(f)[len("primary_"):-len("_pip.csv.gz")]
            arm, L = lab.rsplit("_L", 1)
            base = arm.replace(suf, "") if suf else arm
            s = summarise(f, f.replace("_pip.csv.gz", "_cs_purity.csv.gz"))
            dom = DOMAIN[V]["SuSiEx-ours" if base == "SuSiEx-nomultistep" else base]
            sw.append(dict(variant=V, arm=base, L=int(L), **s, rows_eq_domain=s["rows"] == dom))
    swt = pd.DataFrame(sw).sort_values(["variant", "arm", "L"])
    swt.to_csv(os.path.join(T, "gate7_lsweep.csv"), index=False)
    out.append("## L sweep (primary L reused, not refitted)\n\n" + swt.to_markdown(index=False))
    open(os.path.join(T, "gate7_tables.md"), "w").write("\n\n".join(out) + "\n")
    print("\n\n".join(out))


if __name__ == "__main__":
    main()
