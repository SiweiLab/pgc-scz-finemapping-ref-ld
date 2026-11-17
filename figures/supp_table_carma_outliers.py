#!/usr/bin/env python3
"""Supplementary Table: variants CARMA flagged as inconsistent with the merged LD (seeded v2 no-palindromic run).

Sheet 1, per locus (the loci with >= 1 flagged variant): flagged count, variants fitted, merged-LD s (v2), EUR and
EAS s (susieR estimate_s_rss with n), whether the benchmark lead is flagged, and CARMA's / the ancestry-aware methods' in-CS PIP at the lead.
Sheet 2, the flagged benchmark leads in detail: per method, the in-CS PIP at the lead and whether it is the method's top
credible-set variant. Sheet 3, every flagged variant. Writes output/supp_table_carma_outliers.{xlsx,csv}.
  python3 supp_table_carma_outliers.py
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bench_sensfdr import BENCH, bench_cs_loci, BENCH_RAW  # noqa: E402
from plot_cs_split import disp  # noqa: E402

A = os.path.join(HERE, "..", "results")
EXP = os.path.join(A, "finemap", "main")
ARMS = ["SuSiEx-ours", "MultiSuSiE", "MESuSiE", "FINEMAP-merged", "SuSiE-merged", "RSparsePro", "CARMA"]
PER_ANC = ARMS[:3]
from paths import V2  # noqa: E402
OUT = os.path.join(V2, "supp_table_carma_outliers")


def main():
    o = pd.read_csv(os.path.join(A, "carma_outliers", "carma_outliers_seeded_noambig.tsv"), sep="\t", dtype={"rsid": str})
    flagged = o[o.rsid.notna()]
    s = pd.read_csv(os.path.join(A, "s_statistic", "s_susier_per_ancestry.csv")).pivot(  # real susieR s, with n
        index="locus", columns="ancestry", values="s_susier_n")
    lam = pd.read_csv(os.path.join(A, "lambda_merged", "lambda_merged_noambig.csv"),
                      index_col=0)["lambda"]

    full, _ = bench_cs_loci(BENCH_RAW)
    b = pd.read_csv(BENCH, dtype={"rsid": str})
    b = b[b.locus.isin(full)]
    lead = b.loc[b.groupby("locus").pip.idxmax()].set_index("locus")
    # the lead as Fig. 1c defines it (Yuan's published max-PIP variant); differs from idxmax only at a tie (chr19_19324669)
    pub = pd.read_csv(os.path.join(V2, "fig1c_193_C_per_locus.csv")).drop_duplicates("locus").set_index("locus")
    for L, r in pub.iterrows():
        if L in lead.index and lead.at[L, "rsid"] != r.pub_rsid:
            lead.loc[L, "rsid"] = r.pub_rsid

    incs, top, nvar = {}, {}, {}
    for arm in ARMS:
        d = pd.read_csv(os.path.join(EXP, f"primary_{arm}-noambig_cs_purity.csv.gz"), dtype={"rsid": str}, low_memory=False)
        if arm == "CARMA":
            nvar = d.groupby("locus").size().to_dict()
        c = d[d.cs_ids_all.notna()]
        incs[arm] = {(l, r): p for l, r, p in zip(c.locus, c.rsid, c.pip)}
        top[arm] = c.loc[c.groupby("locus").pip.idxmax()].set_index("locus").rsid.to_dict()

    rows = []
    for L, g in flagged.groupby("locus"):
        ld = lead.rsid.get(L)
        r = dict(locus=L, flagged=len(g), variants_fitted=nvar.get(L), merged_s=lam.get(L),
                 EUR_s=s.EUR.get(L), EAS_s=s.EAS.get(L), benchmark_formed_CS=L in full,
                 benchmark_lead=ld, benchmark_lead_PIP=lead.pip.get(L),
                 lead_flagged=bool(ld is not None and ld in set(g.rsid)),
                 CARMA_inCS_PIP_at_lead=incs["CARMA"].get((L, ld), 0.0) if ld else None,
                 CARMA_top_is_lead=top["CARMA"].get(L) == ld if ld else None,
                 ancestry_aware_top_is_lead=sum(top[a].get(L) == ld for a in PER_ANC) if ld else None)
        rows.append(r)
    t1 = pd.DataFrame(rows).sort_values(["lead_flagged", "flagged"], ascending=[False, False])

    det = []
    for _, r in t1[t1.lead_flagged].iterrows():
        L, ld = r.locus, r.benchmark_lead
        z = flagged[(flagged.locus == L) & (flagged.rsid == ld)].z.iloc[0]
        d = dict(locus=L, benchmark_lead=ld, benchmark_PIP=round(r.benchmark_lead_PIP, 3), CARMA_z=z)
        for a in ARMS:
            d[f"{disp(a)}_inCS_PIP"] = round(incs[a].get((L, ld), 0.0), 3)
            d[f"{disp(a)}_top"] = top[a].get(L) == ld
        det.append(d)
    t2 = pd.DataFrame(det)

    t3 = flagged.merge(lead.rsid.rename("benchmark_lead").reset_index(), on="locus", how="left")
    t3["is_benchmark_lead"] = t3.rsid == t3.benchmark_lead
    t3 = t3.drop(columns="benchmark_lead")

    readme = pd.DataFrame({"note": [
        "CARMA (seed 42, outlier detection on) run on the trans-ancestry meta-analysis z with the merged EUR+EAS reference "
        "LD, no-palindromic variants (pipeline v2). 'flagged' = variants CARMA reported as inconsistent with the LD "
        "(result$Outliers), which it down-weights.",
        "merged_s / EUR_s / EAS_s: susieR estimate_s_rss on the merged matrix with the meta z, and on each ancestry's own "
        "variant set with its own reference LD. Larger = poorer z/LD agreement.",
        "inCS PIP = the method's PIP for the variant if it is in a credible set, else 0. 'top' = it is the method's "
        "highest-PIP credible-set variant at the locus.",
        f"Totals: {len(flagged)} flagged variants at {flagged.locus.nunique()} of {o.locus.nunique()} loci; "
        f"{int(t1.lead_flagged.sum())} flagged variants are the benchmark's lead."]})
    with pd.ExcelWriter(OUT + ".xlsx") as w:
        readme.to_excel(w, sheet_name="README", index=False)
        t1.to_excel(w, sheet_name="per_locus", index=False)
        t2.to_excel(w, sheet_name="flagged_benchmark_leads", index=False)
        t3.to_excel(w, sheet_name="all_flagged_variants", index=False)
    t1.to_csv(OUT + ".csv", index=False)
    pd.set_option("display.width", 250)
    print(readme.note.iloc[-1]); print(t2.to_string(index=False))
    print(f"[done] {OUT}.xlsx")


if __name__ == "__main__":
    main()
