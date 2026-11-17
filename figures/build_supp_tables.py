#!/usr/bin/env python3
"""All supplementary tables in one workbook: output/Supplementary_Tables.xlsx (pipeline v2, primary analysis).

Built only from committed outputs; no refits. One sheet per table (some tables span two sheets), plus a Contents sheet.
Numbering follows first citation in the manuscript; reorder `sheets` in main() to renumber. Columns are autofitted.
  python3 build_supp_tables.py
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bench_sensfdr import BENCH, bench_cs_loci, BENCH_RAW  # noqa: E402
from plot_cs_split import disp  # noqa: E402
from sensitivity_fdr import load_pips, variant_key  # noqa: E402

from paths import V2  # noqa: E402
A = os.path.join(HERE, "..", "results")
EXP = os.path.join(A, "finemap", "main")
ARMS = ["SuSiEx-ours", "MultiSuSiE", "MESuSiE", "FINEMAP-merged", "SuSiE-merged", "RSparsePro", "CARMA"]
GROUP = {a: ("ancestry-aware (per-ancestry LD)" if a in ARMS[:3] else "merged LD") for a in ARMS}


def auprc(f):
    return pd.read_csv(os.path.join(V2, f)).groupby("arm").auprc.first()


def t1_summary():
    """Per method: AUPRC (Fig 1b, S4), exact recovery (Fig 1c, S2), resolution (S3)."""
    t = pd.DataFrame(index=ARMS)
    for col, f in (("AUPRC_CS_members", "fig1b_pr_curve_inCS_curve.csv"),
                   ("AUPRC_all_variants", "fig1b_pr_curve_allvar_curve.csv"),
                   ("AUPRC_CS_members_with_AT_CG", "fig1b_pr_curve_inCS_withpal_curve.csv"),
                   ("AUPRC_all_variants_with_AT_CG", "fig1b_pr_curve_allvar_withpal_curve.csv")):
        t[col] = auprc(f).reindex(ARMS).round(3)
    for st in ("122", "175", "193"):
        for sc, lab in (("E", "common_purity"), ("C", "native_purity")):
            r = pd.read_csv(os.path.join(V2, f"fig1c_{st}_{sc}.csv")).set_index("arm").reindex(ARMS)
            t[f"exact_recovery_{st}_{lab}"] = (100 * r.exact_rate).round(1)
            if sc == "E":
                t[f"n_loci_one_CS_{st}"] = r.n_cs_1
    res = pd.read_csv(os.path.join(V2, "supp_resolution_pip05.csv")).set_index("arm").reindex(ARMS)
    t["resolution41_same_variant"] = res.same
    t["resolution41_different_variant"] = res.different
    t["resolution41_not_resolved"] = res.not_confident
    t.insert(0, "group", [GROUP[a] for a in ARMS])
    t.index = [disp(a) for a in ARMS]
    t.index.name = "method"
    g = pd.read_csv(os.path.join(V2, "supp_fig2_recovery_strata_group_means.csv"))
    g["purity"] = g.scenario.map({"E": "common purity >= 0.5", "C": "each method's native filter"})
    g["group_mean_%"] = (100 * g.group_mean).round(1)
    g = g[["stratum", "purity", "group", "group_mean_%", "wilcoxon_P"]]
    return t.reset_index(), g


def t2_high_conf():
    """In-CS calls at PIP > 0.5 / 0.9 by benchmark support, and every unsupported call at PIP > 0.5."""
    s = pd.read_csv(os.path.join(V2, "fp_in_cs.csv"))
    keep = ["arm", "loci_with_cs", "TP_pip>0.5", "FP_pip>0.5", "TP_pip>0.9", "FP_pip>0.9"]
    s = s[s.arm.isin([disp(a) for a in ARMS])][keep].rename(columns={
        "arm": "method", "loci_with_cs": "loci_with_CS", "TP_pip>0.5": "supported_PIP>0.5",
        "FP_pip>0.5": "unsupported_PIP>0.5", "TP_pip>0.9": "supported_PIP>0.9", "FP_pip>0.9": "unsupported_PIP>0.9"})
    full, _ = bench_cs_loci(BENCH_RAW)
    b = load_pips(BENCH)
    b = b[b.locus.isin(full)].groupby(["locus", "key"], as_index=False).pip.max().rename(columns={"pip": "benchmark_PIP"})
    rows = []
    for a in ARMS:
        d = pd.read_csv(os.path.join(EXP, f"primary_{a}-noambig_cs_purity.csv.gz"), low_memory=False)
        d = d[d.locus.isin(full) & d.cs_ids_all.notna() & (d.pip > 0.5)].copy()
        d["key"] = variant_key(d)
        u = d.merge(b, on=["locus", "key"], how="left")
        u = u[u.benchmark_PIP.fillna(0) <= 0.1]
        for r in u.itertuples():
            rows.append(dict(method=disp(a), locus=r.locus, variant=r.rsid, bp=int(r.bp), PIP=round(r.pip, 3),
                             benchmark_PIP=round(r.benchmark_PIP, 4) if pd.notna(r.benchmark_PIP) else None))
    return s, pd.DataFrame(rows).sort_values(["method", "PIP"], ascending=[True, False])


def t3_carma():
    x = os.path.join(V2, "supp_table_carma_outliers.xlsx")
    return {k: pd.read_excel(x, sheet_name=k) for k in ("per_locus", "flagged_benchmark_leads", "all_flagged_variants")}


def t4_s():
    # real susieR estimate_s_rss with n (the earlier ld_diagnostic estimator was not susieR's)
    s = pd.read_csv(os.path.join(A, "s_statistic", "s_susier_per_ancestry.csv")).rename(columns={"n_snp": "n"})
    w = s.pivot(index="locus", columns="ancestry", values="s_susier_n").rename(columns={"EUR": "s_EUR_UKB", "EAS": "s_EAS_AoU"})
    n = s.pivot(index="locus", columns="ancestry", values="n").rename(columns={"EUR": "n_variants_EUR", "EAS": "n_variants_EAS"})
    m = pd.read_csv(os.path.join(A, "lambda_merged", "lambda_merged_noambig.csv"), index_col=0)["lambda"]
    t = w.join(n).join(m.rename("s_merged"))
    return t[["n_variants_EUR", "s_EUR_UKB", "n_variants_EAS", "s_EAS_AoU", "s_merged"]].reset_index()


def autofit(book, skip_rows):
    """Excel's column autofit (previously done by hand): each column's width from its longest displayed value.
    skip_rows: per sheet, 1-based rows left out of the measure (the one-cell table titles and README notes overflow instead)."""
    from openpyxl.utils import get_column_letter
    for ws in book.worksheets:
        skip = skip_rows.get(ws.title, set())
        width = {}
        for row in ws.iter_rows():
            for c in row:
                if c.value is None or c.row in skip:
                    continue
                v = c.value
                t = f"{v:.10g}"[:11] if isinstance(v, float) else str(v)   # Excel's General shows <= 11 characters
                n = sum(1.0 if ch.isdigit() else 0.9 for ch in t) * (1.08 if c.font.b else 1.0)  # digits = 1 width unit
                width[c.column] = max(width.get(c.column, 0), n)
        for col, n in width.items():
            ws.column_dimensions[get_column_letter(col)].width = min(n + 1.5, 255)


def main():
    t1, t1g = t1_summary()
    t2, t2l = t2_high_conf()
    t3 = t3_carma()
    t4 = t4_s()
    # Table order follows first citation (renumbered 2026-10-14): the s statistic (cited in Results, TMTC1 locus) is ST2,
    # the CARMA tables (Methods s-statistic paragraph) ST3, the high-confidence calls ST4. Figure numbers below are the
    # 2026-10-14 renumbering.
    sheets = [
        ("ST1_concordance", t1, "Supplementary Table 1a | Concordance summary per method: AUPRC (Fig. 1b; Supplementary Fig. 2), exact "
         "top-variant recovery (%) on the 122/175/193 loci under both purity conventions (Fig. 1c; Supplementary Fig. 4), and resolution "
         "at the 41 loci with a benchmark top-variant PIP > 0.5 (Supplementary Fig. 7)."),
        ("ST1_group_tests", t1g, "Supplementary Table 1b | Group mean exact recovery and Wilcoxon signed-rank P (per-locus group scores, "
         "paired by locus) for each locus set and purity convention (Supplementary Fig. 4)."),
        ("ST2_s_statistic", t4, "Supplementary Table 2 | Agreement of association statistics with reference LD (s, susieR "
         "estimate_s_rss) at each of the 250 loci: EUR z with UK Biobank LD, EAS z with All of Us LD, trans-ancestry z with the merged matrix."),
        ("ST3_CARMA_per_locus", t3["per_locus"], "Supplementary Table 3a | Loci at which CARMA flagged at least one variant as inconsistent "
         "with the merged LD (seeded run): counts, s statistics, and whether the benchmark top variant was flagged."),
        ("ST3_CARMA_flagged_tops", t3["flagged_benchmark_leads"], "Supplementary Table 3b | The benchmark top variants CARMA flagged: each "
         "method's credible-set PIP at that variant and whether it is that method's top variant."),
        ("ST3_CARMA_all_flagged", t3["all_flagged_variants"], "Supplementary Table 3c | Every variant CARMA flagged."),
        ("ST4_high_confidence", t2, "Supplementary Table 4a | Credible-set calls at PIP > 0.5 and > 0.9, 193 loci, split by benchmark "
         "support (supported = benchmark PIP > 0.1)."),
        ("ST4_unsupported_calls", t2l, "Supplementary Table 4b | Every credible-set call at PIP > 0.5 that the benchmark does not support "
         "(benchmark PIP <= 0.1 or absent)."),
    ]
    contents = pd.DataFrame({"sheet": [s[0] for s in sheets], "title": [s[2] for s in sheets],
                             "rows": [len(s[1]) for s in sheets]})
    note = pd.DataFrame({"note": [
        "Primary analysis: strand-ambiguous (A/T, C/G) variants excluded; UK Biobank (GCTB) European and All of Us East "
        "Asian reference LD; benchmark = published in-sample SuSiEx (Yuan et al. 2024), 193 loci with a credible set.",
        "Methods: SuSiEx, MultiSuSiE, MESuSiE (ancestry-aware, per-ancestry LD); FINEMAP, SuSiE-RSS, RSparsePro, CARMA (trans-ancestry "
        "meta-analysis with one merged LD matrix). Variant IDs are PGC3 rsIDs; positions GRCh37."]})
    out = os.path.join(V2, "Supplementary_Tables.xlsx")
    with pd.ExcelWriter(out, engine="openpyxl") as w:
        # Layout (manuscript copy, 2026-10-01): one README sheet, contents table on top, notes four rows below it
        contents.to_excel(w, sheet_name="README", index=False)
        note.to_excel(w, sheet_name="README", index=False, header=False, startrow=len(contents) + 5)
        for name, df, title in sheets:
            pd.DataFrame({title: []}).to_excel(w, sheet_name=name, index=False)      # title in row 1
            df.to_excel(w, sheet_name=name, index=False, startrow=2)
        from openpyxl.styles import Alignment
        for name, _, _ in sheets:   # the one-cell table title: left-aligned, not pandas' centred header
            w.book[name]["A1"].alignment = Alignment(horizontal="left", vertical="top")
        autofit(w.book, skip_rows={name: {1} for name, _, _ in sheets} | {"README": set(range(len(contents) + 2, 10**6))})
    print(contents[["sheet", "rows"]].to_string(index=False))
    print(t1.iloc[:, :6].to_string(index=False))
    print(f"[done] {out}")


if __name__ == "__main__":
    main()
