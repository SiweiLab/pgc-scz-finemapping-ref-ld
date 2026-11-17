#!/usr/bin/env python3
"""Supplementary figure: FINEMAP and SuSiE-RSS with merged EUR+EAS LD on the trans-ancestry statistics, against the
same two methods on EUR statistics with UK Biobank LD alone. It separates adding EAS data (as merged) from merging.

a  in-CS precision-recall, 193 loci, label benchmark PIP > 0.1 -- scored by plot_sensfdr_curve.py itself (run here for
   these four arms), so the merged AUPRCs equal Fig. 1b's.
b  exact lead recovery from the Fig. 1c exports (common_purity >= 0.5, scenario E) on 122 / 175 / 193 loci, with in-CS
   FP / TP at PIP > 0.9 from fp_in_cs.csv underneath. Layout as Supplementary Fig. 6.
  python3 plot_supp_eur_only.py --out output/supp_eur_only.png
"""
import argparse
import os
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold"})
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plot_cs_split import ARM_HUE, disp  # noqa: E402
from plot_sensitivity_fdr import INK, INK_MUTED, GRID  # noqa: E402

from paths import V2  # noqa: E402
# (method, merged arm, EUR-only arm)
PAIRS = [("FINEMAP", "FINEMAP-merged", "EUR-FINEMAP"), ("SuSiE-RSS", "SuSiE-merged", "EUR-SuSiE-RSS")]
FP_NAME = {a: disp(a) for a in ("FINEMAP-merged", "EUR-FINEMAP", "SuSiE-merged", "EUR-SuSiE-RSS")}  # fp_in_cs.csv uses display names


def curves():
    with tempfile.TemporaryDirectory() as d:
        png = os.path.join(d, "pr.png")
        subprocess.run([sys.executable, os.path.join(HERE, "plot_sensfdr_curve.py"), "--pr", "--strata", "193",
                        "--export", os.path.join(HERE, "..", "results", "finemap", "main"), "--suffix=-noambig",
                        "--exclude", "SuSiEx-ours,MultiSuSiE,MESuSiE,RSparsePro,CARMA", "--out", png],
                       check=True, capture_output=True)
        return pd.read_csv(png.rsplit(".", 1)[0] + "_curve.csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    C = curves()
    auprc = C.groupby("arm").auprc.first().to_dict()
    rec = {st: pd.read_csv(os.path.join(V2, f"fig1c_{st}_E.csv")).set_index("arm") for st in ("122", "175", "193")}
    fp = pd.read_csv(os.path.join(V2, "fp_in_cs.csv")).set_index("arm")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 4.4), dpi=300, gridspec_kw=dict(width_ratios=[1.15, 1]))
    fig.subplots_adjust(wspace=0.32)

    # ---- a: precision-recall -------------------------------------------------------------------
    for meth, mg, eu in PAIRS:
        hue = ARM_HUE[mg]
        for arm, ls, lab in ((mg, "-", "merged EUR+EAS"), (eu, "--", "EUR only")):
            g = C[C.arm == arm].sort_values("s")
            m = g.precision.notna()
            axA.plot(g.sensitivity[m], g.precision[m], color=hue, ls=ls, lw=1.8,
                     label=f"{meth}, {lab}  {auprc[arm]:.2f}")
    axA.set_xlim(0, 1.0); axA.set_ylim(0, 1.03)
    axA.set_xlabel("Recall", fontsize=10.6, color=INK); axA.set_ylabel("Precision", fontsize=10.6, color=INK)
    axA.set_title("Benchmark concordance, 193 loci", fontsize=10.6, color=INK, loc="left")
    axA.legend(title="Method, LD   AUPRC", title_fontproperties={"weight": "bold", "size": 8.4}, alignment="left",
               frameon=False, fontsize=8.4, labelcolor=INK, loc="lower left")

    # ---- b: exact lead recovery ------------------------------------------------------------------
    strata = [("122", "122 confident\nloci"), ("175", "175 single-CS\nloci"), ("193", "all 193\nloci")]
    order = [(mg, False) for _, mg, _ in PAIRS[:1]] + [(PAIRS[0][2], True)] + [(PAIRS[1][1], False), (PAIRS[1][2], True)]
    w = 0.2
    for i, (st, _) in enumerate(strata):
        for j, (arm, eur) in enumerate(order):
            v = 100 * rec[st].at[arm, "exact_rate"]
            x = i + (j - 1.5) * w
            hue = ARM_HUE[arm]
            axB.bar(x, v, width=w * 0.92, color="white" if eur else hue, edgecolor=hue,
                    hatch="////" if eur else None, linewidth=1.2)
            axB.annotate(f"{v:.0f}", (x, v + 1.2), ha="center", va="bottom", fontsize=7.8, color=INK)
    axB.set_xticks(range(len(strata))); axB.set_xticklabels([s[1] for s in strata], fontsize=9.2, color=INK)
    axB.set_ylim(0, 100); axB.set_ylabel("Exact recovery of the benchmark top variant", fontsize=10.6, color=INK)
    axB.yaxis.set_major_formatter(lambda v, _: f"{v:.0f}%")
    axB.set_title("Top-variant recovery", fontsize=10.6, color=INK, loc="left")
    fpt = "   ".join(f"{m} {l}: {fp.at[FP_NAME[arm], 'FP_pip>0.9']}/{fp.at[FP_NAME[arm], 'TP_pip>0.9']}"
                     for m, mg, eu in PAIRS for arm, l in ((mg, "merged"), (eu, "EUR")))
    axB.annotate(f"in-CS FP/TP at PIP > 0.9 —  {fpt}", (0.5, -0.24), xycoords="axes fraction", ha="center",
                 fontsize=7.8, color=INK_MUTED)
    h = [mpatches.Patch(facecolor=INK_MUTED, edgecolor=INK_MUTED, label="merged EUR+EAS"),
         mpatches.Patch(facecolor="white", edgecolor=INK_MUTED, hatch="////", label="EUR only")] + \
        [mpatches.Patch(facecolor=ARM_HUE[mg], label=meth) for meth, mg, _ in PAIRS]
    axB.legend(handles=h, frameon=False, fontsize=8.2, labelcolor=INK, loc="upper right", ncol=2)

    for ax in (axA, axB):
        ax.grid(False)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=INK_MUTED, labelsize=9.0)
    for ax, lab in ((axA, "a"), (axB, "b")):
        ax.annotate(lab, (-0.13, 1.06), xycoords="axes fraction", fontsize=16, fontweight="bold", color=INK)
    fig.savefig(a.out, bbox_inches="tight", facecolor="white")

    r = {arm: rec["122"].at[arm, "exact_rate"] for _, mg, eu in PAIRS for arm in (mg, eu)}
    note = ("FINEMAP and SuSiE-RSS run on the trans-ancestry meta-analysis statistics with one LD matrix merged from the "
            "UK Biobank European and All of Us East Asian references (merged EUR+EAS), against the same methods and "
            "settings on the European statistics with UK Biobank LD alone (EUR only); A/T and C/G variants excluded. "
            "(a) Precision-recall against the benchmark across the 193 loci, scoring credible-set members only, with the "
            "benchmark's variants at PIP > 0.1 as positives, as in Fig. 1b. (b) Exact recovery of the benchmark's lead "
            "variant, as in Fig. 1c (exactly one credible set, purity >= 0.5, top variant = benchmark lead), on the 122 "
            "confident, 175 single-credible-set and all 193 loci. Below: credible-set calls at PIP > 0.9 whose "
            "benchmark PIP is at most 0.1 (FP) or above 0.1 (TP). On the 122 loci, merging in East Asian data raised "
            f"recovery from {100 * r['EUR-FINEMAP']:.0f}% to {100 * r['FINEMAP-merged']:.0f}% (FINEMAP) and from "
            f"{100 * r['EUR-SuSiE-RSS']:.0f}% to {100 * r['SuSiE-merged']:.0f}% (SuSiE-RSS), still below the "
            "ancestry-aware methods (83-92%; Fig. 1c).")
    open(a.out.rsplit(".", 1)[0] + "_caption.txt", "w").write(note + "\n")
    print(auprc); print(r)
    print(f"[done] wrote {a.out}")


if __name__ == "__main__":
    main()
