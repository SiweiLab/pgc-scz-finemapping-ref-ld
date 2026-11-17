#!/usr/bin/env python3
"""Supplementary Fig. 9: SuSiEx with UK Biobank (GCTB-10k) vs All of Us European LD, East Asian LD (AoU) unchanged.

a  in-CS precision-recall against the benchmark (label benchmark PIP > 0.1; as Fig. 1b), both runs, AUPRC in the legend.
b  exact lead recovery (as Fig. 1c, each run gated on SuSiEx's own purity) on the 122 / 175 / 193 loci, paired bars, with
   in-CS false/true positives at PIP > 0.9 under the bars.
Both runs: no-palindromic, identical SuSiEx settings (results/susiex_eur_panel). No refits.
  python3 plot_supp_susiex_eur_panel.py --out output/supp_susiex_eur_panel.png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold"})
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import susiex_eur_panel_stats as S  # noqa: E402
from plot_cs_split import ARM_HUE  # noqa: E402
from plot_sensitivity_fdr import INK, INK_MUTED, GRID  # noqa: E402

STYLE = {"gctb": dict(color=ARM_HUE["SuSiEx-ours"], ls="-", label="UK Biobank (primary)"),
         "euraou": dict(color="#555555", ls="--", label="All of Us")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    res = {}
    for tag, _ in S.RUNS:
        p, c = S.load(tag)
        u = p[p.pb.notna()]
        r, pr = S.pr_points(u)
        k = u.incs & (u.pip > 0.9)
        res[tag] = dict(r=r, pr=pr, auprc=S.auprc(u), exact=S.exact(p, c),
                        fp=int((k & (u.pb <= 0.1)).sum()), tp=int((k & (u.pb > 0.1)).sum()))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 4.4), dpi=300, gridspec_kw=dict(width_ratios=[1.15, 1]))
    fig.subplots_adjust(wspace=0.32)

    # ---- a: precision-recall -------------------------------------------------------------------
    for tag, _ in S.RUNS:
        st = STYLE[tag]; R = res[tag]
        # draw the observed curve only; the zero-precision padding beyond reach is in the AUPRC, not on the plot
        n = len(R["r"]) - (2 if R["pr"][-1] == 0 else 0)
        axA.plot(R["r"][:n], R["pr"][:n], color=st["color"], ls=st["ls"], lw=1.8,
                 label=f"{st['label']}  {R['auprc']:.2f}")
    axA.set_xlim(0, 1.0); axA.set_ylim(0, 1.03)
    axA.set_xlabel("Recall", fontsize=10.6, color=INK); axA.set_ylabel("Precision", fontsize=10.6, color=INK)
    axA.set_title("Benchmark concordance, 193 loci", fontsize=10.6, color=INK, loc="left")
    axA.legend(title="European LD   AUPRC", title_fontproperties={"weight": "bold", "size": 8.6}, alignment="left", frameon=False, fontsize=8.6, labelcolor=INK, loc="center left", bbox_to_anchor=(0.02, 0.42))

    # ---- b: exact lead recovery, paired bars ------------------------------------------------------
    strata = [("122", "122 confident\nloci"), ("175", "175 single-CS\nloci"), ("193", "all 193\nloci")]
    w = 0.36
    for i, (st, lab) in enumerate(strata):
        nL = len(S.strata[st])
        for j, (tag, _) in enumerate(S.RUNS):
            v = 100 * res[tag]["exact"][st] / nL
            x = i + (j - 0.5) * w
            sty = STYLE[tag]
            axB.bar(x, v, width=w * 0.92, color=sty["color"] if tag == "gctb" else "white",
                    edgecolor=sty["color"], hatch=None if tag == "gctb" else "////", linewidth=1.2)
            axB.annotate(f"{v:.0f}%", (x, v + 1.2), ha="center", va="bottom", fontsize=8.8, color=INK)
    axB.set_xticks(range(len(strata))); axB.set_xticklabels([s[1] for s in strata], fontsize=9.2, color=INK)
    axB.set_ylim(0, 105); axB.set_ylabel("Exact recovery of the benchmark top variant", fontsize=10.6, color=INK)
    axB.yaxis.set_major_formatter(lambda v, _: f"{v:.0f}%")
    axB.set_title("Top-variant recovery", fontsize=10.6, color=INK, loc="left")
    fpt = "   ".join(f"{'UKB' if t == 'gctb' else 'AoU'}: {res[t]['fp']} FP / {res[t]['tp']} TP" for t, _ in S.RUNS)
    axB.annotate(f"in-CS calls at PIP > 0.9 —  {fpt}", (0.5, -0.24), xycoords="axes fraction", ha="center",
                 fontsize=8.6, color=INK_MUTED)
    h = [mpatches.Patch(facecolor=STYLE["gctb"]["color"], edgecolor=STYLE["gctb"]["color"], label="UK Biobank (primary)"),
         mpatches.Patch(facecolor="white", edgecolor=STYLE["euraou"]["color"], hatch="////", label="All of Us")]
    axB.legend(handles=h, frameon=False, fontsize=8.6, labelcolor=INK, loc="upper right", ncol=2)

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

    print({t: {k: v for k, v in res[t].items() if k not in ("r", "pr")} for t in res})
    print(f"[done] wrote {a.out}")


if __name__ == "__main__":
    main()
