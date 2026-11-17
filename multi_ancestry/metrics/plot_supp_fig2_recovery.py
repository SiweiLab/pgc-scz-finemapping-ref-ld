#!/usr/bin/env python3
"""Supp. Fig. 2: exact recovery across the three strata, with and without a purity gate.

Fig. 1C is the headline -- one stratum, one gate, the group test. This is the full grid
behind it, so the numbers quoted in the text can be traced to a stratum and a gate:

  rows      scenario C, each arm's own native purity
            scenario E, purity aligned at common_purity >= 0.5
  columns   193 loci (benchmark formed a CS) / 175 (single CS) / 122 (single CS, max PIP > 0.1)

Exact match throughout: the arm's top CS variant must BE the published lead, with no r2
proxy and no condition on our own PIP. An arm that does not resolve exactly one CS is
non-concordant.

Group means are drawn as dashed lines so the per-ancestry vs merged gap can be read off each
panel; the formal test is in Fig. 1C and is not repeated six times here.

Usage:
  python3 plot_supp_fig2_recovery.py --out figures/supp_fig2_recovery_strata.png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                     "axes.titleweight": "bold"})
import matplotlib.patches as mpatches
from paths import V2
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_cs_split import ARM_HUE, disp
from plot_sensitivity_fdr import INK, INK_MUTED, GRID

GROUPS = [("per-ancestry LD", ["SuSiEx-ours", "MultiSuSiE", "MESuSiE"]),
          ("merged LD", ["FINEMAP-merged", "SuSiE-merged", "RSparsePro", "CARMA"])]
ARMS = [a for _, g in GROUPS for a in g]
STRATA = ["193", "175", "122"]
SUB = {"193": "benchmark formed a CS", "175": "single benchmark CS",
       "122": "single CS, top PIP > 0.1"}
SCEN = [("C", "each method's own purity filter", 0.95), ("E", "common purity $\\geq$ 0.5", 0.40)]



def gate_effect():
    """The C -> E change per arm, in loci, from the tables actually plotted (was a hard-coded v1 sentence)."""
    d = {}
    for st in STRATA:
        c = pd.read_csv(os.path.join(ARGS.fig, f"fig1c_{st}_C.csv")).set_index("arm").exact
        e = pd.read_csv(os.path.join(ARGS.fig, f"fig1c_{st}_E.csv")).set_index("arm").exact
        for a in ARMS:
            d.setdefault(a, []).append(int(e[a] - c[a]))
    sg = lambda x: "0" if x == 0 else f"{x:+d}"
    rng = lambda v: sg(min(v)) if min(v) == max(v) else f"{sg(min(v))} to {sg(max(v))}"
    parts = [f"{disp(a)} {rng(v)}" for a, v in d.items()]
    return ("Effect of the purity gate (E minus C, exact matches, across the three strata): " + "; ".join(parts) +
            ". FINEMAP-merged is the arm whose native output is ungated.")

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--fig", default=V2)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    global ARGS
    ARGS = args

    fig, axes = plt.subplots(len(SCEN), len(STRATA),
                             figsize=(5.6 * len(STRATA), 5.6 * len(SCEN)), dpi=300,
                             sharey=True)
    tbl = []
    for r, (sc, sclab, alpha) in enumerate(SCEN):
        for c, st in enumerate(STRATA):
            ax = axes[r, c]
            d = pd.read_csv(os.path.join(args.fig, f"fig1c_{st}_{sc}.csv")).set_index("arm")
            per = pd.read_csv(os.path.join(args.fig, f"fig1c_{st}_{sc}_per_locus.csv"))
            M = per.pivot(index="locus", columns="arm", values="exact").astype(int)
            n = int(d.n_loci.iloc[0])
            at = 0
            for gname, garms in GROUPS:
                for i, a in enumerate(garms):
                    v = float(d.loc[a, "recovery_rate"])
                    ax.bar(at + i, v, width=0.68, color=ARM_HUE[a], alpha=alpha,
                           linewidth=0, zorder=3)
                    ax.annotate(f"{v:.0%}", (at + i, v), xytext=(0, 3),
                                textcoords="offset points", ha="center", va="bottom",
                                fontsize=10.5, color=INK)
                gm = M[garms].mean(axis=1).mean()
                ax.plot([at - 0.44, at + len(garms) - 1 + 0.44], [gm] * 2, color=INK,
                        lw=2.0, ls=(0, (4, 2)), zorder=6)
                # group name and mean under the group, below the tick labels, as in Fig. 1c
                ax.annotate(f"{gname}\nmean {gm:.1%}", (at + (len(garms) - 1) / 2, -0.205),
                            xycoords=("data", "axes fraction"), ha="center", va="top",
                            fontsize=10.5, color=INK)
                tbl.append(dict(stratum=st, n_loci=n, scenario=sc, group=gname,
                                group_mean=gm))
                at += len(garms)
            # the same per-ancestry vs merged test as Fig. 1C, marked per panel
            g1 = M[GROUPS[0][1]].mean(axis=1)
            g2 = M[GROUPS[1][1]].mean(axis=1)
            w = wilcoxon(g1, g2, zero_method="wilcox")
            st_ = ("***" if w.pvalue < 1e-3 else "**" if w.pvalue < 1e-2
                   else "*" if w.pvalue < 0.05 else "n.s.")
            yb = 1.0   # above every bar AND its % label (v2 bars reach 0.92); was 0.93
            xa, xb = (len(GROUPS[0][1]) - 1) / 2, len(GROUPS[0][1]) + (len(GROUPS[1][1]) - 1) / 2
            ax.plot([xa, xa, xb, xb], [yb, yb + .03, yb + .03, yb], color=INK, lw=1.5,
                    zorder=6)
            ax.annotate(st_, ((xa + xb) / 2, yb + .038), ha="center", va="bottom",
                        fontsize=16, color=INK)
            tbl[-1]["wilcoxon_P"] = w.pvalue
            tbl[-2]["wilcoxon_P"] = w.pvalue
            ax.axvline(len(GROUPS[0][1]) - 0.5, color=GRID, lw=1.2, zorder=1)
            ax.set_xticks(range(len(ARMS)))
            ax.set_xticklabels([disp(a) for a in ARMS],
                               fontsize=11, color=INK, rotation=28, ha="right")
            ax.set_xlim(-0.6, len(ARMS) - 0.4)
            ax.set_ylim(0, 1.15)
            ax.set_yticks([0, .2, .4, .6, .8, 1.0])
            ax.set_title(f"{n} loci — {SUB[st]}\n{sclab}", fontsize=12.5,
                         color=INK, loc="left", pad=8)
            ax.grid(False)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            for s in ("left", "bottom"):
                ax.spines[s].set_color(GRID)
            ax.tick_params(colors=INK_MUTED, labelsize=11.5)
            ax.yaxis.set_major_formatter(lambda v, _: f"{v:.0%}")
        axes[r, 0].set_ylabel("Exact top-variant recovery", fontsize=12.5,
                              color=INK)

    h = [mpatches.Patch(facecolor=INK_MUTED, alpha=a, linewidth=0, label=lab)
         for s, lab, a in SCEN]
    h.append(plt.Line2D([], [], color=INK, lw=2.0, ls=(0, (4, 2)), label="group mean"))
    fig.legend(handles=h, frameon=False, fontsize=12, ncol=3, labelcolor=INK,
               loc="lower center", bbox_to_anchor=(0.5, -0.06))

    note = (
        "Exact match throughout: the arm's top CS variant must BE the published lead, with "
        "no r2 proxy and no condition on our own PIP; an arm that does not resolve exactly "
        "one CS is non-concordant. Columns are the three strata, rows the two purity "
        "settings. C uses each arm's own native gate; E aligns all seven at "
        "common_purity >= 0.5, which is SuSiEx's composition (max over ancestries of the "
        "min pairwise |r|) and the stricter of the two conventions -- it is used because "
        "FINEMAP, CARMA and RSparsePro have no comparable native gate. Dashed lines are "
        "group means; the per-ancestry vs merged comparison is "
        "marked per panel as stars from a Wilcoxon signed-rank test paired by locus, the same test as Fig. 1C, with exact P values in the companion supp_fig2_recovery_strata_group_means.csv. " + gate_effect())
    cp = args.out.rsplit(".", 1)[0] + "_caption.txt"
    with open(cp, "w") as fh:
        fh.write(" ".join(note.split()) + "\n")
    fig.tight_layout()
    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    print(f"[done] wrote {args.out} and {cp}\n")
    T = pd.DataFrame(tbl)
    T.to_csv(args.out.rsplit(".", 1)[0] + "_group_means.csv", index=False)
    print(T.pivot(index=["stratum", "n_loci"], columns=["scenario", "group"],
                  values="group_mean").round(3).to_string())


if __name__ == "__main__":
    main()
