#!/usr/bin/env python3
"""Fig. 1C: exact recovery of the published lead, per-ancestry LD vs merged LD.

122 confident loci, purity aligned at common_purity >= 0.5, exact match only -- the arm's top
CS variant must BE the published lead. No r2 proxy, no condition on our own PIP.

The comparison is between the two GROUPS. Each locus gives a per-ancestry score (fraction of
its 3 arms that recover the lead) and a merged score (fraction of its 4), and the two are
compared with a Wilcoxon signed-rank test paired BY LOCUS -- paired because every arm is
scored on the same loci, so the groups cannot be treated as independent samples.

common_purity is SuSiEx's composition (max over ancestries of the min pairwise |r|) applied to
all seven. It is the stricter of the two conventions in use, so it holds MultiSuSiE and
MESuSiE to a gate they would not apply to themselves; it is used because FINEMAP, CARMA and
RSparsePro have no comparable native gate and it is the only way to put all seven on one
footing.

Usage:
  python3 plot_fig1c_recovery.py --out figures/fig1c_recovery.png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                     "axes.titleweight": "bold"})
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


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--fig", default=V2)
    ap.add_argument("--stratum", default="122")
    ap.add_argument("--note", default="", help="sentence appended to the caption, e.g. which run the exports come from")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rate = pd.read_csv(os.path.join(args.fig, f"fig1c_{args.stratum}_E.csv")).set_index("arm")
    per = pd.read_csv(os.path.join(args.fig, f"fig1c_{args.stratum}_E_per_locus.csv"))
    n = int(rate.n_loci.iloc[0])
    M = per.pivot(index="locus", columns="arm", values="exact").astype(int)
    pa, mg = M[GROUPS[0][1]].mean(axis=1), M[GROUPS[1][1]].mean(axis=1)
    w = wilcoxon(pa, mg, zero_method="wilcox")
    gm = {GROUPS[0][0]: pa.mean(), GROUPS[1][0]: mg.mean()}

    fig, ax = plt.subplots(figsize=(8.6, 6.0), dpi=300)
    at = 0
    for gname, garms in GROUPS:
        for i, a in enumerate(garms):
            v = float(rate.loc[a, "recovery_rate"])
            ax.bar(at + i, v, width=0.66, color=ARM_HUE[a], alpha=0.95, linewidth=0,
                   zorder=3)
            ax.annotate(f"{v:.0%}", (at + i, v), xytext=(0, 4),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=12.5, color=INK)
        x0, x1 = at - 0.42, at + len(garms) - 1 + 0.42
        ax.plot([x0, x1], [gm[gname]] * 2, color=INK, lw=2.2, ls=(0, (4, 2)), zorder=6)
        ax.annotate(f"{gname}\nmean {gm[gname]:.1%}", ((x0 + x1) / 2, -0.185),
                    xycoords=("data", "axes fraction"), ha="center", va="top",
                    fontsize=12, color=INK)
        at += len(garms)
    ax.axvline(len(GROUPS[0][1]) - 0.5, color=GRID, lw=1.3, zorder=1)

    # significance bracket between the two group means, placed above the tallest bar AND its % label
    # (a fixed y = 0.90 ran through the labels once v2's bars passed 0.9)
    vmax = max(float(rate.loc[a, "recovery_rate"]) for _, g in GROUPS for a in g)
    y = max(0.90, vmax + 0.075)
    xa = (len(GROUPS[0][1]) - 1) / 2
    xb = len(GROUPS[0][1]) + (len(GROUPS[1][1]) - 1) / 2
    ax.plot([xa, xa, xb, xb], [y, y + .035, y + .035, y], color=INK, lw=1.6, zorder=6)
    # only the stars on the plot; the test, P and effect size live in the caption
    ax.annotate(stars(w.pvalue), ((xa + xb) / 2, y + .045), ha="center", va="bottom",
                fontsize=20, color=INK)

    ax.set_xticks(range(len(ARMS)))
    ax.set_xticklabels([disp(a) for a in ARMS], fontsize=12, color=INK, rotation=20,
                       ha="right")
    ax.set_xlim(-0.6, len(ARMS) - 0.4)
    ax.set_ylim(0, max(1.04, y + 0.13))
    ax.set_ylabel("Exact recovery of the benchmark top variant", fontsize=14.5, color=INK)
    ax.set_title(f"{n} loci, single benchmark credible set, top PIP > 0.1",
                 fontsize=13.5, color=INK, loc="left", pad=10)
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, labelsize=13.5)
    ax.yaxis.set_major_formatter(lambda v, _: f"{v:.0%}")

    note = (
        f"{n} confident loci; purity aligned at common_purity >= 0.5. Exact match only: the "
        "arm's top CS variant must BE the published lead, with no r2 proxy and no condition "
        "on our own PIP. An arm that does not resolve exactly one CS is non-concordant. "
        "Dashed lines are group means. The test compares the two GROUPS: each locus gives a "
        "per-ancestry score (fraction of its 3 arms recovering the lead) and a merged score "
        "(fraction of its 4), compared by Wilcoxon signed-rank paired by locus, since every "
        f"arm is scored on the same loci: per-ancestry mean {pa.mean():.1%} against merged "
        f"{mg.mean():.1%}, mean per-locus difference +{(pa - mg).mean():.1%}, "
        f"W = {w.statistic:.0f}, P = {w.pvalue:.1e} ({stars(w.pvalue)}); "
        f"{int((pa > mg).sum())} loci favour per-ancestry, {int((pa < mg).sum())} favour "
        f"merged, {int((pa == mg).sum())} are tied. common_purity is SuSiEx's composition applied to all "
        "seven -- the stricter of the two conventions, used because FINEMAP, CARMA and "
        "RSparsePro have no comparable native gate.")
    if args.note:
        note = note.rstrip() + " " + args.note
    cp = args.out.rsplit(".", 1)[0] + "_caption.txt"
    with open(cp, "w") as fh:
        fh.write(" ".join(note.split()) + "\n")
    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    print(f"[done] wrote {args.out} and {cp}\n")
    print(f"per-ancestry mean {pa.mean():.3f} | merged mean {mg.mean():.3f} | "
          f"diff {(pa - mg).mean():+.3f}")
    print(f"loci favouring per-ancestry {int((pa > mg).sum())}, merged "
          f"{int((pa < mg).sum())}, tied {int((pa == mg).sum())}")
    print(f"Wilcoxon W={w.statistic:.0f}, P={w.pvalue:.3g}")


if __name__ == "__main__":
    main()
