#!/usr/bin/env python3
"""Supp. Fig. 3a/3b: the L sweep as sensitivity/FDR curves, one panel per arm.

Same metric as Fig. 1B -- the truth set is fixed at benchmark PIP > 0.1 and only the arm's
own threshold sweeps, 0.01 to 0.99 -- so a reader who has read Fig. 1B already knows how to
read this. What changes is that each arm appears four times, once per L, in its own hue at
four opacities. Faceting by arm rather than drawing 28 curves on one axes is the only way
this stays readable; the comparison the figure supports is within an arm, across L.

  2a  in-CS calls only   (matches Fig. 1B)
  2b  all calls, in and out of credible sets

SuSiEx appears as SuSiEx (multi-step off): it ignores --n_sig entirely when mult-step is on,
so the sweep had to disable it. That arm is therefore NOT the primary run's SuSiEx and its
L=5 point is not the primary point.

Usage:
  python3 plot_Lsweep_curve.py --out figures/supp_fig2a_Lsweep_curve_inCS.png
  python3 plot_Lsweep_curve.py --all-variants \
      --out figures/supp_fig2b_Lsweep_curve_allvar.png
"""
import argparse
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                     "axes.titleweight": "bold"})
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_sensfdr import bench_cs_loci, BENCH, BENCH_RAW
from sensitivity_fdr import load_pips, variant_key
from plot_cs_split import ARM_HUE
from plot_sensitivity_fdr import INK, INK_MUTED, GRID

L_VALUES = [1, 3, 5, 10]
ALPHA = {1: 0.30, 3: 0.52, 5: 0.74, 10: 0.96}
HUE_OF = {"SuSiEx-nomultistep": ARM_HUE["SuSiEx-ours"], "MultiSuSiE": ARM_HUE["MultiSuSiE"],
          "MESuSiE": ARM_HUE["MESuSiE"], "FINEMAP-merged": ARM_HUE["FINEMAP-merged"],
          "SuSiE-merged": ARM_HUE["SuSiE-merged"], "RSparsePro": ARM_HUE["RSparsePro"],
          "CARMA": ARM_HUE["CARMA"]}
ORDER = ["SuSiEx-nomultistep", "MultiSuSiE", "MESuSiE", "FINEMAP-merged", "SuSiE-merged",
         "RSparsePro", "CARMA"]
DISPLAY = {"SuSiEx-nomultistep": "SuSiEx (multi-step off)", "FINEMAP-merged": "FINEMAP", "SuSiE-merged": "SuSiE-RSS"}
PER_ANC = ["SuSiEx-nomultistep", "MultiSuSiE", "MESuSiE"]
LEGEND = "__legend__"
MARKS = [0.1, 0.5, 0.9]
MARK_STYLE = {0.1: ("o", 16, False), 0.5: ("o", 24, True), 0.9: ("s", 18, False)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--export", default=os.path.join(here, "..", "results", "finemap", "L_sweep"))
    ap.add_argument("--all-variants", action="store_true")
    ap.add_argument("--label-pip", type=float, default=0.1)
    ap.add_argument("--suffix", default="-noambig", help="v2 label suffix, e.g. --suffix=-noambig")
    ap.add_argument("--primary-export", default=None,
                    help="v2: the main-run export dir; each arm's PRIMARY L (10, FINEMAP-merged k = 5) is read from here, "
                         "since the sweep reuses it rather than refitting")
    ap.add_argument("--pr", action="store_true", help="precision-recall orientation, AUPRC per L (as Fig. 1B)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    keep, _ = bench_cs_loci(BENCH_RAW)
    b = load_pips(BENCH)
    b = b[b.locus.isin(keep)][["locus", "key", "pip"]].rename(columns={"pip": "pb"})
    b = b.groupby(["locus", "key"], as_index=False).pb.max()

    found, primary_L = {}, {}
    for f in sorted(glob.glob(os.path.join(args.export, "primary_*_cs_purity.csv.gz"))):
        m = re.match(r"primary_(.+)_L(\d+)_cs_purity\.csv\.gz$", os.path.basename(f))
        if m:
            arm = m.group(1)[:-len(args.suffix)] if args.suffix and m.group(1).endswith(args.suffix) else m.group(1)
            found[(arm, int(m.group(2)))] = f
    if args.primary_export:  # the primary fit IS one sweep point for every arm except SuSiEx (mult-step off here)
        for arm, L in (("MultiSuSiE", 10), ("MESuSiE", 10), ("SuSiE-merged", 10), ("RSparsePro", 10), ("CARMA", 10),
                       ("FINEMAP-merged", 5)):
            f = os.path.join(args.primary_export, f"primary_{arm}{args.suffix}_cs_purity.csv.gz")
            if os.path.exists(f) and (arm, L) not in found:
                found[(arm, L)] = f
                primary_L[arm] = L
    arms = [a for a in ORDER if any(k[0] == a for k in found)]
    grid = np.unique(np.concatenate([np.arange(0.01, 1.0, 0.01), np.array(MARKS)]))

    rows = []
    for arm in arms:
        for L in L_VALUES:
            f = found.get((arm, L))
            if f is None:
                print(f"[warn] {arm} L={L} missing, skipped")
                continue
            d = pd.read_csv(f, low_memory=False)
            d = d[d.locus.isin(keep)].copy()
            d["key"] = variant_key(d)
            u = b.merge(d[["locus", "key", "pip", "cs_ids_all"]], on=["locus", "key"],
                        how="inner")
            pos = (u.pb > args.label_pip).to_numpy()
            nP = int(pos.sum())
            incs = u.cs_ids_all.notna().to_numpy()
            pc = u.pip.to_numpy()
            for s in grid:
                call = (pc > s) if args.all_variants else ((pc > s) & incs)
                nC = int(call.sum()); tp = int((call & pos).sum())
                rows.append(dict(arm=arm, L=L, s=float(s), n_calls=nC, tp=tp, n_PB=nP,
                                 sensitivity=tp / nP if nP else np.nan,
                                 fdr=(nC - tp) / nC if nC else np.nan,
                                 precision=tp / nC if nC else np.nan))
    R = pd.DataFrame(rows)

    def _auprc(g):  # Fig. 1B's convention: precision 0 beyond the arm's reach, trapezoid over recall in [0, 1]
        g = g[g.precision.notna() & g.sensitivity.notna()]
        x, y = g.sensitivity.to_numpy(), g.precision.to_numpy()
        if len(x) < 2:
            return np.nan
        o = np.argsort(x); x, y = x[o], y[o]
        return float(np.trapezoid(np.r_[y[0], y, 0.0], np.r_[0.0, x, x[-1]]))
    AU = {(a, L): _auprc(R[(R.arm == a) & (R.L == L)]) for a in arms for L in L_VALUES}
    R["auprc"] = [AU[(r.arm, r.L)] for r in R.itertuples()]
    ycol = "precision" if args.pr else "fdr"

    # row 1: the ancestry-aware methods, then the L / PIP key; row 2: the merged-reference-LD methods
    ncol = 4
    slots = [a for a in PER_ANC if a in arms] + [LEGEND] + [a for a in arms if a not in PER_ANC]
    nrow = int(np.ceil(len(slots) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 4.3 * nrow), dpi=300,
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, arm in zip(axes, slots):
        if arm == LEGEND:
            continue
        for L in L_VALUES:
            g = R[(R.arm == arm) & (R.L == L)].sort_values("s")
            if not len(g):
                continue
            m = g[ycol].notna()
            ax.plot(g.sensitivity[m], g[ycol][m], color=HUE_OF[arm], alpha=ALPHA[L], lw=2.0,
                    zorder=3, solid_capstyle="round")
            for sv in MARKS:
                r = g[np.isclose(g.s, sv)]
                if len(r) and pd.notna(r[ycol].iloc[0]):
                    mk, sz, fc = MARK_STYLE[sv]
                    ax.scatter(r.sensitivity, r[ycol], s=sz, marker=mk, alpha=ALPHA[L],
                               facecolor=HUE_OF[arm] if fc else "white",
                               edgecolor=HUE_OF[arm], linewidths=1.1, zorder=5)
        ax.set_title(DISPLAY.get(arm, arm), fontsize=13, color=INK, loc="left", pad=8)
        ax.set_xlim(0, 1.0); ax.set_ylim(0, 1.0)
        ax.grid(False)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=INK_MUTED, labelsize=12, labelbottom=True, labelleft=True)
        ax.set_xlabel("Recall" if args.pr else "Sensitivity", fontsize=13, color=INK)
        ax.set_ylabel("Precision" if args.pr else "FDR", fontsize=13, color=INK)
        if args.pr:
            txt = "AUPRC\n" + "\n".join(f"L={L}{' (primary)' if primary_L.get(arm) == L else ''}  {AU[(arm, L)]:.2f}"
                                         for L in L_VALUES if np.isfinite(AU.get((arm, L), np.nan)))
            ax.text(0.03, 0.03, txt, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.5, color=INK)
    for ax in axes[len(slots):]:
        ax.set_visible(False)

    hL = [mlines.Line2D([], [], color=INK_MUTED, alpha=ALPHA[L], lw=2.4, label=f"L = {L}")
          for L in L_VALUES]
    hM = [mlines.Line2D([], [], color=INK_MUTED, marker=MARK_STYLE[t][0], ls="",
                        markerfacecolor=INK_MUTED if MARK_STYLE[t][2] else "white",
                        markeredgewidth=1.2, markersize=6 if MARK_STYLE[t][2] else 5,
                        label=f"PIP > {t:g}") for t in MARKS]
    spare = axes[slots.index(LEGEND)]
    if spare is not None:
        spare.set_visible(True)
        spare.axis("off")
        # the FIRST legend is the one that needs add_artist; the second call replaces it
        lg1 = spare.legend(handles=hL, frameon=False, fontsize=13, labelcolor=INK,
                           loc="upper left", bbox_to_anchor=(0.02, 0.98),
                           alignment="left", title="L", title_fontsize=13)
        spare.add_artist(lg1)
        spare.legend(handles=hM, frameon=False, fontsize=12.5, labelcolor=INK,
                     loc="upper left", bbox_to_anchor=(0.02, 0.52), alignment="left")

    fig.tight_layout()
    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    print(f"[done] wrote {args.out}\n")
    R.to_csv(args.out.rsplit(".", 1)[0] + "_curve.csv", index=False)
    print(R[R.s == 0.1].pivot(index="arm", columns="L",
                              values=["sensitivity", "fdr"]).reindex(arms).round(3)
          .to_string())


if __name__ == "__main__":
    main()
