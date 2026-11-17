#!/usr/bin/env python3
"""Supplementary Fig. 8a, b: susieR kriging (kriging_rss) at the Fig. 2b locus, EUR (GCTB-10k) and EAS (AoU) separately.

Reads analyses/kriging.R's tables (results/kriging/<L>_{EUR,EAS}_kriging.tsv): SuSiEx's exact no-palindromic inputs, n = 58,749 /
14,425. Rows = ancestry; columns = observed z against the LD-predicted z, and the standardized difference along the locus.
susieR's own outlier rule (logLR > 2 and |z| > 2) is drawn in red. logLR compares the sign-flipped z with the observed z
(kriging_rss), so it is an allele-switch test: a variant far from its prediction but with the same sign, like
12:29935396_C_A in EAS, is not flagged; the variants SuSiEx placed at PIP > 0.9 are marked.
s is not printed here; kriging.R's s (estimate_s_rss with n) equals s_per_ancestry.py's, which the Methods and Supp Table 2 quote.
  python3 plot_supp_kriging.py --out output/supp_kriging_chr12_29898169.png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold"})
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plot_sensitivity_fdr import INK, INK_MUTED, GRID  # noqa: E402

KRIGING = os.path.join(HERE, "..", "results", "kriging")
OK_C, OUT_C, HL_C = "#9ec5f4", "#e34948", "#1c5cab"
ANC = [("EUR", "European: EUR z with UK Biobank (GCTB) LD"), ("EAS", "East Asian: EAS z with All of Us LD")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--locus", default="chr12_29898169_29948269")
    ap.add_argument("--highlight", default="12:29935396_C_A,rs654825")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    hl = [x for x in a.highlight.split(",") if x]

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.4), dpi=300, gridspec_kw=dict(width_ratios=[1, 1.35]))
    fig.subplots_adjust(hspace=0.42, wspace=0.26)
    rows = []
    for i, (anc, title) in enumerate(ANC):
        d = pd.read_csv(os.path.join(KRIGING, f"{a.locus}_{anc}_kriging.tsv"), sep="\t", dtype={"snp": str})
        d["outlier"] = d.outlier.astype(str).str.upper() == "TRUE"
        lim = float(np.nanmax(np.abs(np.r_[d.z, d.condmean]))) * 1.08
        for j, (X, Y) in enumerate((("condmean", "z"), ("bp", "z_std_diff"))):
            ax = axes[i, j]
            if j == 0:
                ax.plot([-lim, lim], [-lim, lim], color=INK_MUTED, lw=0.8, zorder=1)
            else:
                ax.axhline(0, color=INK_MUTED, lw=0.7, zorder=1)
                ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1e6:.2f}"))
            ok = d[~d.outlier]
            ax.scatter(ok[X], ok[Y], s=14, color=OK_C, linewidths=0, zorder=2)
            o = d[d.outlier]
            ax.scatter(o[X], o[Y], s=24, color=OUT_C, linewidths=0, zorder=3)
            h = d[d.snp.isin(hl)]
            ax.scatter(h[X], h[Y], s=80, marker="D", facecolors="none", edgecolors=HL_C, linewidths=1.5, zorder=4)
            for k, (_, r) in enumerate(h.sort_values(X).iterrows()):
                ax.annotate(r.snp, (r[X], r[Y]), xytext=(7, 9 if k % 2 == 0 else -13), textcoords="offset points",
                            fontsize=7.8, color=HL_C)
            if j == 0:
                ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
                ax.set_xlabel("z predicted from neighbours via LD", fontsize=9.4, color=INK)
                ax.set_ylabel("Observed z", fontsize=9.4, color=INK)
            else:
                ax.set_xlabel("chr12 position (Mb, GRCh37)", fontsize=9.4, color=INK)
                ax.set_ylabel("Standardized difference", fontsize=9.4, color=INK)
            ax.grid(color=GRID, lw=0.5, zorder=0); ax.set_axisbelow(True)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            ax.tick_params(colors=INK_MUTED, labelsize=8.4)
        axes[i, 0].set_title(f"{title}\n{len(d)} variants; susieR allele-switch flags: {int(d.outlier.sum())}", fontsize=9.6,
                             color=INK, loc="left")
        for _, r in d[d.snp.isin(hl)].iterrows():
            rows.append(dict(ancestry=anc, rsid=r.snp, z=round(r.z, 2), predicted=round(r.condmean, 2),
                             std_diff=round(r.z_std_diff, 2), logLR=round(r.logLR, 2), outlier=bool(r.outlier)))
    for ax, lab in ((axes[0, 0], "a"), (axes[1, 0], "b")):
        ax.annotate(lab, (-0.24, 1.18), xycoords="axes fraction", fontsize=16, fontweight="bold", color=INK)
    h = [plt.Line2D([], [], marker="o", ls="", color=OK_C, label="not flagged"),
         plt.Line2D([], [], marker="o", ls="", color=OUT_C, label="susieR allele-switch flag (logLR > 2 and |z| > 2)"),
         plt.Line2D([], [], marker="D", ls="", markerfacecolor="none", markeredgecolor=HL_C,
                    label="SuSiEx PIP > 0.9, unsupported by the benchmark")]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=8.6, bbox_to_anchor=(0.5, -0.02))
    fig.savefig(a.out, bbox_inches="tight", facecolor="white")
    t = pd.DataFrame(rows)
    t.to_csv(a.out.rsplit(".", 1)[0] + ".csv", index=False)
    print(t.to_string(index=False))
    print(f"[done] wrote {a.out}")


if __name__ == "__main__":
    main()
