#!/usr/bin/env python3
"""Supplement: does the EUR reference panel change the answer? Three panels, B53.

GCTB-10k (UK Biobank subsample), AoU-267k (All of Us EUR, in-sample-style hard calls) and
PolyFun/UKBB-337k (British-ancestry UK Biobank; see ESTABLISHED_FACTS sec 25 on the label).
Every fit uses the same EUR z over the same non-palindromic variants per locus (the common set),
single-ancestry susie_rss at the primary EUR-SuSiE-RSS parameters. The panel is the only variable.

  a  distribution of the s statistic per panel (ECDF, log x). Lower s = z more consistent with that LD.
  b  s per locus, paired against GCTB-10k; above the diagonal = that panel fits the locus worse.
  c-e  PIP agreement for each pair of panels, variant by variant.

Usage:
  python3 plot_eur_panels_b53.py --b53-dir ../aou_ld/task_b53_out --out figures/eur_panels_b53.png
"""
import argparse
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INK, INK_MUTED, GRID, SHADE = "#1a1a19", "#6b6a63", "#e3e2dc", "#f1f0ec"
# s from B71 (real susieR estimate_s_rss). Values below 1e-5 are drawn in the shaded "s ≈ 0" band: a display floor, chosen
# when B53's old estimator could not resolve s below ~1e-5; the susieR port resolves to ~1e-8, so it is now a choice, not a limit.
S_RES = 1e-5
PANELS = {"eur_gctb": ("GCTB-10k", "#2a78d6"),
          "eur_aou": ("AoU-267k", "#eb6834"),
          "eur_polyfun": ("PolyFun-337k", "#1baf7a")}


def style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(INK_MUTED)
    ax.tick_params(colors=INK_MUTED, labelsize=8, labelcolor=INK)
    ax.grid(color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)


def label(ax, letter):
    ax.text(-0.14, 1.04 if letter not in "a" else 1.22, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", color=INK)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--b53-dir", default="../aou_ld/task_b53_out")
    ap.add_argument("--s-table", default="../aou_ld/task_b71_out/s_stats_susier.csv",
                    help="per-panel s (B71: real susieR estimate_s_rss); B53's old table is <b53-dir>/s_stats.csv with --s-col s")
    ap.add_argument("--s-col", default="s_susier_n")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    D = a.b53_dir
    S = pd.read_csv(a.s_table)
    W = S.pivot(index="locus", columns="panel", values=a.s_col).dropna()
    pip = {p: pd.read_csv(os.path.join(D, "export", f"pip_{p}.csv.gz"), usecols=["locus", "vkey", "pip"])
           for p in PANELS}

    fig = plt.figure(figsize=(11, 7.6), dpi=200)
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 1], hspace=0.42, wspace=1.1)
    ax_a, ax_b = fig.add_subplot(gs[0, 0:3]), fig.add_subplot(gs[0, 3:6])
    ax_pip = [fig.add_subplot(gs[1, 2 * i:2 * i + 2]) for i in range(3)]

    # a: ECDF of s, log x, direct-labelled with the median
    lo = S_RES / 2.5
    for ax_ in (ax_a,):
        ax_.axvspan(lo, S_RES, color=SHADE, zorder=0, lw=0)
    ax_a.text(np.sqrt(lo * S_RES), 0.5, "s ≈ 0", ha="center", va="center", rotation=90, fontsize=7.5, color=INK_MUTED)
    for p, (name, col) in PANELS.items():
        x = np.sort(np.clip(W[p].to_numpy(), np.sqrt(lo * S_RES), None))
        y = np.arange(1, len(x) + 1) / len(x)
        n0 = int((W[p] < S_RES).sum())
        ax_a.step(x, y, where="post", color=col, lw=2,
                  label=f"{name}: median {np.median(W[p]):.1e}, s ≈ 0 at {n0}")
    ax_a.set_xscale("log")
    ax_a.set_xlim(lo, W.values.max() * 1.3)
    ax_a.set_xlabel("s statistic (log scale)", fontsize=9, color=INK)
    ax_a.set_ylabel("Fraction of loci", fontsize=9, color=INK)
    ax_a.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="lower left", bbox_to_anchor=(0.0, 1.01),
                borderaxespad=0, handlelength=1.6)
    style(ax_a)
    label(ax_a, "a")
    ax_a.title.set_visible(False)

    # b: s per locus against GCTB-10k
    hi = W.values.max() * 1.3
    ax_b.axvspan(lo, S_RES, color=SHADE, zorder=0, lw=0)
    ax_b.axhspan(lo, S_RES, color=SHADE, zorder=0, lw=0)
    ax_b.plot([lo, hi], [lo, hi], color=INK_MUTED, lw=1, zorder=2)
    for p in ("eur_aou", "eur_polyfun"):
        name, col = PANELS[p]
        worse = int((W[p] > W["eur_gctb"]).sum())
        ax_b.scatter(W["eur_gctb"].clip(lower=np.sqrt(lo * S_RES)), W[p].clip(lower=np.sqrt(lo * S_RES)), s=14, color=col, alpha=0.6, linewidths=0.4, edgecolors="white", zorder=3,
                     label=f"{name}: s above GCTB at {worse}/{len(W)} loci")
    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xlim(lo, hi)
    ax_b.set_ylim(lo, hi)
    ax_b.set_xlabel("s, GCTB-10k (log)", fontsize=9, color=INK)
    ax_b.set_ylabel("s, other panel (log)", fontsize=9, color=INK)
    ax_b.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="upper left", bbox_to_anchor=(0.1, 1.0))
    ax_b.text(0.99, 0.01, "shaded: s < 1e-5 (≈ 0)", transform=ax_b.transAxes,
              ha="right", va="bottom", fontsize=7, color=INK_MUTED)
    style(ax_b)
    label(ax_b, "b")

    # c-e: pairwise PIP agreement
    for ax, (p, q), letter in zip(ax_pip, itertools.combinations(PANELS, 2), "cde"):
        j = pip[p].merge(pip[q], on=["locus", "vkey"], suffixes=("_p", "_q"))
        r = np.corrcoef(j["pip_p"], j["pip_q"])[0, 1]
        cross = int(((j["pip_p"] > 0.5) != (j["pip_q"] > 0.5)).sum())
        ax.plot([0, 1], [0, 1], color=INK_MUTED, lw=1, zorder=2)
        ax.scatter(j["pip_p"], j["pip_q"], s=5, color=INK, alpha=0.35, linewidths=0, zorder=3, rasterized=True)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect("equal")
        ax.set_xlabel(f"PIP, {PANELS[p][0]}", fontsize=9, color=INK)
        ax.set_ylabel(f"PIP, {PANELS[q][0]}", fontsize=9, color=INK)
        ax.text(0.03, 0.97, f"r = {r:.3f}\nn = {len(j):,}\ncross 0.5: {cross}", transform=ax.transAxes,
                va="top", fontsize=8, color=INK)
        style(ax)
        label(ax, letter)

    fig.suptitle(f"EUR reference panel: same z, same variants, three LD panels ({len(W)} loci)",
                 fontsize=11, color=INK, x=0.5, y=0.985)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out, bbox_inches="tight", facecolor="white")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
