#!/usr/bin/env python3
"""Fig. 2D: the 13 highest-confidence benchmark variants against all nine arms.

Rows are the 13 variants with benchmark PIP > 0.9 — one per locus, and the same |P_B| = 13
as the alpha = 0.9 column of the credible-set-split figure, so the two share a denominator.

ONE THING IS COLOUR: our in-CS PIP at the benchmark's lead variant (0 outside every CS), on a single light-to-dark
ramp. Everything else that used to be encoded in the fill is a mark beside the number, so
the reader decodes one visual channel instead of three:

  *   our top IN-CS variant is NOT the benchmark's lead (out-of-CS PIPs are not calls)
  †   the lead sits OUTSIDE every credible set the arm formed
  §   the arm did NOT resolve exactly one credible set (none, or two or more)

§ covers both non-commitments because the metric treats them alike. For the SuSiE family
"none" is usually the PURITY GATE, not the model finding nothing: at chr4:698,290-828,700
all four SuSiE-family arms hold rs35734242 at PIP 0.949-0.972 and lose the set to the gate,
while FINEMAP-merged and CARMA keep one only because they have no native gate -- theirs sit
at common_purity 0.48 and 0.44, below the 0.5 the others are held to.

The marks are independent: an arm can hold the lead at high PIP inside a credible set and
still not have it on top. FINEMAP routinely reports two or more credible sets while leading
on the right variant, which is reporting granularity rather than error.

Rows are ordered by CHROMOSOME and position. The first column is the benchmark's own PIP on
the same colour scale. EUR-only arms are excluded: they are out of the primary analysis.
count the runtime figure brackets on).

Usage:
  python3 plot_fig2d.py --per-locus figures/fig1c_193_C_per_locus.csv \
      --out figures/fig2d_top13.png
"""
import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                     "axes.titleweight": "bold"})
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerPatch
from paths import V2
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_sensfdr import BENCH, bench_cs_loci, KAI
from plot_cs_split import RAN_ORDER, EUR_ONLY, disp
from plot_sensitivity_fdr import INK, INK_MUTED, GRID

ABSENT = "#D8D8D8"
# Single-hue sequential ramp (the validated reference palette's blue, steps 100 -> 700, from a near-white start so a
# PIP near 0 still reads as a value, distinct from the grey "absent" cell). Lightness alone carries the value, so it is
# colour-blind safe; replaced matplotlib's Purples.
from matplotlib.colors import LinearSegmentedColormap
RAMP = LinearSegmentedColormap.from_list("blue_seq", ["#f4f8fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6",
                                                      "#1c5cab", "#104281", "#0d366b"])
MARK_NOT_LEAD, MARK_OUT_CS, MARK_NO_SINGLE_CS = "*", "†", "§"


def chrom_key(locus):
    m = re.match(r"chr(\w+)_(\d+)_", locus)
    c = m.group(1)
    return (int(c) if c.isdigit() else 99, int(m.group(2)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--export", default=os.path.join(here, "..", "aou_ld", "task_v2_out", "noambig"))
    ap.add_argument("--per-locus",
                    default=os.path.join(V2, "fig1c_193_C_per_locus.csv"))
    ap.add_argument("--cell", type=float, default=0.40, help="inches per cell")
    ap.add_argument("--suffix", default="-noambig", help="export label suffix, e.g. --suffix=-noambig for task_v2_out/noambig")
    ap.add_argument("--transpose", action="store_true",
                    help="methods as rows, loci as columns (benchmark the top row), for a full-width row under Fig. 2a-c")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    keep, bmem = bench_cs_loci(KAI)
    bench_ncs = bmem.groupby("locus").cs.nunique().to_dict()
    b = pd.read_csv(BENCH)
    b = b[b.locus.isin(keep)]
    top = b[b.pip > 0.9].copy()
    top["k"] = top.locus.map(chrom_key)
    top = top.sort_values("k")                       # by chromosome, not by PIP
    loci = list(top.locus)
    lead = dict(zip(top.locus, top.rsid))
    lead_bp = dict(zip(top.locus, top.bp.astype(int)))
    bench_pip = dict(zip(top.locus, top.pip.astype(float)))

    P = pd.read_csv(args.per_locus)
    # EUR-only arms are out of the primary analysis
    arms = [a for a in RAN_ORDER if a in set(P.arm) and a not in EUR_ONLY]
    ncs = P.set_index(["arm", "locus"]).n_cs.to_dict()

    size = {}
    pv = pd.read_csv(os.path.join(args.export, f"primary_MultiSuSiE{args.suffix}_pip.csv.gz"),
                     low_memory=False)
    size = pv[pv.locus.isin(loci)].groupby("locus").size().to_dict()

    pip, seen, incs, armtop = {}, {}, {}, {}
    for a in arms:
        d = pd.read_csv(os.path.join(args.export, f"primary_{a}{args.suffix}_cs_purity.csv.gz"),
                        low_memory=False)
        d = d[d.locus.isin(loci)]
        pip[a] = {(l, int(bp)): float(p)
                  for l, bp, p in zip(d.locus, d.bp.astype(int), d.pip.astype(float))}
        seen[a] = {l: set(g.bp.astype(int)) for l, g in d.groupby("locus")}
        m = d[d.cs_ids_all.notna()]
        incs[a] = {(l, int(bp)) for l, bp in zip(m.locus, m.bp.astype(int))}
        # the arm's own top variant per locus, among its CREDIBLE-SET members only: a
        # high PIP the arm never put in a set is not a call. Ties with the lead count as the lead
        armtop[a] = m.groupby("locus").pip.max().to_dict()

    # the benchmark is column 0, so its PIP is read on the same colour scale as ours
    # instead of sitting in the row label as a bracketed number
    # "1.00" printed as "1"; everything else two decimals, no leading zero
    def fmt(v):
        t = f"{v:.2f}"
        return "1" if t == "1.00" else t.lstrip("0")

    cols = ["SuSiEx-benchmark"] + arms
    M = np.full((len(loci), len(cols)), np.nan)
    txt = [["" for _ in cols] for _ in loci]
    for i, l in enumerate(loci):
        M[i, 0] = bench_pip[l]
        # the benchmark gets § on the same rule as our arms (chr18: three sets)
        txt[i][0] = fmt(bench_pip[l]) + ("" if bench_ncs.get(l, 0) == 1 else MARK_NO_SINGLE_CS)
        for j, a in enumerate(arms, start=1):
            bp = lead_bp[l]
            if bp not in seen[a].get(l, set()):
                txt[i][j] = "–"
                continue
            # in-CS PIP: a lead the arm left outside every set scores 0 (and carries the dagger)
            v = pip[a].get((l, bp), 0.0) if (l, bp) in incs[a] else 0.0
            M[i, j] = v
            k = int(ncs.get((a, l), 0))
            marks = ""
            # read off the PIPs, not Fig 1C's `exact`: that is only scored at a single CS, so
            # it flagged every multi-CS arm here even when its top variant was the lead
            if v < armtop[a].get(l, 0.0):
                marks += MARK_NOT_LEAD
            if k >= 1 and (l, bp) not in incs[a]:
                marks += MARK_OUT_CS
            if k != 1:
                marks += MARK_NO_SINGLE_CS
            txt[i][j] = fmt(v) + marks

    rowlab = ["Benchmark (in-sample)"] + [disp(a) for a in arms]
    loclab = []
    for l in loci:
        c, s0, s1 = l.split("_")
        loclab.append(f"{c}:{int(s0)/1e6:.1f}–{int(s1)/1e6:.1f} Mb   {lead[l]}")
    txt0 = txt                      # locus-major, for the printout below
    if args.transpose:   # methods as rows (benchmark on top), loci as columns
        M = M.T
        txt = [list(r) for r in zip(*txt)]
    nr, nc = M.shape
    if args.transpose:
        fig, ax = plt.subplots(figsize=(args.cell * 1.55 * nc + 3.6, args.cell * nr + 3.2), dpi=300)
        aspect = "auto"
    else:
        fig, ax = plt.subplots(figsize=(args.cell * len(arms) + 6.0, args.cell * len(loci) + 3.4), dpi=300)
        aspect = "equal"
    ax.imshow(np.ma.masked_invalid(M), cmap=RAMP, vmin=0, vmax=1, aspect=aspect, interpolation="nearest")
    for i in range(nr):
        for j in range(nc):
            # absent (NaN) is masked in imshow, so it shows as a blank cell of exactly the
            # same size; no patch, whose edge stroke would spill into the neighbours
            ax.annotate(txt[i][j], (j, i), ha="center", va="center", fontsize=8.5,
                        color="white" if (not np.isnan(M[i, j]) and M[i, j] > 0.62)
                        else INK, zorder=4)
    if args.transpose:
        ax.axhline(0.5, color=INK, lw=2.2, zorder=7)      # benchmark row | our methods
        ax.set_xticks(range(nc))
        # start position + lead only: the full span does not fit at an angle without colliding
        short = [f"{l.split('_')[0]}:{int(l.split('_')[1]) / 1e6:.1f} Mb\n{lead[l]}" for l in loci]
        ax.set_xticklabels(short, fontsize=9.6, color=INK, rotation=40, ha="left")
        ax.set_yticks(range(nr))
        ax.set_yticklabels(rowlab, fontsize=11, color=INK)
    else:
        ax.axvline(0.5, color=INK, lw=2.2, zorder=7)      # benchmark | our arms
        ax.set_xticks(range(nc))
        ax.set_xticklabels(rowlab, fontsize=12, color=INK, rotation=32, ha="left")
        ax.set_yticks(range(nr))
        ax.set_yticklabels(loclab, fontsize=10.5, color=INK)
    ax.xaxis.set_ticks_position("top")
    ax.tick_params(colors=INK_MUTED, length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)

    sm = plt.cm.ScalarMappable(cmap=RAMP, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, ax=ax, fraction=0.020, pad=0.02)
    cb.set_label("PIP at the benchmark top variant\n(0 if outside a credible set)", fontsize=11, color=INK)
    cb.ax.tick_params(colors=INK_MUTED, labelsize=10)
    cb.outline.set_visible(False)

    class _Key(HandlerPatch):
        """Draw the key in the handle column -- a mark symbol, or the dash for an absent lead -- so every
        label starts at the same x."""
        def __init__(self, sym, cell=False):
            super().__init__()
            self.sym, self.cell = sym, cell

        def create_artists(self, legend, orig, xd, yd, w, h, fs, trans):
            side = min(w, h) * 1.25
            x0, y0 = xd + (w - side) / 2, yd + (h - side) / 2
            out = []
            if self.cell:
                out.append(mpatches.Rectangle((x0, y0), side, side, facecolor=ABSENT, edgecolor="#9A9A9A", lw=0.8,
                                              transform=trans))
            out.append(plt.Text(x0 + side / 2, y0 + side / 2, self.sym, ha="center", va="center", fontsize=fs,
                                fontweight="bold", color=INK, transform=trans))
            return out

    keys = [(MARK_NOT_LEAD, "top credible-set variant is a different variant", False),
            (MARK_OUT_CS, "top variant outside every credible set", False),
            (MARK_NO_SINGLE_CS, "no single credible set (none, or two or more)", False),
            ("–", "top variant absent from the method's input", False)]
    h = [mpatches.Patch(label=lab) for _, lab, _ in keys]
    ax.legend(handles=h, handler_map={hh: _Key(sym, cell) for hh, (sym, _, cell) in zip(h, keys)}, frameon=False,
              fontsize=11, labelcolor=INK, ncol=1, loc="upper left",
              # transposed: one key per line to the right of the colourbar; otherwise under the grid
              bbox_to_anchor=(1.115, 1.0) if args.transpose else (0.0, -0.012), handlelength=1.0, handletextpad=0.5)

    note = (
        "Rows are the 13 variants with benchmark PIP > 0.9, one per locus — the same "
        "|P_B| = 13 as the alpha = 0.9 column of the credible-set-split figure. Ordered by "
        "CHROMOSOME and position, not by PIP. The FIRST COLUMN is the benchmark's own PIP, with § on the same rule as ours, on "
        "the same colour scale, so our columns are read against it. EUR-only arms are "
        "excluded: they are out of the primary analysis. "
        "COLOUR IS ONE THING ONLY: our in-credible-set PIP at the benchmark's lead variant, "
        "not at whatever we picked; a lead the arm left outside every set scores 0. Variants "
        "are the harmonised master without palindromic (A/T, C/G) variants. Everything else is a mark beside "
        "the number: * means our top credible-set member is not the lead (a higher PIP the "
        "arm never put in a set does not count), † that the lead sits outside every credible "
        "set the arm formed, "
        "and § that it did not resolve exactly one credible set — none, or two or "
        "more, which the metric treats alike. For the SuSiE family 'none' is usually the PURITY GATE rather than the model finding nothing (a set below min_abs_corr 0.5 is dropped); FINEMAP-merged and CARMA have no native gate. An "
        "arm can hold the lead at high PIP inside a credible set and still not have it on "
        "top; FINEMAP routinely reports two or more credible sets while still leading on the "
        "right variant. An en dash means the lead is absent from that arm's input entirely.")
    cp = args.out.rsplit(".", 1)[0] + "_caption.txt"
    with open(cp, "w") as fh:
        fh.write(" ".join(note.split()) + "\n")
    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    print(f"[done] wrote {args.out}\n")
    for i, l in enumerate(loci):
        print(f"{l:28s} {lead[l]:12s} n={size.get(l,0):5,d}  " +
              " ".join(f"{txt0[i][j]:>8s}" for j in range(len(arms))))


if __name__ == "__main__":
    main()
