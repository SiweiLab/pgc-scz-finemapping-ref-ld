#!/usr/bin/env python3
"""Variant counts by stage, for each of the three input pipelines.

2026-09-28: a picture of `variant_domain_counts.md`, to track where variants are lost
and how much of each stage is palindromic.

Form: horizontal stacked bars, one row per stage, one panel per pipeline, length = variants.
Each bar carries its non-palindromic count in the solid segment and its palindromic count
and share in the pale one, so the two add to the total printed at the right. There is no
separate "after ambiguous-free rerun" row: that number IS the non-palindromic count
(228,205 - 34,620 = 193,585, and so on), so a fourth bar per pipeline drew the same
quantity twice.

THE STAGE PANELS ARE GREY ON PURPOSE. Their hues would otherwise be PALETTE[0..2],
which already mean MultiSuSiE, MESuSiE and SuSiEx in every other figure -- an input
pipeline tinted like an arm invites exactly the wrong reading. Only the per-arm panel uses
the arm palette, where those hues mean what they always mean.
Horizontal because the stage names are long; stacked because the palindromic share IS the
story -- it is 15.1-15.4% at every stage in every pipeline, which is why the no-palindromic rerun
moves every arm the same way.

Every stage is now measured. A cluster-side count read the two EAS stages the laptop could not (no
allele columns here): the panel is 33,051 of 218,220 (15.15%) from the panel's own
hail_A1/hail_A2, and harmonised is 23,650 of 153,980 (15.36%). It also confirmed, by reading
the panel directly, the 34,620 the laptop had only inferred for EUR from the SuSiEx export.

EAS harmonised (153,980) and MESuSiE's fitted domain (153,978) are drawn as ONE row. That count
showed they are not the same set -- v1 joins the EAS panel on MarkerName and v2-union on
SNP_ID, so two variants with a blank MarkerName survive only in v2-union -- but two variants
in 153,980 is 0.001%, invisible at this width, and two bars would imply a difference the
reader can see. The caption carries the distinction instead.

Yuan et al.'s domain is a dashed rule rather than a bar: it is a different pipeline's endpoint,
not a stage of ours, and drawing it as a bar in our stack would imply a lineage it does not
have.

Usage:
  python3 plot_variant_domains.py --out figures/supp_variant_domains.png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                     "axes.titleweight": "bold"})
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_sensitivity_fdr import INK, INK_MUTED, GRID, PALETTE
from plot_cs_split import ARM_HUE, EUR_ONLY, METHOD_ORDER, disp

# (stage label, total, palindromic). `None` would draw a flat grey "not measurable" bar;
# now nothing needs it, but the branch is kept for a stage added before its alleles
# are in hand -- guessing a share is worse than saying it is unmeasured.
STAGE = INK_MUTED          # one neutral for every stage bar; hue is reserved for arms
PIPELINES = [
    ("EUR", STAGE, [
        ("GWAS sumstats", 231_713, 35_188),
        ("reference LD panel", 228_205, 34_620),
        ("fitted — SuSiEx, EUR arms", 228_205, 34_620),
    ]),
    ("EAS", STAGE, [
        ("GWAS sumstats", 202_739, 31_303),
        ("reference LD panel", 218_220, 33_051),
        ("harmonised = fitted, MESuSiE", 153_980, 23_650),
    ]),
    ("trans-ancestry meta", STAGE, [
        ("core meta sumstats", 228_754, 34_776),
        ("fitted — 4 merged arms", 215_341, 32_641),
    ]),
]

# Panel 4: what each arm actually fitted, with the LD it was fitted ON -- the count alone
# does not say that, and two arms can share a count while using different matrices.
#
# "merged" is not a fourth panel above because it is not a fourth source: merge_eur_eas_ld.py
# builds it from the SAME two per-ancestry R matrices, elementwise,
#     R_merged = (w_EUR * R_EUR + w_EAS * R_EAS) / (w_EUR + w_EAS),  diag forced to 1,
# with w = each ancestry's median NEFF from its own sumstats (LDmergeFM's METAL weighting,
# computed directly rather than through its file format). Its domain is therefore the two
# panels intersected, then intersected again with the core meta sumstats.
#
# Panel 4: what each arm actually fitted. The three panels above are organised by INPUT
# pipeline, and MultiSuSiE has no row in any of them -- its domain is the union of both
# ancestries' valid sets restricted to LD present in EVERY ancestry
# (harmonize_loci_v2_union.py:203-208), a fourth lineage. We could not find its number,
# which is the point of this panel: every arm, its own count, in the per-arm colours the
# rest of the figures use.
ARMS = [
    ("SuSiEx-ours",    228_205, 34_620, "EUR ∪ EAS own-valid sets, each LD kept separate"),
    ("MultiSuSiE",     217_154, 32_892, "EUR ∪ EAS, but LD required in BOTH"),
    ("MESuSiE",        153_978, 23_650, "EUR ∩ EAS"),
    ("FINEMAP-merged", 215_341, 32_641, "ESS-weighted mean of EUR & EAS R"),
    ("SuSiE-merged",   215_341, 32_641, "ESS-weighted mean of EUR & EAS R"),
    ("RSparsePro",     215_341, 32_641, "ESS-weighted mean of EUR & EAS R"),
    ("CARMA",          215_341, 32_641, "ESS-weighted mean of EUR & EAS R"),
    ("EUR-FINEMAP",    228_205, 34_620, "EUR panel only"),
    ("EUR-SuSiE-RSS",  228_205, 34_620, "EUR panel only"),
]
# "et al." in italic, bold to match every other label: mathtext's italic slot is pointed at DejaVu Sans bold italic
plt.rcParams.update({"mathtext.fontset": "custom", "mathtext.it": "DejaVu Sans:italic:bold"})
YUAN = r"Benchmark (Yuan $\mathit{et\ al.}\!$)"
YUAN_N = 231_272          # his union across the 250 loci, zero palindromic

# Pipeline v2 (2026-10-01): every count from the v2 run. Sumstats and GCTB: master flags (the master table's source flags).
# Arm domains and the per-ancestry valid sets: domains logged by harmonize.py, palindromic = with-pal minus no-pal domain (exact).
# The AoU EAS panel's palindromic count was counted on the cluster from the stored vkeys (40,554 of 265,899).
V2_PIPELINES = [
    ("EUR", STAGE, [
        ("GWAS sumstats", 231_713, 35_188),
        ("reference LD panel (GCTB-10k)", 235_941, 35_838),
        ("sumstats ∩ LD (SuSiEx EUR, EUR arms)", 228_853, 34_718),
    ]),
    ("EAS", STAGE, [
        ("GWAS sumstats", 202_739, 31_303),
        ("reference LD panel (AoU EAS)", 265_899, 40_554),   # from the panel's stored vkeys
        ("sumstats ∩ LD (SuSiEx EAS)", 198_611, 30_702),
    ]),
    ("trans-ancestry meta", STAGE, [
        ("core meta sumstats", 228_754, 34_776),
        ("fitted — 4 merged arms", 221_830, 33_705),
    ]),
]
V2_ARMS = [
    ("SuSiEx-ours",    269_466, 41_069, "EUR ∪ EAS own-valid sets, each LD kept separate"),
    ("MultiSuSiE",     231_220, 35_143, "EUR ∪ EAS, LD required in BOTH; missing z censored"),
    ("MESuSiE",        157_998, 24_351, "EUR ∩ EAS"),
    ("FINEMAP-merged", 221_830, 33_705, "ESS-weighted mean of EUR & EAS R"),
    ("SuSiE-merged",   221_830, 33_705, "ESS-weighted mean of EUR & EAS R"),
    ("RSparsePro",     221_830, 33_705, "ESS-weighted mean of EUR & EAS R"),
    ("CARMA",          221_830, 33_705, "ESS-weighted mean of EUR & EAS R"),
    ("EUR-FINEMAP",    228_853, 34_718, "EUR panel only"),
    ("EUR-SuSiE-RSS",  228_853, 34_718, "EUR panel only"),
]
# The palindromic count and its percentage sit inside the pale segment on one line, as
# "23,650 (15.4%)" -- see _pal_label for why that needs regular weight at 5.4pt and a
# wider figure to clear the narrowest segment.
#
# The totals are at a FIXED x rather than offset from each bar's end: the bars span
# 154k-232k, so anything placed relative to a bar end staggers down the panel.
TOTAL_X, XMAX = 238_000, 272_000
# xlim runs past the data to hold those columns, so the axis is trimmed to where the bars
# stop -- an axis drawn under the label columns reads as plot area a bar could occupy.
AXIS_END = 232_000


def _pal_label(ax, x, y, pal, tot):
    """Palindromic count with its percentage in parentheses, on one line.

    Set in regular weight, unlike every other number here. That is partly what was requested
    for and partly what makes it fit: at 5.4pt regular the string is ~0.60 in against
    ~0.69 in of room in the NARROWEST pale segment (MESuSiE / EAS harmonised, 23,650
    variants). Bold, or a point larger, overflows there -- which is also why the figure
    widened from 10.4 to 11.6 in. The percentage keeps one decimal: 15.1 / 15.2 / 15.4 is
    the distinction the figure exists to show, and at whole percent every row reads 15%.
    """
    ax.annotate(f"{pal:,} ({pal / tot * 100:.1f}%)", (x, y), va="center", ha="center",
                fontsize=5.4, color=INK, fontweight="normal")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--v2", action="store_true", help="pipeline v2 counts (2026-10-01) instead of v1's")
    ap.add_argument("--arms-only", action="store_true",
                    help="only the 'what each arm fitted' block, trans-ancestry arms, sized for Fig 1A beside the flowchart")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    global PIPELINES, ARMS, TOTAL_X, XMAX, AXIS_END
    if args.v2:
        PIPELINES, ARMS = V2_PIPELINES, V2_ARMS
        TOTAL_X, XMAX, AXIS_END = 276_000, 312_000, 270_000

    if args.arms_only:
        return arms_only(args)
    n = [len(p[2]) for p in PIPELINES] + [len(ARMS) + 0.8]
    fig, axes = plt.subplots(len(PIPELINES) + 1, 1, figsize=(11.6, 11.2), dpi=300,
                             gridspec_kw={"height_ratios": n}, sharex=True)
    fig.subplots_adjust(hspace=0.40, right=0.80, top=0.955, bottom=0.07)

    for ax, (name, hue, stages) in zip(axes, PIPELINES):
        y = list(range(len(stages)))[::-1]
        for yi, (lab, tot, pal) in zip(y, stages):
            if pal is None:
                ax.barh(yi, tot, height=0.58, color=GRID, edgecolor="white", linewidth=1.4)
                ax.annotate("palindromic share not measured (per-variant list on the cluster)",
                            (tot * 0.012, yi), va="center", ha="left",
                            fontsize=7.6, color=INK_MUTED, style="italic")
            else:
                keep = tot - pal
                ax.barh(yi, keep, height=0.58, color=hue, edgecolor="white", linewidth=1.4)
                ax.annotate(f"{keep:,}", (keep - 3000, yi), va="center", ha="right",
                            fontsize=6.8, color="white")
                if pal:
                    ax.barh(yi, pal, left=keep, height=0.58, color=hue, alpha=0.30,
                            edgecolor="white", linewidth=1.4)
                    _pal_label(ax, keep + pal / 2, yi, pal, tot)
            ax.annotate(f"{tot:,}", (TOTAL_X, yi), va="center", ha="left",
                        fontsize=9.2, color=INK)
        ax.axvline(YUAN_N, color=INK_MUTED, lw=1.1, ls="--", zorder=1)
        ax.set_yticks(y)
        ax.set_yticklabels([s[0] for s in stages], fontsize=9.4, color=INK)
        ax.set_title(name, fontsize=11.5, color=INK, loc="left", pad=5)
        ax.set_xlim(0, XMAX)
        ax.grid(False)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color(GRID)
        ax.spines["bottom"].set_bounds(0, AXIS_END)
        ax.set_xticks(range(0, (250_001 if args.v2 else 200_001), 50_000))
        ax.tick_params(colors=INK_MUTED, labelsize=9.4)
    # in axes coords vertically, so it sits clear above the top bar instead of on it
    axes[0].annotate(YUAN + f"  {YUAN_N:,}", (YUAN_N, 1.04), xycoords=("data", "axes fraction"),
                     ha="right", va="bottom", fontsize=8.8, color=INK_MUTED)
    # ---- panel 4: per arm -------------------------------------------------------------
    axA = axes[-1]
    ya = list(range(len(ARMS)))[::-1]
    for yi, (arm, tot, pal, lineage) in zip(ya, ARMS):
        hue = ARM_HUE[arm]
        axA.barh(yi, tot - pal, height=0.62, color=hue, edgecolor="white", linewidth=1.4)
        axA.barh(yi, pal, left=tot - pal, height=0.62, color=hue, alpha=0.30,
                 edgecolor="white", linewidth=1.4)
        axA.annotate(f"{tot - pal:,}", (tot - pal - 3000, yi), va="center", ha="right",
                     fontsize=6.8, color="white")
        _pal_label(axA, tot - pal / 2, yi, pal, tot)
        axA.annotate(f"{tot:,}", (TOTAL_X, yi), va="center", ha="left",
                     fontsize=9.0, color=INK)
        axA.annotate(lineage, (2600, yi), va="center", ha="left", fontsize=6.8,
                     color="white" if tot - pal > 60_000 else INK_MUTED)
    axA.axvline(YUAN_N, color=INK_MUTED, lw=1.1, ls="--", zorder=1)
    axA.set_yticks(ya)
    axA.set_yticklabels([disp(a[0]) for a in ARMS], fontsize=9.4, color=INK)
    axA.set_title("What each method fitted", fontsize=11.5, color=INK, loc="left", pad=5)
    axA.set_xlim(0, XMAX)
    axA.grid(False)
    for sp in ("top", "right", "left"):
        axA.spines[sp].set_visible(False)
    axA.spines["bottom"].set_color(GRID)
    axA.spines["bottom"].set_bounds(0, AXIS_END)
    axA.set_xticks(range(0, (250_001 if args.v2 else 200_001), 50_000))
    axA.tick_params(colors=INK_MUTED, labelsize=9.4)

    axes[-1].set_xlabel("Variants across the 250 loci", fontsize=11.5, color=INK)
    axes[-1].xaxis.set_major_formatter(lambda v, _: f"{v/1000:,.0f}k")

    h = [mpatches.Patch(facecolor=INK_MUTED, label="not palindromic"),
         mpatches.Patch(facecolor=INK_MUTED, alpha=0.30, label="palindromic (A/T, C/G)"),
         mlines.Line2D([], [], color=INK_MUTED, lw=1.1, ls="--",
                       label="Benchmark (non-palindromic only)")]
    lg = fig.legend(handles=h, frameon=False, fontsize=9.4, labelcolor=INK,
                    loc="center left", bbox_to_anchor=(0.815, 0.5), alignment="left")
    fig.savefig(args.out, bbox_inches="tight", facecolor="white",
                bbox_extra_artists=[lg])
    print(f"[done] wrote {args.out}")


def arms_only(args):
    """Fig 1A inset: the seven trans-ancestry arms' fitted domains (EUR-only arms are out of the primary analysis).
    The sumstats and LD-panel counts are in the hand-drawn flowchart, so only the per-arm block is drawn here."""
    arms = [r for r in ARMS if r[0] not in EUR_ONLY]
    fig, ax = plt.subplots(figsize=(7.4, 3.7), dpi=300)
    ya = list(range(len(arms)))[::-1]
    for yi, (arm, tot, pal, lineage) in zip(ya, arms):
        hue = ARM_HUE[arm]
        ax.barh(yi, tot - pal, height=0.66, color=hue, edgecolor="white", linewidth=1.4)
        ax.barh(yi, pal, left=tot - pal, height=0.66, color=hue, alpha=0.30, edgecolor="white", linewidth=1.4)
        ax.annotate(f"{tot - pal:,}", (tot - pal - 3000, yi), va="center", ha="right", fontsize=7.6, color="white")
        if pal < 30_000:   # too narrow for its label at this width (MESuSiE): set it just past the bar
            ax.annotate(f"{pal:,} ({pal / tot * 100:.1f}%)", (tot + 1500, yi), va="center", ha="left",
                        fontsize=5.4, color=INK, fontweight="normal")
        else:
            _pal_label(ax, tot - pal / 2, yi, pal, tot)
        ax.annotate(f"{tot:,}", (TOTAL_X, yi), va="center", ha="left", fontsize=9.0, color=INK)
        ax.annotate(lineage, (2600, yi), va="center", ha="left", fontsize=7.2, color="white")
    ax.axvline(YUAN_N, color=INK_MUTED, lw=1.1, ls="--", zorder=1)
    # starts at the line and runs right
    ax.annotate(YUAN + f"  {YUAN_N:,}", (YUAN_N, 1.0), xycoords=("data", "axes fraction"), xytext=(3, 0),
                textcoords="offset points", ha="left", va="bottom", fontsize=8.4, color=INK_MUTED)
    ax.set_yticks(ya)
    ax.set_yticklabels([disp(a[0]) for a in arms], fontsize=9.6, color=INK)
    ax.set_title("Variants each method fitted (250 loci)", fontsize=11.5, color=INK, loc="left", pad=14)
    ax.set_xlim(0, XMAX)
    ax.grid(False)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.spines["bottom"].set_bounds(0, AXIS_END)
    ax.set_xticks(range(0, 250_001, 50_000))
    ax.xaxis.set_major_formatter(lambda v, _: f"{v/1000:,.0f}k")
    ax.tick_params(colors=INK_MUTED, labelsize=9.0)
    h = [mpatches.Patch(facecolor=INK_MUTED, label="not palindromic"),
         mpatches.Patch(facecolor=INK_MUTED, alpha=0.30, label="palindromic (A/T, C/G)"),
         mlines.Line2D([], [], color=INK_MUTED, lw=1.1, ls="--", label="Benchmark (non-palindromic only)")]
    ax.legend(handles=h, frameon=False, fontsize=8.6, labelcolor=INK, ncol=3, loc="upper center",
              bbox_to_anchor=(0.45, -0.13))
    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    print(f"[done] wrote {args.out}")


if __name__ == "__main__":
    main()
