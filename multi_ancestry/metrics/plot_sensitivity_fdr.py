#!/usr/bin/env python3
"""Sensitivity and false discovery against the PIP threshold, one panel each.

Ran asked for a sensitivity-versus-FDR curve. We built one and dropped it: FDR is a
ratio whose numerator and denominator both shrink with the threshold, so the parametric
curve is non-monotone and folds back on itself, and at these sample sizes most of its
apparent shape is a handful of variants changing state. Plotting each metric against
the threshold says the same thing without the false structure.

Vertical reference lines mark PIP 0.5 and 0.95, the two thresholds the manuscript
reports. An ordinal sweep -- every series sharing one base method, as in the rho sweep --
is drawn on a sequential light-to-dark ramp rather than categorical hues.

Usage:
  python3 plot_sensitivity_fdr.py --curves curves.csv --out fig.png [--title ...]
"""
import argparse
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Categorical slots 1-8 in fixed order (dataviz reference palette, light mode).
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
           "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SERIES_COLOR = {
    "MultiSuSiE":     PALETTE[0],
    "MESuSiE":        PALETTE[1],
    "SuSiEx-ours":    PALETTE[2],
    "SuSiE-merged":   PALETTE[3],
    "FINEMAP-merged": PALETTE[4],
}
# The order the arms appear in every figure and legend (Qin, 2026-09-18): the three
# ancestry-aware arms first, then the two merged-LD ones, so the grouping the paper argues
# about is the grouping the eye sees. Colour stays bound to the arm above, never to its
# position, so this ordering never repaints anything.
ARM_ORDER = ["SuSiEx-ours", "MultiSuSiE", "MESuSiE", "SuSiE-merged", "FINEMAP-merged"]


def arm_sort_key(name):
    """Sort key placing known arms in ARM_ORDER; unknown ones keep their own order after.

    Matches on the leading base name, so swept variants (`MultiSuSiE_L3`,
    `SuSiE-merged_mawson`) sort with their base arm and then among themselves.

    An `_L<n>` suffix sorts NUMERICALLY. Sorting it as text puts L=10 between L=1 and
    L=3, which silently scrambles any ordinal light-to-dark ramp built from this order.
    """
    m = re.search(r"_L(\d+)", name)
    within = (int(m.group(1)), name) if m else (0, name)
    for i, base in enumerate(ARM_ORDER):
        if name == base or name.startswith(base + "_"):
            return (i,) + within
    return (len(ARM_ORDER),) + within


def order_arms(names):
    """Unique arm names in canonical order."""
    return sorted(dict.fromkeys(names), key=arm_sort_key)
# Sequential ramp for an ordinal sweep: one hue, light to dark. Ordinal levels are ends
# of a scale, not unrelated series, so they must not get categorical hues.
SEQ_RAMP = ["#bcd7f2", "#7bb0e4", "#3f83cf", "#14508f"]
INK, INK_MUTED, GRID = "#1a1a19", "#6b6a63", "#e3e2dc"
REF_THRESHOLDS = (0.5, 0.95)


def ordinal_ramp(order):
    """Map an ordered series list onto the sequential ramp, or None if it doesn't fit."""
    if not (2 <= len(order) <= len(SEQ_RAMP)):
        return None
    return {m: SEQ_RAMP[len(SEQ_RAMP) - len(order) + i] for i, m in enumerate(order)}


def colour_for(method, order):
    """Fixed slot for a primary arm; otherwise the next palette slot in appearance order.
    Never cycles -- a 9th series needs small multiples, not a reused hue."""
    if method in SERIES_COLOR:
        return SERIES_COLOR[method]
    i = order.index(method)
    if i >= len(PALETTE):
        raise SystemExit(f"[error] more than {len(PALETTE)} series; facet instead of cycling hues")
    return PALETTE[i]


def plot_vs_threshold(d, out, title):
    """Companion: each metric against the threshold directly.

    The sensitivity-vs-FDR curve is what Ran asked for, but FDR is a ratio whose
    numerator and denominator both shrink with the threshold, so it is not monotone
    and the parametric curve saws back and forth. Plotting each metric against the
    threshold separately is monotone, uncluttered, and lets a reader pick an
    operating point -- read the two together.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.6), dpi=200, sharex=True)
    for ax, col, label in ((axes[0], "sensitivity", "Sensitivity  →  better"),
                           (axes[1], "fdr", "False discovery rate  →  worse")):
        order = order_arms(d["method"])
        ramp = (ordinal_ramp(order)
                if len({m.split("_")[0] for m in order}) == 1 else None)
        for method in (order if ramp else order_arms(d["method"])):
            g = d[d["method"] == method].sort_values("threshold")
            per_anc = g["arm"].iloc[0] == "per-ancestry"
            ax.plot(g["threshold"], g[col], "-" if per_anc else "--", lw=2,
                    color=(ramp[method] if ramp else colour_for(method, order)), zorder=3,
                    label=f"{method} ({'per-ancestry' if per_anc else 'merged LD'})")
        ax.set_xlabel("PIP threshold", fontsize=10, color=INK)
        ax.set_ylabel(label, fontsize=10, color=INK)
        for rt in REF_THRESHOLDS:
            ax.axvline(rt, color=GRID, lw=1, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=INK_MUTED, labelsize=9)
        ax.set_ylim(-0.02, 1.02)
    axes[0].legend(frameon=False, fontsize=8, loc="lower left", labelcolor=INK)
    fig.suptitle(title, fontsize=12.5, color=INK, x=0.008, ha="left", y=1.02)
    note = ""
    if "benchmark" in d.columns:
        r = d.iloc[0]
        bt = r["benchmark_threshold"]
        mode = ("both sides gated at the same threshold (Ran's matched sweep)"
                if pd.isna(bt) else f"benchmark positives fixed at PIP>{bt:g}")
        note = (f"Ground truth: {r['benchmark']} only — {int(r['benchmark_loci'])} loci, "
                f"{mode}, noise below {r['noise_threshold']:g}; "
                f"{r.get('matching', 'exact variant match')}.")
    fig.text(0.01, -0.03, note, fontsize=8, color=INK_MUTED, va="top")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    print(f"[done] wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--curves", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="Sensitivity and false discovery by PIP threshold")
    args = ap.parse_args()
    plot_vs_threshold(pd.read_csv(args.curves), args.out, args.title)


if __name__ == "__main__":
    main()
