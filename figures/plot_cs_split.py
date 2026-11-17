#!/usr/bin/env python3
"""FDR and recall split by whether the CALL sits inside the arm's own credible set.

Motivation, from the 2026-09-22 diagnosis. The published benchmark is credible-set
restricted *by construction*: across the 193 loci that formed a credible set, no variant
outside one exceeds PIP 0.0077, and the 57 loci that formed none are all exactly 0. Our own
SuSiEx-ours reproduces that cap (0.0075), so it is a property of SuSiEx's output format,
not a quirk of Yuan et al.'s run -- posterior mass held by a purity-filtered component is discarded
rather than redistributed (`KNOWN_ISSUE_susiex_cs_restricted_pips.md`).

MultiSuSiE, MESuSiE and SuSiE-merged export COMPLETE posteriors, with up to PIP 0.97
outside any credible set. Scored against a CS-restricted benchmark, every one of those
calls is a false positive almost automatically. That is a format mismatch, not a
disagreement about biology, and it falls on three of the five arms drawn here.

Rather than gate it away -- which deletes calls and silently changes the metric's
definition -- this figure SEPARATES the two populations so the size of the effect is
visible:

  row 1  FDR,    paired bars: calls inside the arm's own CS vs calls outside any of them
  row 2  recall, paired bars: same split. The two groups partition the true positives, so
         they SUM to the arm's total recall, printed above each pair.

Read the OUTSIDE bars with the asymmetry in mind. An out-of-CS call is not *structurally
barred* from being a true positive -- it can land on a variant inside the BENCHMARK's
credible set, and some do -- but it is heavily disadvantaged, so a high FDR there is partly
guaranteed and is not independent evidence that those calls are wrong.

**This is the PRIMARY run, not the L sweep.** Credible-set membership was only ever
exported for `results/benchmark/pip_vectors_full250_qcexcl.csv.gz` (whose `pip` column
is byte-identical to it). The `_bars` figures are built from `pip_vectors_L_sweep_v2` at
L=5 and give materially different numbers for the same arm and threshold. Do not quote the
two side by side.

Pinning is matched (s = alpha) throughout, so each column asks a symmetric question and
recall is on a 0-1 scale.

Usage:
  python3 plot_cs_split.py --out figures/primary_benchfull_matched_193_csplit.png
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# every label bold: ticks and legends inherit font.weight
plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                     "axes.titleweight": "bold"})
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, ".")
from bench_sensfdr import bench_cs_loci, BENCH, BENCH_RAW
from sensitivity_fdr import load_pips, variant_key
from plot_sensitivity_fdr import (INK, INK_MUTED, GRID, SERIES_COLOR, order_arms,
                                  PALETTE as PALETTE_)

# arm -> file holding its credible-set membership. LEGACY PATH ONLY: the first export covered
# four arms, and None here meant the arm had no membership and was drawn as one
# undifferentiated bar. That is no longer the state of the figure. The 2026-09 primary
# export supplies credible-set membership for all NINE arms, FINEMAP included, and under
# --export this dict is overwritten per arm from what the export actually contains. So
# nothing is drawn unsplit any more; the legend is built from what was drawn rather than
# from a fixed list, so it cannot drift back out of step with the data.
SOURCES = {
    "SuSiEx-ours": None,
    "MultiSuSiE": None,
    "MESuSiE": None,
    "SuSiE-merged": None,
    "FINEMAP-merged": None,
    # present only in the 2026-09 combined export (--export); no per-arm legacy file exists
    "RSparsePro": None, "CARMA": None, "EUR-SuSiE-RSS": None, "EUR-FINEMAP": None,
    "MACHINE": None,
}
# Ready-made alternatives for --source. SuSiE-merged's primary export is the
# min_abs_corr=0.1 run; a rerun at 0.5, which is what the manuscript's concordance
# numbers use and what the published benchmark itself was gated at. The choice moves 13
# of 122 loci from one credible set to ZERO, so it is not cosmetic.
ALIASES = {}
PRIMARY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "benchmark",
                       "pip_vectors_full250_qcexcl.csv.gz")

# The grouping of the method families. Arm order follows this list and a
# divider is drawn wherever the group changes. MACHINE is deliberately absent: across the 250
# loci its median top-CL is 0.155 and only 14 variants reach 0.5, so it is reported in the
# text as a dose-response finding rather than placed beside arms in their normal range.
# The two in-sample references (PGC-FINEMAP, Mawson SuSiE-RSS) are absent for a different
# reason -- their published tables contain ONLY credible-set members, so they have no
# out-of-CS population and cannot be split. They belong in the recovery-rate figure.
METHOD_GROUPS = [
    ("per-ancestry ref LD", ["SuSiEx-ours", "MultiSuSiE", "MESuSiE"]),
    ("FINEMAP",             ["FINEMAP-merged", "EUR-FINEMAP"]),
    ("SuSiE-RSS",           ["SuSiE-merged", "EUR-SuSiE-RSS"]),
    ("mismatch / latent",   ["RSparsePro", "CARMA"]),
]
GROUP_OF = {a: g for g, _arms in METHOD_GROUPS for a in _arms}
# the two long ones overprint their neighbours at this figure width
GROUP_WRAP = {"per-ancestry ref LD": "per-ancestry\nref LD",
              "mismatch / latent": "mismatch /\nlatent"}

# Option B encoding: HUE = TOOL, TEXTURE = LD SOURCE.
# The validated categorical theme has exactly eight slots and the skill's rule is that a
# ninth is never an invented hue, so nine distinct certified hues do not exist. Instead the
# two EUR-only arms share their tool's hue and are marked by hatch. That is not a
# workaround: it encodes the comparison the panel is making -- same tool, different LD.
# Fill OPACITY is already taken (solid = inside a credible set, washed = outside), which is
# why LD source goes on texture rather than on fill.
ARM_HUE = {
    "SuSiEx-ours":    PALETTE_[2],   # aqua
    "MultiSuSiE":     PALETTE_[0],   # blue
    "MESuSiE":        PALETTE_[1],   # orange
    "FINEMAP-merged": PALETTE_[4],   # magenta -- FINEMAP tool
    "EUR-FINEMAP":    PALETTE_[4],   # magenta + hatch
    "SuSiE-merged":   PALETTE_[3],   # yellow  -- SuSiE-RSS tool
    "EUR-SuSiE-RSS":  PALETTE_[3],   # yellow + hatch
    "RSparsePro":     PALETTE_[5],   # green
    "CARMA":          PALETTE_[7],   # red
    "MACHINE":        PALETTE_[6],   # violet (not in the 9-arm figure)
}
EUR_ONLY = {"EUR-FINEMAP", "EUR-SuSiE-RSS"}
HATCH = "///"
# Display names: the LD source reads better as a suffix than a prefix.
# Manuscript names: the merged-LD runs carry the tool's own name, as in the text; the grouping shows the LD.
DISPLAY = {"EUR-FINEMAP": "FINEMAP-EUR", "EUR-SuSiE-RSS": "SuSiE-RSS-EUR",
           "SuSiEx-ours": "SuSiEx", "FINEMAP-merged": "FINEMAP",
           "SuSiE-merged": "SuSiE-RSS"}
def disp(a):
    return DISPLAY.get(a, a)
METHOD_ORDER = [a for _, _arms in METHOD_GROUPS for a in _arms]


def boot_ratio(num_by_locus, den_by_locus, W):
    """Percentile CI of a ratio of sums under a bootstrap over LOCI (bench_sensfdr.py's
    design: variants inside a locus are in LD, so a binomial interval would be far too
    narrow)."""
    n, d = W @ num_by_locus, W @ den_by_locus
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(d > 0, n / np.where(d > 0, d, 1), np.nan)
    if np.all(np.isnan(r)):
        return np.nan, np.nan
    return (float(np.nanpercentile(r, 2.5)), float(np.nanpercentile(r, 97.5)))


def per_locus(df, loci_index, col):
    """Sum `col` within each locus, aligned to loci_index."""
    s = df.groupby("locus")[col].sum()
    return s.reindex(loci_index, fill_value=0).to_numpy(dtype=float)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--thresholds", default="0.1,0.5,0.9")
    ap.add_argument("--loci", default=None, help="restrict to a stratum file")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--export", default=None, metavar="DIR",
                    help="read the 2026-09 primary run's combined export from DIR "
                         "(primary_cs_membership.csv.gz + primary_machine_cl_matrix.csv.gz) "
                         "instead of the per-arm legacy files. MACHINE comes from the CL matrix "
                         "and is drawn UNSPLIT: run_machine_diag_b29.R never wrote fit@CS, so "
                         "it has no credible-set membership at all.")
    ap.add_argument("--arms", default=None,
                    help="comma-separated subset, in this order")
    ap.add_argument("--source", action="append", default=[], metavar="ARM=PATH",
                    help="override one arm's credible-set source; PATH may be an alias "
                         f"({', '.join(ALIASES)}). Repeatable. e.g. "
                         "--source 'SuSiE-merged=p05'")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    for spec in args.source:
        arm, path = spec.split("=", 1)
        if arm not in SOURCES:
            raise SystemExit(f"--source: unknown arm {arm!r}; known: {list(SOURCES)}")
        SOURCES[arm] = ALIASES.get(path, path)
        print(f"[info] {arm} credible sets from {SOURCES[arm]}")

    thresholds = [float(x) for x in args.thresholds.split(",")]
    keep, _ = bench_cs_loci(BENCH_RAW)
    if args.loci:
        keep = keep & set(open(args.loci).read().split())

    b = load_pips(BENCH)
    b = b[b.locus.isin(keep)][["locus", "key", "pip"]].rename(columns={"pip": "pb"})
    b = b.groupby(["locus", "key"], as_index=False).pb.max()

    # --- combined-export loader (2026-09 primary run) -------------------------------
    EXPORT = {}
    if args.export:
        # The export is PER ARM as of the revised export: primary_<arm>_cs_purity.csv.gz
        # carries membership (a superset of the old combined primary_cs_membership.csv.gz,
        # adding common_purity), and MACHINE's M x K matrix now rides in its own pip file as
        # CL_EUR/CL_EAS rather than a separate matrix file. Globbed rather than listed so a
        # new arm appears without editing this loader.
        for f in sorted(glob.glob(os.path.join(args.export, "primary_*_cs_purity.csv.gz"))):
            arm_ = os.path.basename(f)[len("primary_"):-len("_cs_purity.csv.gz")]
            g = pd.read_csv(f, low_memory=False)
            EXPORT[arm_] = g[["locus", "rsid", "bp", "a1", "a2", "pip", "cs_ids_all"]].copy()
        mp = os.path.join(args.export, "primary_MACHINE_pip.csv.gz")
        if os.path.exists(mp):
            mm = pd.read_csv(mp)
            # The row max is OUR aggregation of MACHINE's M x K output, not the paper's --
            # it reports CL_EUR and CL_EAS separately and states no pooling rule. Flagged in
            # the caption; do not let it harden into "MACHINE's PIP".
            mm["pip"] = mm[["CL_EUR", "CL_EAS"]].max(axis=1)
            mm["cs_ids_all"] = pd.NA          # no fit@CS was ever written -> unsplit
            EXPORT["MACHINE"] = mm[["locus", "rsid", "bp", "a1", "a2", "pip", "cs_ids_all"]]
        if not EXPORT:
            raise SystemExit(f"no primary_*_cs_purity.csv.gz under {args.export}")

    if args.export:
        avail = list(EXPORT)
        arms = ([a.strip() for a in args.arms.split(",")] if args.arms
                else [a for a in METHOD_ORDER if a in avail])
        bad = [a for a in arms if a not in EXPORT]
        if bad:
            raise SystemExit(f"not in the export: {bad}; available: {avail}")
        # unsplit = no credible-set membership at all
        for a in arms:
            if EXPORT[a].cs_ids_all.notna().sum() == 0:
                SOURCES[a] = None
            else:
                SOURCES[a] = "export"
    else:
        arms = order_arms(list(SOURCES))
    rng = np.random.default_rng(args.seed)
    rows = []
    primary = None
    for arm in arms:
        if args.export:
            d = EXPORT[arm].copy()
            d = d[d.locus.isin(keep)]
            d["key"] = variant_key(d)
        elif SOURCES[arm] is None:
            if primary is None:
                primary = load_pips(PRIMARY)
            d = primary[primary.method == arm][["locus", "key", "pip"]].copy()
            d = d[d.locus.isin(keep)]
            d["cs_ids_all"] = pd.NA          # unused; this arm is not split
        elif True:
            d = pd.read_csv(SOURCES[arm], low_memory=False)
            d = d[d.locus.isin(keep)].copy()
            d["key"] = variant_key(d)
        u = b.merge(d[["locus", "key", "pip", "cs_ids_all"]], on=["locus", "key"],
                    how="inner")
        loci = np.sort(u.locus.unique())
        W = rng.multinomial(len(loci), np.full(len(loci), 1 / len(loci)),
                            size=args.n_boot).astype(float)
        for t in thresholds:
            u = u.assign(is_pos=(u.pb > t).astype(int))
            nPB_loc = per_locus(u, loci, "is_pos")
            nPB = int(nPB_loc.sum())
            # None == "not split": every call in one group (FINEMAP-merged)
            groups = [None] if SOURCES.get(arm) is None else [True, False]
            for inside in groups:
                sel = (u[u.pip > t] if inside is None
                       else u[(u.pip > t) & (u.cs_ids_all.notna() == inside)]).copy()
                sel["tp"] = (sel.pb > t).astype(int)
                sel["fp"] = 1 - sel.tp
                sel["n"] = 1
                nc = len(sel)
                tp = int(sel.tp.sum())
                fdr = (nc - tp) / nc if nc else np.nan
                rec = tp / nPB if nPB else np.nan
                if nc:
                    flo, fhi = boot_ratio(per_locus(sel, loci, "fp"),
                                          per_locus(sel, loci, "n"), W)
                    rlo, rhi = boot_ratio(per_locus(sel, loci, "tp"), nPB_loc, W)
                else:
                    flo = fhi = rlo = rhi = np.nan
                rows.append(dict(arm=arm, threshold=t, inside=inside, n_calls=nc,
                                 tp=tp, n_PB=nPB, fdr=fdr, fdr_lo=flo, fdr_hi=fhi,
                                 recall=rec, recall_lo=rlo, recall_hi=rhi))
    R = pd.DataFrame(rows)

    n_th = len(thresholds)
    fig, axes = plt.subplots(2, n_th, figsize=(3.6 * n_th + 1.4, 7.0), dpi=200,
                             sharey="row")
    if n_th == 1:
        axes = axes.reshape(2, 1)
    w = 0.34
    for ci, t in enumerate(thresholds):
        aF, aR = axes[0, ci], axes[1, ci]
        totals = []
        for ai, arm in enumerate(arms):
            col = ARM_HUE.get(arm, "#5b6b7a")
            sub = R[(R.arm == arm) & (R.threshold == t)]
            unsplit = sub.inside.isna().all()
            if unsplit:
                # FINEMAP-merged: one centred bar, no split available
                g = sub.iloc[0]
                for ax, val, lo, hi in ((aF, g.fdr, g.fdr_lo, g.fdr_hi),
                                        (aR, g.recall, g.recall_lo, g.recall_hi)):
                    ax.bar(ai, val, width=w * 1.5, color=col, alpha=0.95, zorder=3,
                           linewidth=0, hatch="\\\\")
                    if not np.isnan(lo):
                        ax.errorbar(ai, val,
                                    yerr=[[max(val - lo, 0)], [max(hi - val, 0)]],
                                    fmt="none", ecolor=INK, elinewidth=0.9, capsize=2.2,
                                    capthick=0.9, zorder=5)
                aF.annotate(f"{int(g.n_calls)}", (ai, max(g.fdr_hi, g.fdr)),
                            xytext=(0, 3), textcoords="offset points", ha="center",
                            va="bottom", fontsize=6.2, color=INK_MUTED)
                aR.annotate(f"total {g.recall:.2f}", (ai, max(g.recall_hi, g.recall)),
                            xytext=(0, 4), textcoords="offset points", ha="center",
                            va="bottom", fontsize=6.6, color=INK_MUTED)
                continue
            gi = sub[sub.inside == True].iloc[0]
            go = sub[sub.inside == False].iloc[0]

            # --- FDR: paired bars, inside solid / outside washed out + outlined
            for k, (g, dx, kw) in enumerate((
                    (gi, -w / 2, dict(color=col, alpha=0.95, linewidth=0,
                                    hatch=HATCH if arm in EUR_ONLY else None,
                                    edgecolor="white" if arm in EUR_ONLY else None)),
                    (go, +w / 2, dict(color=col, alpha=0.30, edgecolor=col, linewidth=1.1,
                                    hatch=HATCH if arm in EUR_ONLY else None)))):
                if g.n_calls == 0:
                    aF.annotate("no calls", (ai + dx, 0.02), ha="center", va="bottom",
                                fontsize=6.2, color=INK_MUTED, rotation=90)
                    continue
                aF.bar(ai + dx, g.fdr, width=w * 0.9, zorder=3, **kw)
                if not np.isnan(g.fdr_lo):
                    aF.errorbar(ai + dx, g.fdr,
                                yerr=[[max(g.fdr - g.fdr_lo, 0)], [max(g.fdr_hi - g.fdr, 0)]],
                                fmt="none", ecolor=INK, elinewidth=0.9, capsize=2.2,
                                capthick=0.9, zorder=5)
                aF.annotate(f"{int(g.n_calls)}", (ai + dx, max(g.fdr_hi, g.fdr)),
                            xytext=(0, 3), textcoords="offset points", ha="center",
                            va="bottom", fontsize=6.2, color=INK_MUTED)

            # --- recall: paired bars, same encoding as the FDR row. The two groups
            # still partition the true positives, so they SUM to the arm's total recall,
            # printed above the pair since side-by-side bars no longer show it directly.
            tops = []
            for g, dx, kw in ((gi, -w / 2, dict(color=col, alpha=0.95, linewidth=0,
                                    hatch=HATCH if arm in EUR_ONLY else None,
                                    edgecolor="white" if arm in EUR_ONLY else None)),
                              (go, +w / 2, dict(color=col, alpha=0.30, edgecolor=col,
                                                linewidth=1.1))):
                if g.n_calls == 0:
                    aR.annotate("no calls", (ai + dx, 0.02), ha="center", va="bottom",
                                fontsize=6.2, color=INK_MUTED, rotation=90)
                    continue
                aR.bar(ai + dx, g.recall, width=w * 0.9, zorder=3, **kw)
                if not np.isnan(g.recall_lo):
                    aR.errorbar(ai + dx, g.recall,
                                yerr=[[max(g.recall - g.recall_lo, 0)],
                                      [max(g.recall_hi - g.recall, 0)]],
                                fmt="none", ecolor=INK, elinewidth=0.9, capsize=2.2,
                                capthick=0.9, zorder=5)
                    tops.append(max(g.recall_hi, g.recall))
                else:
                    tops.append(g.recall)
            tot = gi.recall + (go.recall if not np.isnan(go.recall) else 0)
            # stagger alternate labels: with nine arms the totals collide whenever two
            # neighbouring bars are the same height, which at alpha=0.9 is most of them
            totals.append((ai, max(tops) if tops else 0.0, tot))

        for k, (ai_, top_, tot_) in enumerate(totals):
            aR.annotate(f"{tot_:.2f}", (ai_, top_), xytext=(0, 4 + 11 * (k % 2)),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=6.6, color=INK_MUTED)
        for ax in (aF, aR):
            # The grouping: a divider wherever the family changes, label above. Without it
            # nine arms read as one undifferentiated row and the comparison the paper makes
            # (per-ancestry vs merged vs EUR vs mismatch-aware) is invisible.
            prev = None
            for i, _a in enumerate(arms):
                g = GROUP_OF.get(_a)
                if prev is not None and g != prev:
                    ax.axvline(i - 0.5, color=GRID, lw=1.2, zorder=1)
                prev = g
            runs = []
            for _a in arms:
                g = GROUP_OF.get(_a, "")
                if not runs or runs[-1][0] != g:
                    runs.append([g, 0])
                runs[-1][1] += 1
            at = 0
            for g, n in runs:
                if ax is aR and g:
                    ax.annotate(GROUP_WRAP.get(g, g), (at + (n - 1) / 2, -0.34),
                                xycoords=("data", "axes fraction"), ha="center",
                                va="top", fontsize=7.6, color=INK)
                at += n
            ax.set_xticks(range(len(arms)))
            ax.set_xticklabels([disp(a_) for a_ in arms] if ax is aR else [],
                               fontsize=7.6, color=INK, rotation=20, ha="right")
            ax.set_xlim(-0.65, len(arms) - 0.35)
            ax.grid(False)
            ax.set_axisbelow(True)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            for s in ("left", "bottom"):
                ax.spines[s].set_color(GRID)
            ax.tick_params(colors=INK_MUTED, labelsize=8)
        # |P_B| is the BENCHMARK's own positive count, not any arm's join-restricted
        # view of it, so the denominator is one fixed number per column. The arms differ
        # in how much of it they can see, which is reported as the join range.
        npb = int((b.pb > t).sum())
        lo_ = int(R[R.threshold == t].n_PB.min()); hi_ = int(R[R.threshold == t].n_PB.max())
        rng = f"{lo_}" if lo_ == hi_ else f"{lo_}\u2013{hi_}"
        aF.set_title(f"PIP > {t:g}    |P_B| = {npb}   (arm joins {rng})",
                     fontsize=9.5, color=INK, pad=6)
    axes[0, 0].set_ylabel("FDR", fontsize=10.5, color=INK)
    axes[1, 0].set_ylabel("Sensitivity", fontsize=10.5, color=INK)
    axes[0, 0].set_ylim(0, 1.08)
    axes[1, 0].set_ylim(0, 1.08)

    h = [plt.Line2D([], [], color=INK_MUTED, lw=7, alpha=0.95,
                    label="call INSIDE the arm's own CS"),
         plt.Line2D([], [], color=INK_MUTED, lw=7, alpha=0.30,
                    label="call OUTSIDE any of its CS")]
    import matplotlib.patches as mpatches
    # a Line2D cannot show the hatch the FINEMAP bars actually carry; use a Patch so the
    # legend swatch matches what is drawn
    # Built from what was ACTUALLY drawn, not from a fixed list. The "not split" entry
    # described the old state where FINEMAP exported no credible-set membership; the
    # 2026-09 export supplies it for all nine arms, so that category no longer appears
    # here and a hardcoded entry would be describing a bar that is not on the figure.
    h = [mpatches.Patch(facecolor=INK_MUTED, alpha=0.95, linewidth=0,
                        label="call INSIDE the arm's own CS"),
         mpatches.Patch(facecolor=INK_MUTED, alpha=0.30, edgecolor=INK_MUTED,
                        linewidth=1.1, label="call OUTSIDE any of its CS")]
    if any(a in EUR_ONLY for a in arms):
        h.append(mpatches.Patch(facecolor=INK_MUTED, alpha=0.95, linewidth=0, hatch=HATCH,
                                edgecolor="white",
                                label="hatched = EUR-only reference LD (hue = tool)"))
    if any(SOURCES.get(a) is None for a in arms):
        h.append(mpatches.Patch(facecolor=INK_MUTED, alpha=0.95, linewidth=0, hatch="\\\\",
                                label="not split (no CS membership)"))
    fig.legend(handles=h, frameon=False, fontsize=7.8, ncol=len(h), labelcolor=INK,
               loc="lower center", bbox_to_anchor=(0.5, -0.045))

    n_loci = len(keep)
    nb = int(R.n_boot.iloc[0]) if "n_boot" in R else 2000
    ci = 0.95
    fig.suptitle("FDR and sensitivity split by credible-set membership of the call  ·  "
                 f"{n_loci} loci  ·  primary run  ·  s = α, matched",
                 fontsize=12.5, color=INK, x=0.004, ha="left", y=1.005)
    note = (
        "Sensitivity = TP/(TP+FN), the fraction of the benchmark's positives we recover; "
        "higher is better. FDR = FP/(TP+FP), the fraction of our calls the benchmark does "
        "not support; lower is better.\n"
        "BENCHMARK: the published cross-ancestry SuSiEx fine-mapping of these same PGC3-SCZ "
        "loci (Yuan et al.), reconstructed per variant from its per-locus LogBF output. A "
        "benchmark positive is a variant with benchmark PIP > \u03b1, using the same \u03b1 "
        "applied to our arm (matched pinning), so each column is a symmetric question. The 57 "
        "loci where the benchmark formed no credible set are excluded: every PIP there is "
        "exactly 0 by construction, not a judgement of noise.\n"
        "THE SPLIT: each arm's calls are separated by whether the call lies inside one of that "
        "arm's OWN credible sets. The benchmark is credible-set restricted by construction "
        "(nothing outside one of its sets exceeds PIP 0.0077), so out-of-CS calls are scored "
        "against a reference with almost no way to agree \u2014 a format mismatch affecting only "
        "the arms that export full posteriors. "
        "COLOUR: hue is the tool, hatch is EUR-only reference LD, so FINEMAP-merged and "
        "FINEMAP-EUR share a hue, as do SuSiE-merged and SuSiE-EUR. The validated "
        "colourblind-safe theme has eight slots and there are nine arms; rather than "
        "invent an uncertified ninth hue, hue is bound to the TOOL, which needs only "
        "seven, and the two EUR-only arms are separated from their merged-reference "
        "counterparts by texture instead. "
        "THE FINEMAP ARMS: their credible sets are only ~1.3x the SuSiE family's here "
        "(median 36 vs 26-29), but they form one at EVERY locus (250/250 vs 140-198), so "
        "29% of their variants sit inside a CS against ~5% for the others. Inside-a-CS is "
        "a weaker statement for them, which is why their out-of-CS bars read no calls.\n"
        f"Error bars are {ci:.0%} percentile intervals from {nb:,} bootstrap resamples of the "
        f"{n_loci} loci (loci, not variants \u2014 variants within a locus are in LD). Numbers "
        "above the FDR bars are call counts; watch them where a ratio rests on <20 calls. "
        "Primary run; the `_bars` figures use the L=5 sweep and differ.")
    fig.text(0.004, -0.085, note, fontsize=8, color=INK_MUTED, va="top", wrap=True)
    fig.tight_layout()
    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    print(f"[done] wrote {args.out}")
    out = args.out.rsplit(".", 1)[0] + "_summary.csv"
    R.to_csv(out, index=False)
    print(f"[done] wrote {out}\n")
    print(R.to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
