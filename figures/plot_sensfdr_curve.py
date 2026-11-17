#!/usr/bin/env python3
"""Sensitivity/FDR curve with the label FIXED at benchmark PIP > 0.1, in-CS calls only.

Different question from the bar figures. Those use MATCHED pinning -- a benchmark positive
is a variant above the same alpha we apply to ourselves -- so both sides move together and
each column is a self-contained symmetric comparison. Here the truth set is fixed once
(benchmark PIP > 0.1) and only OUR threshold sweeps, which is the ordinary
precision-recall setup: one label, one caller, a curve.

Fixing the label is what makes the arms comparable along their whole range instead of at
three points, and it is why the curve can cross -- an arm can be better at high threshold
and worse at low, which a bar at a single alpha cannot show.

ONLY IN-CREDIBLE-SET CALLS ARE SCORED. Out-of-CS calls are a different object: the
benchmark is credible-set restricted by construction (nothing outside one of its sets
exceeds PIP 0.0077), so scoring our out-of-CS tail against it measures a format mismatch
rather than a disagreement. The sensitivity denominator stays every benchmark positive in
the joined set, so an arm that forms fewer credible sets is charged for it rather than
excused.

Hue is the tool, as in every other figure; the two EUR-only arms take their tool's hue with
a dashed line, so hue = tool and linestyle = LD source. Markers sit at our PIP = 0.1, 0.5
and 0.9 so the bar figures' columns can be located on the curve.

Usage:
  python3 plot_sensfdr_curve.py --out figures/sensfdr_curve_inCS_193.png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# every label bold: ticks and legends inherit font.weight
plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                     "axes.titleweight": "bold"})
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_sensfdr import bench_cs_loci, BENCH, BENCH_RAW
from sensitivity_fdr import load_pips, variant_key
from plot_cs_split import ARM_HUE, METHOD_ORDER, EUR_ONLY, disp
from plot_sensitivity_fdr import INK, INK_MUTED, GRID

MARKS = [0.1, 0.5, 0.9]
# The 122 stratum is defined with >= 0.1 (derive_confident_stratum.py:81) but is LABELLED
# as > 0.1 everywhere, for one consistent form across the figures. The two are the
# same 122 loci: the lowest max PIP among the 122 is 0.1006 and the highest among the 53
# excluded is 0.0973, so no locus sits on the boundary.
# two lines: on one line the 193 title ran into the 175 panel's title
STRATUM_TITLE = {
    "193": "193 loci — benchmark formed a CS\n136 carry a variant at PIP > 0.1",
    "175": "175 loci — single benchmark CS",
    "122": "122 loci — single CS, top PIP > 0.1",
}
# marker, size, filled -- three DISTINCT glyphs; hollow-circle for both 0.1 and 0.9 made
# the two legend entries identical and therefore useless
MARK_STYLE = {0.1: ("o", 22, False), 0.5: ("o", 34, True), 0.9: ("s", 24, False)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--export", default=os.path.join(here, "..", "results", "finemap", "main"))
    ap.add_argument("--no-marks", dest="marks", action="store_false",
                    help="drop the PIP 0.1/0.5/0.9 glyphs from each curve")
    ap.add_argument("--strata", default="193",
                    help="comma-separated; one stacked panel each, in the order given")
    ap.add_argument("--all-variants", action="store_true",
                    help="score EVERY call, not only those inside the arm's credible set")
    ap.add_argument("--exclude", default="EUR-FINEMAP,EUR-SuSiE-RSS",
                    help="arms to leave out; the two EUR-only arms by default, since they "
                         "belong with the published single-ancestry comparators rather than "
                         "with the seven multi-ancestry/merged arms")
    ap.add_argument("--label-pip", type=float, default=0.1,
                    help="a benchmark variant above this is a positive; FIXED, not swept")
    ap.add_argument("--pr", action="store_true",
                    help="plot precision (1 - FDR) against recall (= sensitivity) instead "
                         "of FDR against sensitivity, and label each arm with its AUPRC. "
                         "Same numbers, the requested orientation (MultiSuSiE Fig 2c).")
    ap.add_argument("--caption", action="store_true",
                    help="print the caption under the figure; default writes it to a "
                         "sidecar .txt so captions live in one place")
    ap.add_argument("--suffix", default="-noambig",
                    help="export label suffix, e.g. -noambig for the no-palindromic tree (results/finemap/main)")
    ap.add_argument("--s-flag", type=float, default=None,
                    help="drop every locus whose z/LD mismatch s (susieR estimate_s_rss with n, "
                         "s_susier_per_ancestry.csv) exceeds this in EITHER ancestry -- one uniform "
                         "rule, applied to all loci (review round 3)")
    ap.add_argument("--s-table", default=os.path.join(here, "..", "results", "s_statistic",
                                                      "s_susier_per_ancestry.csv"))
    ap.add_argument("--s-col", default="s_susier_n",
                    help="column of --s-table holding s (s_legacy: the earlier estimator)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    noambig = "noambig" in args.suffix

    strata = [x.strip() for x in args.strata.split(",") if x.strip()]
    full, _ = bench_cs_loci(BENCH_RAW)
    n_flagged = 0
    if args.s_flag is not None:
        smax = pd.read_csv(args.s_table).groupby("locus")[args.s_col].max()
        flagged = set(smax[smax > args.s_flag].index)
        n_flagged = len(set(full) & flagged)
        full = type(full)(l for l in full if l not in flagged)
        for st in list(STRATUM_TITLE):
            STRATUM_TITLE[st] = (f"{len(full)} loci — benchmark formed a CS\n"
                                 f"{n_flagged} loci with z/LD mismatch s > {args.s_flag:g} removed")
    b_all = load_pips(BENCH)
    b_all = b_all[b_all.locus.isin(full)][["locus", "key", "pip"]].rename(
        columns={"pip": "pb"})
    b_all = b_all.groupby(["locus", "key"], as_index=False).pb.max()

    def stratum_loci(st):
        if st == "193":
            return set(full)
        f = os.path.join(here, "..", "inputs", {"122": "loci_confident122.txt",
                                "175": "loci_single175.txt"}[st])
        return set(full) & {l.strip() for l in open(f) if l.strip()}

    drop = {x.strip() for x in args.exclude.split(",") if x.strip()}
    arms = [a for a in METHOD_ORDER
            if a not in drop
            and os.path.exists(os.path.join(args.export, f"primary_{a}{args.suffix}_cs_purity.csv.gz"))]
    grid = np.unique(np.concatenate([np.arange(0.01, 1.0, 0.01), np.array(MARKS)]))
    raw = {}
    for arm in arms:
        d = pd.read_csv(os.path.join(args.export, f"primary_{arm}{args.suffix}_cs_purity.csv.gz"),
                        low_memory=False)
        d = d[d.locus.isin(full)].copy()
        d["key"] = variant_key(d)
        raw[arm] = d[["locus", "key", "pip", "cs_ids_all"]]
    rows = []
    for st in strata:
        keep_st = stratum_loci(st)
        b = b_all[b_all.locus.isin(keep_st)]
        for arm in arms:
            u = b.merge(raw[arm][raw[arm].locus.isin(keep_st)], on=["locus", "key"],
                        how="inner")
            pos = (u.pb > args.label_pip).to_numpy()
            nP = int(pos.sum())
            incs = u.cs_ids_all.notna().to_numpy()
            pc = u.pip.to_numpy()
            for s in grid:
                call = (pc > s) if args.all_variants else ((pc > s) & incs)
                nC = int(call.sum()); tp = int((call & pos).sum())
                rows.append(dict(stratum=st, n_loci=len(keep_st), arm=arm, s=float(s),
                                 n_calls=nC, tp=tp, n_PB=nP,
                                 sensitivity=tp / nP if nP else np.nan,
                                 fdr=(nC - tp) / nC if nC else np.nan,
                                 precision=tp / nC if nC else np.nan))
    R = pd.DataFrame(rows)

    # AUPRC, two conventions, both written to the CSV; the legend shows `auprc`.
    #   auprc          -- trapezoid over recall in [0, 1], precision taken as 0 beyond the
    #                     arm's reach. An arm that forms few credible sets cannot buy a
    #                     high score by stopping early, which matches the sensitivity
    #                     denominator: fewer sets costs an arm rather than excusing it.
    #   auprc_observed -- trapezoid over the recall range the arm actually spans, divided
    #                     by that range. Read it as "precision while it is still calling",
    #                     not as an area, and never compare it across arms with different
    #                     reach.
    def _auprc(g):
        g = g.sort_values("s", ascending=False)
        m = g.precision.notna() & g.sensitivity.notna()
        x, y = g.sensitivity[m].to_numpy(), g.precision[m].to_numpy()
        if len(x) < 2:
            return np.nan, np.nan
        o = np.argsort(x); x, y = x[o], y[o]
        full = np.trapezoid(np.r_[y[0], y, 0.0], np.r_[0.0, x, x[-1]])
        span = x[-1] - x[0]
        return full, (np.trapezoid(y, x) / span if span > 0 else np.nan)

    A = {}
    for st in strata:
        for arm in arms:
            A[(st, arm)] = _auprc(R[(R.stratum == st) & (R.arm == arm)])
    R["auprc"] = [A[(r.stratum, r.arm)][0] for r in R.itertuples()]
    R["auprc_observed"] = [A[(r.stratum, r.arm)][1] for r in R.itertuples()]

    # the legend lives outside the axes; reserve figure width for it explicitly rather
    # than relying on the tight bbox, which was clipping the longer method names
    # NOT sharex: each panel is meant to stand alone, with its own axis labels and its own
    # legend, so a single panel can be lifted out of the figure without losing its key
    horiz = len(strata) > 1
    if horiz:
        fig, axes = plt.subplots(1, len(strata), figsize=(6.2 * len(strata), 5.6), dpi=300,
                                 sharex=True, sharey=True)
        fig.subplots_adjust(wspace=0.22, top=0.90)
    else:
        fig, axes = plt.subplots(1, 1, figsize=(9.8, 5.4), dpi=300)
        fig.subplots_adjust(right=0.745, top=0.95)
    axes = np.atleast_1d(axes)
    for ax, st in zip(axes, strata):
        Rs = R[R.stratum == st]
        for arm in arms:
            g = Rs[Rs.arm == arm].sort_values("s")
            ycol = "precision" if args.pr else "fdr"
            m = g[ycol].notna()
            ax.plot(g.sensitivity[m], g[ycol][m], color=ARM_HUE[arm], lw=1.7,
                    ls="--" if arm in EUR_ONLY else "-", zorder=3, solid_capstyle="round")
            if args.marks:
                for sv in MARKS:
                    r = g[np.isclose(g.s, sv)]
                    if len(r) and pd.notna(r[ycol].iloc[0]):
                        mk, sz, fc = MARK_STYLE[sv]
                        ax.scatter(r.sensitivity, r[ycol], s=sz, marker=mk,
                                   facecolor=ARM_HUE[arm] if fc else "white",
                                   edgecolor=ARM_HUE[arm], linewidths=1.3, zorder=5)
        # the BENCHMARK's own positive count, not any arm's join-restricted view of it
        npb = int((b_all[b_all.locus.isin(stratum_loci(st))].pb > args.label_pip).sum())
        lo_, hi_ = int(Rs.n_PB.min()), int(Rs.n_PB.max())
        rng = f"{lo_}" if lo_ == hi_ else f"{lo_}\u2013{hi_}"
        ax.set_title(STRATUM_TITLE.get(st, f"{int(Rs.n_loci.iloc[0])} loci"),
                     fontsize=13.5, color=INK, loc="left", pad=8)
        ax.set_ylabel("Precision" if args.pr else "FDR", fontsize=14.5, color=INK)
        ax.set_xlim(0, 1.0); ax.set_ylim(0, 1.0)
        # NO NO-SKILL LINE, deliberately. A first version drew one at n_PB divided by the
        # call count at the loosest threshold, which is not the prevalence of anything --
        # the numerator counts positives in the whole joined set and the denominator only
        # the variants called at PIP > 0.01. It read 0.12; the real positive rate is
        # 0.0028 over the joined set. Drawn correctly the line would hug the axis and say
        # nothing, and there is no single right value anyway: a caller restricted to
        # credible sets faces a prevalence of 0.0095 (FINEMAP-merged) to 0.058
        # (RSparsePro) depending on how widely the arm scopes its sets, so one averaged
        # line would flatter the arms with the tightest sets. The number is in the caption
        # instead.
        ax.grid(False)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=INK_MUTED, labelsize=13.5)
    # Each panel stands alone: its own x label and its own legends, so a single panel can
    # be lifted out of the figure without losing its key.
    any_eur = any(a in EUR_ONLY for a in arms)
    extra = []
    for ax in axes:
        # the axis titles are just the quantity; what defines them lives in the caption
        ax.set_xlabel("Recall" if args.pr else "Sensitivity", fontsize=14.5, color=INK)
        ax.tick_params(labelbottom=True)
        st_ = strata[list(axes).index(ax)]
        # precision-recall legends rank methods by AUPRC, high to low; ties keep method order
        leg_arms = sorted(arms, key=lambda a: -round(A[(st_, a)][0], 2)) if args.pr else arms
        h = [mlines.Line2D([], [], color=ARM_HUE[a], lw=1.9,
                           ls="--" if a in EUR_ONLY else "-",
                           label=(f"{disp(a)}  {A[(st_, a)][0]:.2f}" if args.pr
                                  else disp(a))) for a in leg_arms]
        if not horiz:
            lg = ax.legend(handles=h, frameon=False, fontsize=11.5, labelcolor=INK,
                           loc="center left", bbox_to_anchor=(1.02, 0.5),
                           bbox_transform=ax.transAxes, borderaxespad=0, alignment="left",
                           title=("Method" + (" (dashed = EUR-only LD)" if any_eur else "")
                                  + ("     AUPRC" if args.pr else "")),
                           title_fontsize=12.2)
            ax.add_artist(lg)
            extra.append(lg)
        if args.marks:
            hm = [mlines.Line2D([], [], color=INK_MUTED, marker=MARK_STYLE[t][0], ls="",
                                markerfacecolor=INK_MUTED if MARK_STYLE[t][2] else "white",
                                markeredgewidth=1.3, markersize=6 if MARK_STYLE[t][2] else 5,
                                label=f"PIP > {t:g}") for t in MARKS]
            # upper-left is empty in the FDR orientation and full of curves in the
            # precision-recall one, so the marker key moves to the lower left there --
            # lifted clear of the no-skill line, which sits near 0.12
            kw = (dict(loc="lower left", bbox_to_anchor=(0.012, 0.13),
                       bbox_transform=ax.transAxes, borderaxespad=0)
                  if args.pr else dict(loc="upper left", borderaxespad=0.8))
            extra.append(ax.legend(handles=hm, frameon=False, fontsize=11.0,
                                   labelcolor=INK, alignment="left", **kw))
    if horiz:
        extra.append(fig.legend(handles=h, frameon=False, fontsize=12.5, labelcolor=INK,
                                ncol=len(arms), loc="lower center",
                                bbox_to_anchor=(0.5, -0.13), alignment="left",
                                title="Method" + (" (dashed = EUR-only LD)"
                                                  if any_eur else ""),
                                title_fontsize=13))


    # caption numbers, computed from this run rather than carried over from v1
    st0 = strata[0]
    bst = b_all[b_all.locus.isin(stratum_loci(st0))]
    n_pos_loci = int(bst[bst.pb > args.label_pip].locus.nunique())
    n_loci_st = len(stratum_loci(st0))
    prev = {}
    for arm in arms:
        u = bst.merge(raw[arm], on=["locus", "key"], how="inner")
        pos = u.pb > args.label_pip
        inc = u.cs_ids_all.notna()
        prev[arm] = (pos.mean(), (pos & inc).sum() / max(inc.sum(), 1))
    p_join = float(np.mean([v[0] for v in prev.values()]))
    lo_a = min(prev, key=lambda a: prev[a][1]); hi_a = max(prev, key=lambda a: prev[a][1])
    note = (
        ("Precision-recall as each arm's PIP threshold sweeps, with the label fixed at "
         if args.pr else
         "Sensitivity vs FDR as each arm's PIP threshold sweeps, with the label fixed at ")
        + f"benchmark PIP > {args.label_pip:g}"
        + (", scoring all calls in and out of credible sets. " if args.all_variants
           else ", scoring credible-set calls only. ")
        + "Each line is one arm, swept over ITS OWN PIP threshold from 0.01 to 0.99; the truth "
        f"set is fixed throughout at benchmark PIP > {args.label_pip:g}. That is the "
        "difference from the bar figures, which pin the benchmark at whatever alpha they "
        "apply to us: matched pinning makes each column symmetric, a fixed label makes the "
        "arms comparable along their whole range and lets curves cross where a bar at one "
        "alpha cannot show it. "
        + ("EVERY CALL IS SCORED here, in or out of a credible set. Read the out-of-CS part "
           "with care: the benchmark is credible-set restricted by construction (nothing "
           "outside one of its sets exceeds PIP 0.0077), so a call of ours out there is "
           "scored against a reference with almost no way to agree, and the FDR it earns is "
           "partly a format mismatch. The in-CS companion figure is the cleaner comparison. "
           if args.all_variants else
           "ONLY IN-CREDIBLE-SET CALLS are scored: the benchmark is credible-set restricted "
           "by construction, so our out-of-CS tail measures a format mismatch rather than a "
           "disagreement. "
           + ("PALINDROMIC VARIANTS: this run excludes A/T and C/G variants from every arm's inputs, as the benchmark "
              "does (Yuan et al. removed them from their sumstats upstream), so both sides fit the same kind of variant set. "
              if noambig else
              "STRAND-AMBIGUOUS VARIANTS: the benchmark contains NO palindromic (A/T, C/G) variant at all -- Yuan et al. "
              "stripped them from their sumstats upstream -- while about 15% of ours are palindromic. The FDR axis is "
              "unaffected (the join is an inner join, so no palindromic call is ever scored), but sensitivity can be "
              "depressed where a palindromic variant absorbs posterior mass; the no-palindromic run is the primary analysis. "))
        + (f"Note that only {n_pos_loci} of the {n_loci_st} loci carry a benchmark variant above "
           f"{args.label_pip:g} at all; the other {n_loci_st - n_pos_loci} contribute no positive and can "
           "therefore only add false discoveries. ")
        + (f"Positives are rare: a benchmark variant above the label threshold is {100 * p_join:.2f}% "
           f"of the joined variant set, and {100 * prev[lo_a][1]:.1f}% ({disp(lo_a)}) to {100 * prev[hi_a][1]:.1f}% "
           f"({disp(hi_a)}) of the credible-set rows an arm can actually call, the spread reflecting how widely each "
           "arm scopes its credible sets rather than any property of the truth. No "
           "no-skill line is drawn, because at that rate it would hug the axis and "
           "because no single value applies across arms. "
           if args.pr else "")
        + ("Precision is 1 - FDR and recall is sensitivity, so this is the same scoring "
           "as the FDR-axis version, re-oriented. AUPRC in the legend is the trapezoid "
           "over recall in [0, 1] with precision taken as zero beyond the arm's reach, so "
           "an arm that forms few credible sets cannot buy a high score by stopping "
           "early; `auprc_observed` in the companion CSV instead averages precision over "
           "the recall range the arm actually spans and is NOT comparable across arms. "

           if args.pr else "")
        + "The sensitivity denominator is every benchmark positive in the joined set, so "
        "forming fewer credible sets costs an arm rather than excusing it, which is why "
        "several curves stop short of the right-hand edge. "
        + ("The two EUR-only arms are excluded; they belong with the published "
           "single-ancestry comparators rather than with these seven. " if drop else "")
        + (f"LD-MISMATCH FILTER: the {n_flagged} benchmark loci whose z/LD mismatch parameter s "
           f"(susieR estimate_s_rss, computed per ancestry on SuSiEx's inputs: EUR z vs GCTB, EAS z vs "
           f"AoU EAS) exceeds {args.s_flag:g} in either ancestry are removed, by one rule applied to "
           "every locus. " if args.s_flag is not None else "")
        + "Hue is the tool."
        + (" Markers: hollow circle = PIP > 0.1, filled = > 0.5, hollow square = > 0.9."
           if args.marks else ""))
    if args.caption:
        fig.text(0.004, -0.02, note, fontsize=9.7, color=INK_MUTED, va="top", wrap=True,
                 transform=fig.transFigure)
    else:
        # captions live in figures/CAPTIONS.md; this sidecar is what that file is built from
        cp = args.out.rsplit(".", 1)[0] + "_caption.txt"
        with open(cp, "w") as fh:
            fh.write(" ".join(note.split()) + "\n")
        print(f"[done] wrote {cp}")
    fig.savefig(args.out, bbox_inches="tight", facecolor="white",
                bbox_extra_artists=extra)
    print(f"[done] wrote {args.out}")
    out = args.out.rsplit(".", 1)[0] + "_curve.csv"
    R.to_csv(out, index=False)
    piv = (R[R.s.isin(MARKS)].pivot(index=["stratum", "arm"], columns="s",
                                    values=["sensitivity", "fdr"]))
    print(piv.round(3).to_string())


if __name__ == "__main__":
    main()
