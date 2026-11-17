#!/usr/bin/env python3
"""Supplementary: published in-sample merged-LD fine-mapping vs ours on merged reference LD.

Three published analyses of the PGC3-SCZ trans-ancestry meta-analysis with in-sample LD --
PGC-FINEMAP (Trubetskoy et al. 2022, ST11d), Mawson et al. 2026 FINEMAP and SuSiE -- scored
against the same SuSiEx benchmark as Fig. 1, next to our FINEMAP and SuSiE-RSS on the merged
reference LD. Capacity is matched to the published runs: k = 5 / L = 5. Our
FINEMAP primary run already is k = 5; our SuSiE-RSS is the v2 L = 5 sweep run, NOT the L = 10
primary of Fig. 1, which deliberately keeps each tool's default.

  a  precision-recall on the 193 loci, as Fig. 1b: positives = benchmark PIP > 0.1, calls =
     credible-set members above the arm's own PIP threshold, variants absent from the
     benchmark not scored. ONE DEPARTURE: the published tables list credible-set members
     only, so an inner join would silently drop every positive they left out of a credible
     set from the recall denominator. Here every arm, ours included, is scored against all
     389 positives, and a positive missing from an arm's table is a miss.
  b  exact lead recovery on the 122 loci, as Fig. 1c, under each arm's NATIVE credible sets:
     the published tables have no LD, so the common-purity gate of Fig. 1c cannot be applied.
     PGC-FINEMAP has no credible-set column; a region counts as one credible set when its
     expected number of causal variants is <= 1.5 (in Mawson's FINEMAP, which has the column,
     all 92 of the 122 with k <= 1.5 have one credible set, and 1 with k > 1.5 does).

Not matched, and stated in the caption: coverage (published 0.95, ours 0.99) and n (Mawson
161,405 total N; ours effective N).

Usage:
  python3 plot_supp_published_vs_ref.py --out output/supp_published_vs_ref.png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                     "axes.titleweight": "bold"})
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_sensfdr import bench_cs_loci, BENCH, BENCH_RAW
from sensitivity_fdr import load_pips, variant_key
from plot_cs_split import ARM_HUE
from plot_sensitivity_fdr import INK, INK_MUTED, GRID

HERE = os.path.dirname(os.path.abspath(__file__))
A = os.path.join(HERE, "..", "results")
ST11 = os.path.join(HERE, "..", "inputs", "external", "Trubetskoy2022_SupplementaryTable11.xlsx")
FM, SU = ARM_HUE["FINEMAP-merged"], ARM_HUE["SuSiE-merged"]
# (arm, label, colour, linestyle, hatch) -- dashed/hatched = in-sample LD
ARMS = [
    ("PGC-FINEMAP",    "FINEMAP, in-sample (PGC3)",        FM, (0, (5, 2)), "///"),
    ("Mawson-FINEMAP", "FINEMAP, in-sample (Mawson)",      FM, (0, (1, 1.4)), "..."),
    ("FINEMAP-ref",    "FINEMAP, reference (ours, k = 5)", FM, "-", None),
    ("Mawson-SuSiE",   "SuSiE-RSS, in-sample (Mawson)",    SU, (0, (5, 2)), "///"),
    ("SuSiE-ref",      "SuSiE-RSS, reference (ours, L = 5)", SU, "-", None),
]
AA_RANGE = (0.83, 0.92)   # Fig. 1c, ancestry-aware methods, common purity
# context only, in grey: the three ancestry-aware arms at their Fig. 1 (default) settings,
# scored by THIS figure's rule (389-positive denominator), so their curves are comparable here
AA_ARMS = ["SuSiEx-ours", "MultiSuSiE", "MESuSiE"]
TICK = ["FINEMAP\nPGC3", "FINEMAP\nMawson", "FINEMAP\nours", "SuSiE-RSS\nMawson", "SuSiE-RSS\nours"]


def load_arms(loci):
    """Each arm as rows (locus, key, bp, pip, cs); cs is NaN outside a credible set."""
    out = {}
    pub = load_pips(os.path.join(A, "benchmark", "published_pip_vectors.csv.gz"))
    pgc = pub[pub.method == "PGC-FINEMAP-published"].drop(columns="pip")
    st = pd.read_excel(os.path.expanduser(ST11), sheet_name="ST11d 95% Credible Sets PIP")
    st = (st.rename(columns={"position": "bp",
                             "finemap_marginal_posterior_inclusion_probability": "pip"})
            [["rsid", "bp", "pip", "expected_causals_k"]].drop_duplicates(["rsid", "bp"]))
    pgc = pgc.merge(st, on=["rsid", "bp"], how="left")
    assert pgc.pip.notna().all(), "PGC-FINEMAP rows without an ST11d PIP"
    # no CS id: one CS when expected k <= 1.5, else >= 2. Members of the >= 2 sets are unknown,
    # so every row there carries two placeholder ids; lead recovery then sees n_cs = 2
    pgc["cs"] = np.where(pgc.expected_causals_k <= 1.5, "1", "k>1.5a;k>1.5b")
    out["PGC-FINEMAP"] = pgc
    for arm, f in (("Mawson-FINEMAP", "mawson_finemap_dedup.csv"),
                   ("Mawson-SuSiE", "mawson_susie_dedup.csv")):
        d = pd.read_csv(os.path.join(A, "benchmark", f)).rename(
            columns={"locus_id": "locus_", "BP": "bp", "PIP": "pip", "CS": "cs"})
        d = d.drop(columns="locus").rename(columns={"locus_": "locus"})
        d["method"], d["a1"], d["a2"] = arm, pd.NA, pd.NA
        d["key"] = variant_key(d)          # alleles backfilled from the v2 master by rsid
        d["cs"] = d.cs.astype(str)
        out[arm] = d
    for arm, f in (("FINEMAP-ref", "main/primary_FINEMAP-merged-noambig_cs_purity.csv.gz"),
                   ("SuSiE-ref", "L_sweep/primary_SuSiE-merged-noambig_L5_cs_purity.csv.gz")):
        d = pd.read_csv(os.path.join(A, "finemap", f), low_memory=False)
        d["key"] = variant_key(d)
        d["cs"] = d.cs_ids_all
        out[arm] = d
    for arm in AA_ARMS:
        d = pd.read_csv(os.path.join(A, "finemap", "main", f"primary_{arm}-noambig_cs_purity.csv.gz"),
                        low_memory=False)
        d["key"] = variant_key(d)
        d["cs"] = d.cs_ids_all
        out[arm] = d
    return {a: d[d.locus.isin(loci)][["locus", "key", "bp", "pip", "cs"]] for a, d in out.items()}


def pr_curves(b, arms, label_pip=0.1):
    grid = np.round(np.arange(0.01, 1.0, 0.01), 2)
    pos = b[b.pb > label_pip][["locus", "key"]]
    nP = len(pos)
    rows = []
    for arm, d in arms.items():
        c = d[d.cs.notna()].groupby(["locus", "key"], as_index=False).pip.max()
        u = c.merge(b, on=["locus", "key"], how="inner")      # calls on benchmark variants
        for s in grid:
            call = u[u.pip > s]
            tp = int((call.pb > label_pip).sum())
            rows.append(dict(arm=arm, s=s, n_calls=len(call), tp=tp, n_PB=nP,
                             recall=tp / nP, precision=tp / len(call) if len(call) else np.nan))
    R = pd.DataFrame(rows)
    auc = {}
    for arm, g in R.groupby("arm"):        # Fig. 1b's convention: precision 0 beyond reach
        g = g.dropna(subset=["precision"]).sort_values("recall")
        x, y = g.recall.to_numpy(), g.precision.to_numpy()
        auc[arm] = float(np.trapezoid(np.r_[y[0], y, 0.0], np.r_[0.0, x, x[-1]])) if len(x) > 1 else np.nan
    R["auprc"] = R.arm.map(auc)
    return R, auc, nP


def lead_recovery(arms, loci122):
    b = pd.read_csv(os.path.join(HERE, BENCH))
    b = b[b.locus.isin(loci122)]
    lead = b.loc[b.groupby("locus").pip.idxmax()].set_index("locus").bp.astype(int)
    rows, per = [], []
    for arm, d in arms.items():
        d = d[d.locus.isin(loci122) & d.cs.notna()].copy()
        d["cs"] = d.cs.astype(str).str.split(r"[;,|]")
        d = d.explode("cs")
        n_ex = 0
        for l in loci122:
            g = d[d.locus == l]
            ncs = g.cs.nunique()
            ok = ncs == 1 and int(g.loc[g.pip.idxmax()].bp) == int(lead[l])
            n_ex += ok
            per.append(dict(arm=arm, locus=l, n_cs=ncs, exact=ok))
        rows.append(dict(arm=arm, n_loci=len(loci122), exact=n_ex, rate=n_ex / len(loci122)))
    return pd.DataFrame(rows), pd.DataFrame(per)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    full, _ = bench_cs_loci(BENCH_RAW)
    loci193 = set(full)
    loci122 = sorted(loci193 & {l.strip() for l in open(os.path.join(HERE, "..", "inputs", "loci_confident122.txt")) if l.strip()})
    b = load_pips(os.path.join(HERE, BENCH))
    b = b[b.locus.isin(loci193)].groupby(["locus", "key"], as_index=False).pip.max().rename(columns={"pip": "pb"})
    arms = load_arms(loci193)
    R, auc, nP = pr_curves(b, arms)
    L, P = lead_recovery({a: d for a, d in arms.items() if a not in AA_ARMS}, loci122)
    L["auprc"] = L.arm.map(auc)
    print(f"positives: {nP}")
    print(L.to_string(index=False))
    print("ancestry-aware AUPRC, this rule:", {a: round(auc[a], 3) for a in AA_ARMS})

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(15.5, 6.0), dpi=300,
                                 gridspec_kw=dict(width_ratios=[1.15, 1]))
    fig.subplots_adjust(wspace=0.28, bottom=0.33, top=0.90)
    for arm in AA_ARMS:
        g = R[(R.arm == arm) & R.precision.notna()].sort_values("s")
        ax.plot(g.recall, g.precision, color="#b9b8b1", lw=1.4, zorder=1)
    aa = [auc[a] for a in AA_ARMS]
    for arm, lab, col, ls, _ in ARMS:
        g = R[(R.arm == arm) & R.precision.notna()].sort_values("s")
        ax.plot(g.recall, g.precision, color=col, ls=ls, lw=2.0, solid_capstyle="round")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("Recall", fontsize=14, color=INK); ax.set_ylabel("Precision", fontsize=14, color=INK)
    ax.set_title(f"193 loci, {nP} benchmark variants at PIP > 0.1", fontsize=13, color=INK, loc="left")

    x = np.arange(len(ARMS))
    bx.axhspan(*AA_RANGE, color="#d9d8d2", alpha=0.6, lw=0, zorder=0)
    bx.annotate("ancestry-aware, reference LD (Fig. 1c)", (len(ARMS) - 0.5, AA_RANGE[1]),
                xytext=(0, 3), textcoords="offset points", ha="right", va="bottom",
                fontsize=10.5, color=INK_MUTED)
    for i, (arm, lab, col, ls, hatch) in enumerate(ARMS):
        v = float(L.set_index("arm").at[arm, "rate"])
        bx.bar(i, v, width=0.66, color=col if hatch is None else "white", edgecolor=col,
               hatch=hatch, linewidth=2.0, zorder=3)
        bx.annotate(f"{v:.0%}", (i, v), xytext=(0, 4), textcoords="offset points",
                    ha="center", va="bottom", fontsize=12.5, color=INK)
    bx.set_xticks(x, TICK)
    bx.set_ylim(0, 1.05)
    bx.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    bx.set_ylabel("Exact recovery of the benchmark top variant", fontsize=14, color=INK)
    bx.set_title("122 loci, single benchmark credible set", fontsize=13, color=INK, loc="left")
    for a_ in (ax, bx):
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            a_.spines[sp].set_color(GRID)
        a_.tick_params(colors=INK_MUTED, labelsize=12.5)
    bx.tick_params(axis="x", labelsize=10.5)   # "SuSiE-RSS" twice; 12.5 overlaps
    h = [mlines.Line2D([], [], color=col, ls=ls, lw=2.2, label=f"{lab}   AUPRC {auc[arm]:.2f}")
         for arm, lab, col, ls, _ in ARMS]
    h.append(mlines.Line2D([], [], color="#b9b8b1", lw=1.6,
                           label="SuSiEx, MultiSuSiE, MESuSiE: ancestry-aware, reference (Fig. 1)"))
    # no AUPRC on this entry: under the 389 denominator MESuSiE drops to ~0.77 (25 positives
    # lie outside its EUR/EAS intersection input), which would contradict Fig. 1b's 0.82
    fig.legend(handles=h, frameon=False, fontsize=12, labelcolor=INK, ncol=2,
               loc="lower center", bbox_to_anchor=(0.5, 0.01),
               title="Dashed / hatched = in-sample merged LD (published); solid = merged reference LD (ours)",
               title_fontsize=11.5)
    for lab, a_ in (("a", ax), ("b", bx)):
        a_.text(-0.12, 1.07, lab, transform=a_.transAxes, fontsize=22, fontweight="bold", color=INK)
    fig.savefig(args.out, bbox_inches="tight")
    stem = os.path.splitext(args.out)[0]
    R.to_csv(stem + "_curve.csv", index=False)
    L.to_csv(stem + "_summary.csv", index=False)
    P.to_csv(stem + "_per_locus.csv", index=False)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
