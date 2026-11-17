#!/usr/bin/env python3
"""Supp Fig 5: one variable at a time -- published in-sample analyses as the benchmark.

Exact recovery of each benchmark's top variant on its single-signal loci (one credible set,
top PIP > 0.1), excluding loci whose benchmark top variant is strand-ambiguous: the published
merged-LD analyses keep A/T and C/G SNPs, ours drop them, as the SuSiEx benchmark does. Our merged-LD runs are scored at the benchmark's credible-set coverage: 0.99 against the
SuSiEx benchmark, 0.95 against the published merged-LD analyses. Every analysis uses its own (native) credible sets, as in Supp Fig 3b: the
published results carry no LD, so the common-purity gate cannot be applied to them, and
gating only our side would add a second difference to a same-model comparison.

  a  PGC3 FINEMAP (Trubetskoy 2022, ST11d) as the benchmark; hatched = in-sample LD.
  b  every published analysis as the benchmark (rows) against every analysis (columns).

Usage:
  python3 plot_supp_published_benchmark.py --out output/supp_published_benchmark.png
"""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold"})
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_supp_published_vs_ref as P
from bench_sensfdr import BENCH, BENCH_RAW, bench_cs_loci
from sensitivity_fdr import load_pips
from plot_cs_split import ARM_HUE
from plot_sensitivity_fdr import INK, INK_MUTED, GRID

OURS = {"SuSiEx-ours": "main/primary_SuSiEx-ours-noambig_cs_purity.csv.gz",
        "MultiSuSiE":  "main/primary_MultiSuSiE-noambig_cs_purity.csv.gz",
        "MESuSiE":     "main/primary_MESuSiE-noambig_cs_purity.csv.gz",
        "FINEMAP-ref": "main/primary_FINEMAP-merged-noambig_cs_purity.csv.gz",
        "SuSiE-ref":   "L_sweep/primary_SuSiE-merged-noambig_L5_cs_purity.csv.gz",
        # coverage 0.95: scored against the published merged-LD analyses, which use 95% sets
        "FINEMAP-ref@0.95": "coverage_95/primary_FINEMAP-merged-noambig_cov095_cs_purity.csv.gz",
        "SuSiE-ref@0.95":   "coverage_95/primary_SuSiE-merged-noambig_L5_cov095_cs_purity.csv.gz"}
# credible-set coverage of each benchmark; our merged-LD runs are scored at the benchmark's
COVERAGE_OF = {"SuSiEx-published": 0.99, "PGC-FINEMAP": 0.95, "Mawson-FINEMAP": 0.95, "Mawson-SuSiE": 0.95}
# (arm, tick label, colour, hatch) -- hatched = in-sample LD, solid = reference LD
COLS = [("SuSiEx-published", "SuSiEx\nbenchmark", ARM_HUE["SuSiEx-ours"], "///"),
        ("SuSiEx-ours", "SuSiEx\nours", ARM_HUE["SuSiEx-ours"], None),
        ("MultiSuSiE", "MultiSuSiE\nours", ARM_HUE["MultiSuSiE"], None),
        ("MESuSiE", "MESuSiE\nours", ARM_HUE["MESuSiE"], None),
        ("PGC-FINEMAP", "FINEMAP\nPGC3", ARM_HUE["FINEMAP-merged"], "///"),
        ("Mawson-FINEMAP", "FINEMAP\nMawson", ARM_HUE["FINEMAP-merged"], "..."),
        ("FINEMAP-ref", "FINEMAP\nours", ARM_HUE["FINEMAP-merged"], None),
        ("Mawson-SuSiE", "SuSiE-RSS\nMawson", ARM_HUE["SuSiE-merged"], "///"),
        ("SuSiE-ref", "SuSiE-RSS\nours", ARM_HUE["SuSiE-merged"], None)]
ROWS = [("SuSiEx-published", "SuSiEx benchmark"), ("PGC-FINEMAP", "FINEMAP, PGC3"),
        ("Mawson-FINEMAP", "FINEMAP, Mawson"), ("Mawson-SuSiE", "SuSiE-RSS, Mawson")]
N_AA = 4   # first four columns are ancestry-aware


def load():
    braw = load_pips(os.path.join(P.HERE, BENCH))
    loci = set(braw.locus)
    arms = P.load_arms(loci)
    _, mem = bench_cs_loci(BENCH_RAW)
    bb = (braw.groupby(["locus", "bp"], as_index=False).pip.max()
              .merge(mem.drop_duplicates(["locus", "bp"]), on=["locus", "bp"], how="left"))
    arms["SuSiEx-published"] = bb.assign(cs=bb.cs.str.replace("|", "_", regex=False))
    X = os.path.join(P.A, "finemap")
    for a, f in OURS.items():
        d = pd.read_csv(os.path.join(X, f), low_memory=False)
        d = d[d.locus.isin(loci) & d.cs_ids_all.notna()]
        arms[a] = d.assign(cs=d.cs_ids_all)[["locus", "bp", "pip", "cs"]]
    S = {}
    for a, d in arms.items():
        d = d[d.cs.notna()].copy()
        d["cs"] = d.cs.astype(str).str.split(r"[;,|]")
        d = d.explode("cs").reset_index(drop=True)
        d["cs"] = d.cs.str.strip()
        g = d.groupby("locus")
        top = d.loc[g.pip.idxmax()].set_index("locus")
        S[a] = pd.DataFrame({"ncs": g.cs.nunique(), "top_bp": top.bp.astype(int), "top_pip": top.pip})
    # strand-ambiguous (A/T, C/G) variants, from PGC3 ST11d's alleles: all three published merged-LD
    # analyses come from the same PGC3 data, and every benchmark top our inputs lack is in ST11d
    st = pd.read_excel(P.ST11, sheet_name="ST11d 95% Credible Sets PIP")
    st = st[st.chromosome.astype(str) != "X"]
    amb = {(int(c), int(b)) for c, b, a1, a2 in zip(st.chromosome, st.position, st.GWAS_effect_allele,
                                                     st.other_allele)
           if {str(a1).upper(), str(a2).upper()} in ({"A", "T"}, {"C", "G"})}
    return S, amb


def score(S, amb):
    rate, cnt, nloc = {}, {}, {}
    for r, _ in ROWS:
        R = S[r]; ss = R[(R.ncs == 1) & (R.top_pip > 0.1)]
        # one variable at a time: drop loci whose benchmark top variant is strand-ambiguous --
        # kept by the published merged-LD analyses, removed from ours as from the SuSiEx benchmark.
        # A top missing from our input for any other reason stays a miss, as in Fig. 1c.
        ss = ss[[(int(l.split("_")[0][3:]), int(ss.top_bp[l])) not in amb for l in ss.index]]
        nloc[r] = len(ss)
        for t, *_ in COLS:
            if t == r:
                continue
            src = f"{t}@0.95" if COVERAGE_OF[r] == 0.95 and f"{t}@0.95" in S else t
            T = S[src].reindex(ss.index)
            k = int(((T.ncs == 1) & (T.top_bp == ss.top_bp)).sum())
            rate[(r, t)], cnt[(r, t)] = k / len(ss), k
    return rate, cnt, nloc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    rate, cnt, nloc = score(*load())

    fig = plt.figure(figsize=(15.5, 12.4), dpi=300)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.78], hspace=0.42)
    ax, hx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    # a -- PGC3 FINEMAP as the benchmark
    bars = [c for c in COLS if c[0] != "PGC-FINEMAP"]
    for i, (a, lab, col, hatch) in enumerate(bars):
        v = rate[("PGC-FINEMAP", a)]
        ax.bar(i, v, width=0.66, color=col if hatch is None else "white", edgecolor=col,
               hatch=hatch, linewidth=2.0, zorder=3)
        ax.annotate(f"{v:.0%}", (i, v), xytext=(0, 4), textcoords="offset points",
                    ha="center", va="bottom", fontsize=12.5, color=INK)
    ax.axvline(N_AA - 0.5, color=GRID, lw=1.3, zorder=1)
    for x, t in (((N_AA - 1) / 2, "ancestry-aware model"),
                 ((N_AA + len(bars) - 1) / 2, "merged model (as the benchmark)")):
        ax.annotate(t, (x, -0.17), xycoords=("data", "axes fraction"), ha="center",
                    va="top", fontsize=12.5, color=INK)
    ax.set_xticks(range(len(bars)), [b[1] for b in bars], fontsize=11.5)
    ax.set_xlim(-0.6, len(bars) - 0.4); ax.set_ylim(0, 1.05)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    ax.set_ylabel("Exact recovery of the\nbenchmark top variant", fontsize=14, color=INK)
    ax.set_title(f"Benchmark: PGC3 FINEMAP, in-sample merged LD; {nloc['PGC-FINEMAP']} loci, "
                 "one credible set, top PIP > 0.1, top not strand-ambiguous", fontsize=13, color=INK, loc="left")

    # b -- every published analysis as the benchmark
    M = np.full((len(ROWS), len(COLS)), np.nan)
    for i, (r, _) in enumerate(ROWS):
        for j, (t, *_) in enumerate(COLS):
            M[i, j] = rate.get((r, t), np.nan)
    cmap = LinearSegmentedColormap.from_list("seq", ["#f2f6fb", "#9dbbe0", "#2a5c9e"])
    cmap.set_bad("#e9e8e4")
    hx.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0.3, vmax=1.0, aspect="auto")
    for i in range(len(ROWS)):
        for j in range(len(COLS)):
            v = M[i, j]
            txt, c = ("–", INK_MUTED) if np.isnan(v) else (f"{v:.0%}", "white" if v > 0.75 else INK)
            hx.text(j, i, txt, ha="center", va="center", fontsize=12.5, color=c)
    hx.axvline(N_AA - 0.5, color="white", lw=4)
    hx.set_xticks(range(len(COLS)), [c[1] for c in COLS], fontsize=11.5)
    hx.set_yticks(range(len(ROWS)), [f"{lab}  (n = {nloc[r]})" for r, lab in ROWS], fontsize=12)
    hx.tick_params(length=0, colors=INK_MUTED)
    for sp in hx.spines.values():
        sp.set_visible(False)
    hx.set_ylabel("Benchmark", fontsize=14, color=INK)
    hx.set_title("Every published in-sample fine-mapping as benchmark; published and our analyses scored in columns",
                 fontsize=13, color=INK, loc="left")
    for x, t in (((N_AA - 1) / 2, "ancestry-aware model"),
                 ((N_AA + len(COLS) - 1) / 2, "merged model")):
        hx.annotate(t, (x, -0.24), xycoords=("data", "axes fraction"), ha="center",
                    va="top", fontsize=12.5, color=INK)

    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, labelsize=12.5)
    ax.text(0.995, 0.985, "Hatched = in-sample LD (published)\nSolid = reference LD (ours)",
            transform=ax.transAxes, ha="right", va="top", fontsize=11.5, color=INK_MUTED)
    for lab, a_ in (("a", ax), ("b", hx)):
        a_.text(-0.08, 1.08, lab, transform=a_.transAxes, fontsize=22, fontweight="bold", color=INK)
    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    stem = os.path.splitext(args.out)[0]
    pd.DataFrame([dict(benchmark=r, n_loci=nloc[r], analysis=t, exact=cnt[(r, t)], rate=rate[(r, t)],
                       our_coverage=(COVERAGE_OF[r] if t in ("FINEMAP-ref", "SuSiE-ref") else None))
                  for (r, t) in rate]).to_csv(stem + "_matrix.csv", index=False)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
