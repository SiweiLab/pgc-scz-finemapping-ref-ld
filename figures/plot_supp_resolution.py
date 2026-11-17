#!/usr/bin/env python3
"""Supplementary figure: fine-mapping resolution and its agreement with the benchmark, 193 loci.

Denominator: the loci the benchmark (Yuan et al.) resolved to a lead with PIP > --pip (default 0.5; 41 loci).
Each method's bar splits those loci by its own top credible-set variant: the same variant at PIP > --pip (solid), a
different variant at PIP > --pip (hatched), or not resolved (grey). v2, no-palindromic, credible-set members only (as Fig. 1b/2d).
  python3 plot_supp_resolution.py --out output/supp_resolution_pip05.png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold"})
import matplotlib.patches as mpatches
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bench_sensfdr import BENCH, bench_cs_loci, BENCH_RAW
from plot_cs_split import ARM_HUE, disp
from plot_sensitivity_fdr import INK, INK_MUTED, GRID

PER_ANC = ["SuSiEx-ours", "MultiSuSiE", "MESuSiE"]
MERGED = ["FINEMAP-merged", "SuSiE-merged", "RSparsePro", "CARMA"]
NOTCONF = "#E4E4E4"   # not resolved: neutral, so only the arm hue carries identity


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export", default=os.path.join(HERE, "..", "results", "finemap", "main"))
    ap.add_argument("--suffix", default="-noambig")
    ap.add_argument("--pip", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    full, _ = bench_cs_loci(BENCH_RAW)
    b = pd.read_csv(BENCH, dtype={"rsid": str})
    b = b[b.locus.isin(full)]
    bt = b.loc[b.groupby("locus").pip.idxmax()].set_index("locus")
    conf_b = bt[bt.pip > a.pip]                  # the benchmark's confidently resolved loci: the denominator
    n = len(conf_b)

    rows, extra = [], []
    for arm in PER_ANC + MERGED:
        d = pd.read_csv(os.path.join(a.export, f"primary_{arm}{a.suffix}_cs_purity.csv.gz"),
                        dtype={"rsid": str}, low_memory=False)
        d = d[d.locus.isin(full) & d.cs_ids_all.notna()]
        top = d.loc[d.groupby("locus").pip.idxmax()].set_index("locus")
        t_ = top.reindex(conf_b.index)
        conf = t_.pip > a.pip
        same = conf & (t_.rsid == conf_b.rsid)
        rows.append(dict(arm=arm, group="per-ancestry LD" if arm in PER_ANC else "merged LD", n=n,
                         same=int(same.sum()), different=int((conf & ~same).sum()), not_confident=int((~conf).sum())))
        # outside the denominator: the arm is confident on the benchmark's lead where the benchmark itself is not
        o = top[(top.pip > a.pip) & ~top.index.isin(conf_b.index)]
        extra.append(int((o.rsid == bt.rsid.reindex(o.index)).sum()))
    t = pd.DataFrame(rows)
    t.to_csv(a.out.rsplit(".", 1)[0] + ".csv", index=False)
    g = t.groupby("group")[["same", "n"]].sum()

    fig, ax = plt.subplots(figsize=(7.6, 4.0), dpi=300)
    ys = [7.6, 6.6, 5.6] + [4.0, 3.0, 2.0, 1.0]    # a gap row between the two families
    for y, r in zip(ys, t.itertuples()):
        hue = ARM_HUE[r.arm]
        ax.barh(y, r.same, height=0.66, color=hue, edgecolor="white", linewidth=1.2)
        # hatched in the arm's hue, so a confident disagreement never reads as the grey "not resolved"
        ax.barh(y, r.different, left=r.same, height=0.66, facecolor="white", edgecolor=hue, hatch="////",
                linewidth=0)
        ax.barh(y, r.not_confident, left=r.same + r.different, height=0.66, color=NOTCONF, edgecolor="white",
                linewidth=1.2)
        ax.annotate(f"{r.same}/{n}  ({100 * r.same / n:.0f}%)", (n + 1.2, y), va="center", ha="left",
                    fontsize=9.0, color=INK)
    for y0, grp in ((6.6, "per-ancestry LD"), (2.5, "merged LD")):
        k, m = g.loc[grp, "same"], g.loc[grp, "n"]
        ax.annotate(f"{grp}\n{k}/{m} ({100 * k / m:.0f}%)", (1.02, y0), xycoords=("axes fraction", "data"),
                    va="center", ha="left", fontsize=9.0, color=INK_MUTED)
    ax.axhline(4.8, color=GRID, lw=1.0)
    ax.set_yticks(ys)
    ax.set_yticklabels([disp(r) for r in t.arm], fontsize=9.6, color=INK)
    ax.set_xlim(0, n + 11)
    ax.spines["bottom"].set_bounds(0, n)
    ax.set_xticks([x for x in range(0, n, 10) if n - x >= 5] + [n])
    ax.set_xlabel(f"The {n} loci the benchmark resolved to a top variant with PIP > {a.pip:g}", fontsize=10.4,
                  color=INK)
    ax.grid(False)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, labelsize=9.2)
    ax.tick_params(axis="y", length=0)
    h = [mpatches.Patch(facecolor=INK_MUTED, label=f"same variant, PIP > {a.pip:g}"),
         mpatches.Patch(facecolor="white", edgecolor=INK_MUTED, hatch="////", linewidth=0,
                        label=f"different variant, PIP > {a.pip:g}"),
         mpatches.Patch(facecolor=NOTCONF, label=f"not resolved (PIP \u2264 {a.pip:g} or no CS)")]
    ax.legend(handles=h, frameon=False, fontsize=8.8, labelcolor=INK, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.17))
    fig.savefig(a.out, bbox_inches="tight", facecolor="white")

    print(t.to_string(index=False)); print("extra (outside the denominator):", dict(zip(t.arm, extra)))
    print(f"[done] wrote {a.out}")


if __name__ == "__main__":
    main()
