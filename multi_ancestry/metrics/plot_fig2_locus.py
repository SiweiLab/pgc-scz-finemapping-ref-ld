#!/usr/bin/env python3
"""Fig. 2A-2C: one locus, stacked panels, from Task B32's exported statistics.

Panel order, top to bottom, as specified:
  1-3.  -log10 p in EUR / EAS / MERGED, each coloured by r2 to the benchmark lead IN THAT
        SAME PANEL. Three association plots of one locus, so the reader sees the evidence
        each arm actually had rather than a single EUR view standing in for all of them.
  4.    PIP: benchmark as stems, our seven arms as a median line with a min-max ribbon.
  5.    gene track, Ensembl GRCh37 via fetch_gene_track.py (this pipeline's bp are GRCh37;
        the first, hand-pulled tables came from Ensembl's GRCh38 default and put chr12's
        lead in PPP1CC instead of ATP2A2, 372 kb away).

EUR-only arms are excluded: they are out of the primary analysis. The benchmark's variants
and their presence in each arm used to be a fifth panel; it was dropped because showing only
the top ~30 of them was more confusing than informative.

COLOUR: r2 uses a sequential single-hue ramp, light to dark, NOT LocusZoom's binned rainbow.
A rainbow is not ordered perceptually and fails for colourblind readers; the convention is
familiar but wrong, and the caption says so rather than quietly following it. Missing r2 --
a variant absent from that panel's LD -- is drawn as a hollow grey ring, never as r2 = 0,
because "not in this panel" and "unlinked" are different statements and the difference is the
substance of these loci.

Usage:
  python3 plot_fig2_locus.py --locus chr1_243512441_244077592 --out figures/fig2c.png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# every label bold (Qin): ticks and legends inherit font.weight
plt.rcParams.update({"font.weight": "bold", "axes.labelweight": "bold",
                     "axes.titleweight": "bold"})
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.ticker
import matplotlib.transforms
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_sensfdr import BENCH
from plot_cs_split import ARM_HUE, RAN_ORDER, EUR_ONLY, disp
from plot_sensitivity_fdr import INK, INK_MUTED, GRID

ABSENT = "#D8D8D8"
BENCH_C = "#3A3A3A"
RIBBON = "#8C8C8C"
FOCUS = "#555555"   # not red: CARMA's line is red in the PIP panel
RAMP = "Blues"   # single hue, light->dark: r2~0 must recede, not shout
PANELS = [("EUR", "z_eur", "p_eur", "r2_eur"),
          ("EAS", "z_eas", "p_eas", "r2_eas"),
          ("merged", "z_merged", "p_merged", "r2_merged")]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--locus", required=True)
    ap.add_argument("--b32", default=os.path.join(here, "..", "aou_ld", "task_b64_out"))
    ap.add_argument("--export", default=os.path.join(here, "..", "aou_ld", "task_v2_out", "noambig"))
    ap.add_argument("--title", default="")
    ap.add_argument("--suffix", default="-noambig", help="v2 label suffix, e.g. --suffix=-noambig")
    ap.add_argument("--focus", default="", help="comma-separated rsids to mark in every panel (e.g. the false positive)")
    ap.add_argument("--pip-dots", type=float, default=None,
                    help="also draw each arm's PIPs above this value as dots in the arm's hue (e.g. 0.1)")
    ap.add_argument("--width", type=float, default=6.0,
                    help="figure width in inches; below 9 the layout goes compact (one column of a "
                         "four-panel row): legends move under the plot, fewer ticks, shorter titles")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    L = args.locus

    assoc = pd.read_csv(os.path.join(args.b32, f"fig2_{L}_assoc.csv.gz"))
    r2 = pd.read_csv(os.path.join(args.b32, f"fig2_{L}_r2_to_lead.csv.gz"))
    d = assoc.merge(r2.drop(columns=["rsid"]), on="bp", how="left")
    lead_rsid = r2.lead_rsid.iloc[0]
    lead_bp = int(d.loc[d.rsid == lead_rsid, "bp"].iloc[0])
    focus = [x for x in args.focus.split(",") if x]
    focus_bp = {f: int(d.loc[d.rsid == f, "bp"].iloc[0]) for f in focus if (d.rsid == f).any()}

    bench = pd.read_csv(BENCH)
    bl = bench[bench.locus == L]
    bench_pip = dict(zip(bl.bp.astype(int), bl.pip.astype(float)))
    # EUR-only arms are out of the primary analysis (Qin, 2026-09-28)
    arms = [a for a in RAN_ORDER
            if a not in EUR_ONLY
            and os.path.exists(os.path.join(args.export, f"primary_{a}{args.suffix}_pip.csv.gz"))]
    grids, seen = {}, {}
    for a in arms:
        g = pd.read_csv(os.path.join(args.export, f"primary_{a}{args.suffix}_pip.csv.gz"),
                        low_memory=False)
        g = g[g.locus == L]
        grids[a] = dict(zip(g.bp.astype(int), g.pip.astype(float)))
        seen[a] = set(g.bp.astype(int))

    lo, hi = int(L.split("_")[1]), int(L.split("_")[2])
    chrom = L.split("_")[0]
    gpath = os.path.join(args.b32, f"fig2_{L}_genes.csv")
    genes = pd.read_csv(gpath) if os.path.exists(gpath) else None
    C = args.width < 9
    # gene rows are packed before the layout so the track's height follows its row count
    g2 = None
    if genes is not None and len(genes):
        # unnamed features are dropped upstream by fetch_gene_track.py, but a stale CSV
        # may still carry raw ENSG accessions, which tell a reader nothing
        genes = genes[~genes.gene.astype(str).str.fullmatch(r"ENSG\d+")]
        if "flank" not in genes.columns:
            genes = genes.assign(flank=0)
        if C:   # a quarter-width column cannot carry non-coding labels too; keep the lead's gene
            genes = genes[(genes.biotype == "protein_coding")
                          | ((genes.start <= lead_bp) & (genes.end >= lead_bp))]
        g2 = genes.sort_values(["flank", "start"]).reset_index(drop=True)
        # pack rows by what each gene OCCUPIES, bar or centred label, whichever is wider, so
        # labels never run into a neighbour's -- matters most when the figure is narrow
        bp_per_in = (hi - lo) / (args.width * 0.72)
        ends, rows = [], []
        for _, g in g2.iterrows():
            labw = (len(str(g.gene)) + 2) * (7.6 if C else 9.5) * 0.6 / 72 * bp_per_in
            mid = (max(g.start, lo) + min(g.end, hi)) / 2
            a0, a1 = min(max(g.start, lo), mid - labw / 2), max(min(g.end, hi), mid + labw / 2)
            r = next((k for k, e in enumerate(ends) if a0 > e + 0.01 * (hi - lo)), len(ends))
            if r == len(ends):
                ends.append(-np.inf)
            ends[r] = a1
            rows.append(r)
        g2["row"] = rows
    n_grow = int(g2.row.max()) + 1 if g2 is not None and len(g2) else 1
    fig, axes = plt.subplots(5, 1, figsize=(args.width, 12.2), dpi=300,
                             gridspec_kw=dict(height_ratios=[1, 1, 1, 1.1, float(np.clip(0.17 * (n_grow + 1), 0.3, 1.1))]))
    fig.subplots_adjust(hspace=0.42, top=0.95, bottom=0.08)

    # ---- 1-3: association, one per LD panel ------------------------------------------
    for k, (name, zc, pc, rc) in enumerate(PANELS):
        ax = axes[k]
        m = d[pc].notna()
        y = -np.log10(d.loc[m, pc].clip(lower=1e-300))
        rr = d.loc[m, rc]
        has = rr.notna()
        ax.scatter(d.loc[m, "bp"][~has], y[~has], s=11, facecolors="none",
                   edgecolors="#B0B0B0", linewidths=0.6, zorder=2)
        sc = ax.scatter(d.loc[m, "bp"][has], y[has], c=rr[has], cmap=RAMP, vmin=-0.12, vmax=1,
                        s=13, linewidths=0, zorder=3)
        if (d.bp == lead_bp).any() and pd.notna(d.loc[d.bp == lead_bp, pc].iloc[0]):
            ly = -np.log10(max(float(d.loc[d.bp == lead_bp, pc].iloc[0]), 1e-300))
            ax.scatter([lead_bp], [ly], marker="D", s=46, facecolor="none",
                       edgecolor=INK, linewidths=1.3, zorder=6)
        ax.axhline(-np.log10(5e-8), color=GRID, lw=0.9, ls="--", zorder=1)
        for fr in focus:  # the variant(s) the panel is about, other than the benchmark lead
            fm = (d.rsid == fr) & d[pc].notna()
            if fm.any():
                fy = -np.log10(max(float(d.loc[fm, pc].iloc[0]), 1e-300))
                ax.scatter(d.loc[fm, "bp"], [fy], marker="s", s=46, facecolor="none", edgecolor=FOCUS,
                           linewidths=1.4, zorder=7)
        n_miss = int((~has).sum())
        ax.set_title(f"{name}: {int(has.sum()):,} with LD" + (f", {n_miss:,} without" if n_miss else "")
                     if C else
                     f"{name} — {int(has.sum()):,} variants with LD"
                     + (f", {n_miss:,} without (hollow)" if n_miss else ""),
                     fontsize=11.0, color=INK, loc="left", pad=4)
        ax.set_ylabel("$-\\log_{10}p$", fontsize=10.9, color=INK)
        cb = fig.colorbar(sc, ax=ax, fraction=0.016, pad=0.008)
        cb.set_label("r$^2$ to top variant" if C else f"r$^2$ to {lead_rsid}", fontsize=9.0, color=INK_MUTED)
        cb.ax.tick_params(colors=INK_MUTED, labelsize=8.3)
        cb.outline.set_visible(False)

    # ---- 5: PIP, benchmark vs the trans-ancestry arms ------------------------------------------
    axP = axes[3]
    # median line + min-max ribbon over the arms (Qin: the thin per-arm lines read as clutter).
    # Neutral grey, so it never reads as one of the arm hues the dots below use.
    allbp = sorted({bp for a in arms for bp in grids[a]})
    A = np.array([[grids[a].get(bp, np.nan) for bp in allbp] for a in arms])
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med, mn, mx = (np.nanmedian(A, axis=0), np.nanmin(A, axis=0), np.nanmax(A, axis=0))
    x = np.array(allbp)
    axP.fill_between(x, np.nan_to_num(mn), np.nan_to_num(mx), color=RIBBON, alpha=0.45, linewidth=0,
                     zorder=2, label=f"min–max, {len(arms)} methods")
    axP.plot(x, np.nan_to_num(med), color=RIBBON, lw=1.1, zorder=3, label=f"median, {len(arms)} methods")
    if args.pip_dots is not None:
        # every arm's PIP above the threshold as a dot in the arm's hue, so who called what is legible
        for a in arms:
            hx = np.array(sorted(b for b, v in grids[a].items() if v > args.pip_dots))
            if len(hx):
                axP.scatter(hx, [grids[a][b] for b in hx], s=24, color=ARM_HUE[a], edgecolor="white",
                            linewidths=0.5, zorder=5, label=disp(a))
    bx = np.array(sorted(bench_pip)); bv = np.array([bench_pip[b] for b in bx])
    mb = bv >= 0.01
    axP.scatter(bx[mb], bv[mb], s=18, color=BENCH_C, zorder=6, label="benchmark")
    for k, (f, fb) in enumerate(sorted(focus_bp.items(), key=lambda x: x[1])):
        for ax in axes[:4]:
            ax.axvline(fb, color=FOCUS, lw=0.8, ls=":", zorder=1)
        # alternate sides of the line so neighbouring focus labels never overlap
        axes[0].annotate(f, (fb, 1.0), xycoords=("data", "axes fraction"), xytext=(-3 if k % 2 == 0 else 3, -2),
                         textcoords="offset points", ha="right" if k % 2 == 0 else "left", va="top", fontsize=8.5,
                         color=FOCUS, rotation=90)
    axP.set_ylim(0, 1.05)
    axP.set_ylabel("PIP", fontsize=10.9, color=INK)
    # outside the axes, on the right where panels 1-3 carry their colourbars, so the
    # legend never sits on top of the data it describes
    if not C:   # compact: these keys join the figure legend under the plot instead
        axP.legend(frameon=False, fontsize=9.0, labelcolor=INK, loc="upper left",
                   bbox_to_anchor=(1.005, 1.0), borderaxespad=0)

    # ---- gene track (Ensembl GRCh37, see fetch_gene_track.py) -----------------------------------
    axG = axes[4]
    if g2 is not None and len(g2):
        span = hi - lo
        # the gene the lead variant falls in, highlighted as in the earlier figures
        g2["hit"] = (g2.flank == 0) & (g2.start <= lead_bp) & (g2.end >= lead_bp)
        for _, g in g2.iterrows():
            coding = g.biotype == "protein_coding"
            x0, x1 = max(g.start, lo), min(g.end, hi)
            if g.flank:
                # outside the window: a stub at the edge it lies beyond, labelled with how
                # far away it is, so a gene desert still tells the reader what is nearby
                right = g.start > hi
                x0, x1 = ((hi - 0.055 * span, hi) if right else (lo, lo + 0.055 * span))
                d_kb = (g.start - hi if right else lo - g.end) / 1000.0
                lab = f"{g.gene}  {d_kb:,.0f} kb"
            else:
                lab = g.gene
            col = (BENCH_C if g.hit else (INK if coding else INK_MUTED))
            axG.plot([x0, x1], [-g.row] * 2,
                     lw=6 if coding or g.hit else 2.6, color=col,   # one height: length is size (Qin)
                     solid_capstyle="butt", zorder=4 if g.hit else 3,
                     alpha=0.45 if g.flank else 1.0)
            axG.annotate(lab, ((x0 + x1) / 2, -g.row - 0.32),
                         ha="right" if (g.flank and g.start > hi) else
                            ("left" if g.flank else "center"),
                         va="top",
                         fontsize=(8.4 if g.hit else 7.4) if C else (10 if g.hit else (9 if coding else 7.6)),
                         fontweight="bold" if g.hit else "normal",
                         color=(BENCH_C if g.hit else
                                (INK_MUTED if (g.flank or not coding) else INK)),
                         style="normal" if coding else "italic")
        axG.set_ylim(-n_grow, 0.6)
        if not ((g2.biotype == "protein_coding") & (g2.flank == 0)).any():
            axG.annotate("no protein-coding gene in this interval; nearest shown at "
                         "the edge" if (g2.flank == 1).any() else
                         "no protein-coding gene in this interval",
                         (lo, 0.45), ha="left", va="top", fontsize=10.5, color=INK_MUTED)
    axG.set_yticks([])
    for sp in axG.spines.values():
        sp.set_visible(False)

    for ax in [axes[0], axes[1], axes[2], axP, axG]:
        ax.set_xlim(lo, hi)
        # enough decimals for the tick spacing: on a 50 kb window the ticks fall at x.xx5 Mb, and 2 decimals mislabelled
        # them by 5 kb (Fig. 2b, Qin)
        ax.xaxis.set_major_formatter(lambda v, _: f"{v / 1e6:.{3 if hi - lo < 300_000 else 2}f}")
        ax.grid(False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=INK_MUTED, labelsize=9.6)
        if C:
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    axG.set_xlabel(f"{chrom} position (Mb)", fontsize=12.2, color=INK)

    h = [mlines.Line2D([], [], marker="D", ls="", markerfacecolor="none",
                       markeredgecolor=INK, markersize=7, label="benchmark top variant"),
         mlines.Line2D([], [], marker="o", ls="", markerfacecolor="none",
                       markeredgecolor="#B0B0B0", markersize=6,
                       label="no r$^2$ in that panel (not r$^2$=0)"),
         ]
    if focus_bp:
        h.append(mlines.Line2D([], [], marker="s", ls=":", color=FOCUS, markerfacecolor="none",
                               markeredgecolor=FOCUS, markersize=7, label="focus: " + ", ".join(focus_bp)))
    if C:
        # two columns (Qin): every arm on the left -- all of them, whether or not it has a dot at
        # this locus, so the key is the same in 2a-2c -- and everything else on the right
        hh = dict(zip(*axP.get_legend_handles_labels()[::-1]))
        right = [hh[k] for k in [f"min–max, {len(arms)} methods", f"median, {len(arms)} methods", "benchmark"]
                 if k in hh] + h
        left = ([mlines.Line2D([], [], marker="o", ls="", color=ARM_HUE[a], markeredgecolor="white",
                               markersize=6.5, label=disp(a)) for a in arms]
                if args.pip_dots is not None else [])
        n = max(len(left), len(right))
        blank = lambda: mlines.Line2D([], [], ls="", label=" ")
        h = left + [blank() for _ in range(n - len(left))] + right + [blank() for _ in range(n - len(right))]
    # anchored to the gene track, just under its x label, so it never floats far from the plot
    fig.legend(handles=h, frameon=False, fontsize=9.7, ncol=2 if C else 3, labelcolor=INK,
               loc="upper center", bbox_to_anchor=(0.5, 0.0),
               # a fixed 0.62 in under the gene axis (clears its ticks and x label) whatever its height
               bbox_transform=matplotlib.transforms.offset_copy(axG.transAxes, fig=fig, y=-0.62, units="inches"))
    fig.suptitle(f"{chrom}:{lo:,}–{hi:,}" + ("\n" if C else "   ·   ") + f"top variant {lead_rsid}"
                 + (f"   ·   {args.title}" if args.title else ""),
                 fontsize=13 if C else 15.4, color=INK, x=0.004, ha="left", y=1.005 if C else 0.997)
    note = (
        "Association statistics and r2 on the non-palindromic variants of the harmonised "
        "master, with r2 from the LD panels the arms were run on (EUR: GCTB 10k; EAS: AoU; "
        "MERGED: the per-locus weighted mix of the two used by the merged arms). Panels top to bottom: "
        "-log10 p in EUR, EAS and MERGED (trans-ancestry meta-analysis), each coloured by r2 "
        "to the benchmark lead in that same panel; and PIP, with the benchmark (Yuan et al.) "
        f"as black points and our {len(arms)} trans-ancestry arms as a grey median line with a "
        "min-max ribbon"
        + (f"; every PIP above {args.pip_dots:g} from any arm is a dot in that arm's colour"
           if args.pip_dots is not None else "")
        + (". Open squares and dotted lines mark " + ", ".join(focus_bp)
           if focus_bp else "")
        + ". r2 uses a "
        "sequential single-hue ramp rather than LocusZoom's binned rainbow: the convention "
        "is familiar but a rainbow is not perceptually ordered and fails for colourblind "
        "readers. A HOLLOW RING is a variant with no r2 in that panel, never r2 = 0 -- 'not "
        "in this panel' and 'unlinked' are different statements and the difference is the "
        "substance of these loci. The dashed line is p = 5e-8. The GENE TRACK is Ensembl "
        "GRCh37, matching this pipeline's coordinates, fetched by "
        "metrics/fetch_gene_track.py, which refuses to write unless the locus lead's "
        "position on that assembly equals the bp our own table carries. Bars are "
        "protein-coding genes, their length the gene's span; non-coding genes, clone-name "
        "placeholders and pseudogenes are omitted. The gene containing the benchmark lead is "
        "labelled in bold. A gene shown "
        "faded at an edge with a distance lies outside the interval and is the nearest "
        "protein-coding gene to a locus that contains none.")
    cp = args.out.rsplit(".", 1)[0] + "_caption.txt"
    with open(cp, "w") as fh:
        fh.write(" ".join(note.split()) + "\n")
    # the colourbars narrow panels 1-3 only; match the PIP and gene panels to them so one
    # bp sits at the same x in every panel (the focus lines must cross the stems they mark)
    p0 = axes[0].get_position()
    for ax in (axP, axG):
        p = ax.get_position()
        ax.set_position([p0.x0, p.y0, p0.width, p.height])
    fig.savefig(args.out, bbox_inches="tight", facecolor="white")
    print(f"[done] wrote {args.out}")


if __name__ == "__main__":
    main()
