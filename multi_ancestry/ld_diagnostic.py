#!/usr/bin/env python3
"""
LD-consistency diagnostic (SuSiE-RSS style) per population, BEFORE running MultiSuSiE.

Detects mismatch between the GWAS z-scores and the LD matrix R — the failure mode when
LD comes from a reference / subset panel (e.g. Epi25 in-sample LD used with a larger
per-pop meta-analysis) rather than the exact GWAS sample. Mismatch silently produces
confident but false credible sets.

For each population it reports, from a single eigendecomposition of R:
  s (lambda)   : susieR estimate_s_rss (null-mle), via the exact port in susier_s.py:
                 z ~ N(0, (1-s)R + sI). s≈0 → z and R consistent; larger s → more mismatch.
                 No published threshold for s exists (Zou et al. 2022; susieR issue #239).
  z_std_diff   : susieR kriging_rss's standardized difference (z - condmean)/sqrt(condvar),
                 condmean_i = -sum_{j!=i} Omega_ij z_j / Omega_ii, condvar_i = 1/Omega_ii,
                 Omega = ((1-s)R + sI)^-1 (same n adjustment of z as susieR). susieR's logLR
                 (allele-switch likelihood ratio, needs mixsqp) is NOT computed here.
  flip flag    : this script's own heuristic, not susieR's: |z_std_diff| > --zsd-thresh AND
                 observed z opposite-sign to the LD-predicted z.

History: until 2026-11-16 estimate_s_and_kriging fitted a free scale (z ~ N(0, sigma^2[(1-s)R+sI]))
and skipped the n adjustment, while being labelled susieR. That estimator is kept, renamed, as
estimate_s_scaled_legacy, only so B22/B53/B67 outputs can be reproduced; it is not susieR's s.

Outputs (in --out):
  {locus}_ld_diagnostic.tsv                 per-pop summary (one row per ancestry)
  {locus}_ld_diagnostic_snps.tsv.gz         per-SNP residuals for the worst outliers (optional)
  {locus}_{ancestry}_ld_diagnostic.png      per-ancestry kriging plot (observed vs LD-
                                             predicted z, s(lambda) annotated) [--plot]

Usage:
  python ld_diagnostic.py chr3_48883650_54214255 \\
      --harmonized-dir <harmonized_dir> \\
      --out <out_dir> [--zsd-thresh 4] [--dump-snps] [--plot]
"""
import os
import gzip
import argparse
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import norm


def estimate_s_and_kriging(z, R, n=None, r_tol=1e-8):
    """susieR estimate_s_rss + kriging_rss (z_std_diff, condmean), susieR 0.14.2 semantics.

    Returns (s, z_std_diff, cond_mean). z_std_diff and cond_mean are on susieR's scale: computed from the
    n-adjusted z when n is given, as kriging_rss does."""
    from susier_s import estimate_s_rss
    s = estimate_s_rss(z, R, n, r_tol)
    z = np.where(np.isfinite(z), np.asarray(z, dtype=np.float64), 0.0)
    lam, U = np.linalg.eigh(np.asarray(R, dtype=np.float64))
    lam = np.where(lam < r_tol, 0.0, lam)
    if n is not None:
        z = np.sqrt((n - 1) / (z ** 2 + n - 2)) * z
    v = (1.0 - s) * lam + s
    dinv = np.divide(1.0, v, out=np.zeros_like(v), where=v > 0)
    omega_z = U @ (dinv * (U.T @ z))      # Omega z
    omega_diag = (U * U) @ dinv           # diag(Omega)
    cond_mean = z - omega_z / omega_diag  # = -Omega_{i,-i} z_{-i} / Omega_ii
    z_std_diff = (z - cond_mean) * np.sqrt(omega_diag)
    return s, z_std_diff, cond_mean


def estimate_s_scaled_legacy(z, R):
    """NOT susieR. The pre-2026-11-16 estimator (free scale profiled out, no n adjustment), kept only to
    reproduce B22/B53/B67. Return (s, z_std_diff, cond_mean) using one symmetric eigendecomposition."""
    lam, U = np.linalg.eigh(R.astype(np.float64))
    lam = np.clip(lam, 0.0, None)
    b = U.T @ z
    p = z.size

    def negll(s):                       # profile -logL over scale sigma^2
        d = np.clip((1.0 - s) * lam + s, 1e-12, None)
        sigma2 = max(np.sum(b * b / d) / p, 1e-300)
        return 0.5 * (np.sum(np.log(d)) + p * np.log(sigma2))

    s = float(minimize_scalar(negll, bounds=(1e-6, 1 - 1e-6), method="bounded").x)
    d = np.clip((1.0 - s) * lam + s, 1e-12, None)
    omega_z = U @ (b / d)               # (Omega z)
    omega_diag = (U * U) @ (1.0 / d)    # diag(Omega)
    z_std_diff = omega_z / np.sqrt(omega_diag)
    cond_mean = z - omega_z / omega_diag   # LD-predicted (kriging) z
    return s, z_std_diff, cond_mean


def plot_kriging_diagnostic(z, cond_mean, zsd, flip, s, corr, anc, locus, zsd_thresh, out_png):
    """Observed z vs LD-predicted (kriging) z, per SNP -- the standard SuSiE-RSS
    consistency plot. Points that agree with their LD sit on the y=x line; the
    further a point strays, the more its z is inconsistent with the reference LD
    at that position (allele flip, strand issue, or a genuinely mismatched panel).

    Color is status, not category: most SNPs are unflagged (muted, low-emphasis);
    "outlier" (|z_std_diff| > threshold) gets the warning step; "flip candidate"
    (outlier AND opposite-signed to its LD-predicted value) gets the critical step
    -- reserved status colors, never reused as a generic series color, always
    paired with a legend label so color never carries the flag alone.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outlier = np.abs(zsd) > zsd_thresh
    normal = ~outlier
    outlier_only = outlier & ~flip

    fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
    lim = float(np.nanmax(np.abs(np.concatenate([z, cond_mean])))) * 1.08 + 1e-6
    ax.plot([-lim, lim], [-lim, lim], color="#c3c2b7", lw=1, ls="--", zorder=1,
            label="_nolegend_")

    ax.scatter(cond_mean[normal], z[normal], s=10, c="#2a78d6", alpha=0.35,
               linewidths=0, zorder=2, label=f"consistent (n={int(normal.sum())})")
    if outlier_only.any():
        ax.scatter(cond_mean[outlier_only], z[outlier_only], s=22, c="#fab219",
                   edgecolors="#0b0b0b", linewidths=0.3, zorder=3,
                   label=f"outlier, |z_std_diff|>{zsd_thresh:g} (n={int(outlier_only.sum())})")
    if flip.any():
        ax.scatter(cond_mean[flip], z[flip], s=28, c="#d03b3b",
                   edgecolors="#0b0b0b", linewidths=0.4, zorder=4,
                   label=f"flip candidate (n={int(flip.sum())})")

    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("LD-predicted z (kriging mean)")
    ax.set_ylabel("Observed GWAS z")
    ax.set_title(f"{locus} — {anc}\ns(λ)={s:.4f}   corr(z, pred)={corr:.3f}   "
                 f"n={z.size}", fontsize=11)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    return out_png


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("locus")
    ap.add_argument("--harmonized-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--zsd-thresh", type=float, default=4.0,
                    help="|z_std_diff| above this flags an outlier SNP (default 4)")
    ap.add_argument("--n", type=float, nargs="+", default=None,
                    help="GWAS sample size, one per ancestry in the file's order (or one for all); "
                         "susieR strongly recommends it")
    ap.add_argument("--dump-snps", action="store_true",
                    help="also write per-SNP residuals for outliers")
    ap.add_argument("--plot", action="store_true",
                    help="also save a per-ancestry kriging diagnostic PNG (observed vs "
                         "LD-predicted z, s(lambda) annotated)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    d = np.load(os.path.join(args.harmonized_dir, f"{args.locus}_harmonized.npz"),
                allow_pickle=True)
    ancs = [str(a) for a in d["ancestries"]]
    P = d["z"].shape[1]
    bonf = norm.isf(0.5 * 0.05 / P)     # per-SNP outlier threshold at 0.05 Bonferroni

    rows, snp_rows = [], []
    for k, anc in enumerate(ancs):
        z = d["z"][k].astype(np.float64)
        R = d["R"][k]
        n = None if args.n is None else (args.n[k] if len(args.n) > 1 else args.n[0])
        s, zsd, cm = estimate_s_and_kriging(z, R, n)
        absz = np.abs(zsd)
        n_out = int((absz > args.zsd_thresh).sum())
        n_bonf = int((absz > bonf).sum())
        # allele-flip candidates: big residual & opposite sign to LD-predicted z
        flip = (absz > args.zsd_thresh) & (np.sign(z) != np.sign(cm)) & (np.abs(cm) > 2)
        corr = float(np.corrcoef(z, cm)[0, 1]) if np.std(cm) > 0 else np.nan
        rows.append(dict(ancestry=anc, n_snps=P, s_lambda=round(s, 4),
                         corr_z_pred=round(corr, 3),
                         max_abs_zsd=round(float(absz.max()), 2),
                         n_outlier_gt_thresh=n_out, n_outlier_bonf=n_bonf,
                         n_flip_candidates=int(flip.sum())))
        print(f"  [{anc}] s(lambda)={s:.4f}  corr(z,pred)={corr:.3f}  "
              f"max|zsd|={absz.max():.1f}  outliers(>{args.zsd_thresh})={n_out}  "
              f"flip_cand={int(flip.sum())}")
        if args.plot:
            png_path = os.path.join(args.out, f"{args.locus}_{anc}_ld_diagnostic.png")
            plot_kriging_diagnostic(z, cm, zsd, flip, s, corr, anc, args.locus,
                                    args.zsd_thresh, png_path)
            print(f"    wrote {png_path}")
        if args.dump_snps:
            idx = np.where(absz > args.zsd_thresh)[0]
            for j in idx:
                snp_rows.append(dict(ancestry=anc, rsid=str(d["rsids"][j]),
                                     chr=str(d["chr"][j]), pos=int(d["pos"][j]),
                                     z=round(float(z[j]), 3),
                                     pred_z=round(float(cm[j]), 3),
                                     z_std_diff=round(float(zsd[j]), 2),
                                     flip=bool(flip[j])))

    pd.DataFrame(rows).to_csv(os.path.join(args.out, f"{args.locus}_ld_diagnostic.tsv"),
                              sep="\t", index=False)
    if args.dump_snps:
        with gzip.open(os.path.join(args.out, f"{args.locus}_ld_diagnostic_snps.tsv.gz"),
                       "wt") as g:
            pd.DataFrame(snp_rows).to_csv(g, sep="\t", index=False)
    print(f"Wrote {args.locus}_ld_diagnostic.tsv  (Bonferroni |zsd| thresh = {bonf:.2f})")


if __name__ == "__main__":
    main()
