#!/usr/bin/env python3
"""Rebuild the published SuSiEx posterior in full, from Kai Yuan's per-locus output.

Why this exists: Suppl. Table 27 is doubly filtered -- variants with PIP>0.1 *and* in a
purity-surviving credible set -- which left 389 variants over 136 loci and made the
metric's FDR term uninterpretable (every counted false discovery was an absence, because
"listed but below the noise floor" cannot exist in a table truncated at 0.1).

Kai shared the run's own per-locus output (2026-09-11). The `.snp` files carry per-
population log Bayes factors per credible set, from which the single-effect alpha vectors
follow exactly:

    alpha_kj = softmax_j( LogBF(CSk,Pop1) + LogBF(CSk,Pop2) )      [uniform prior]
    PIP_j    = 1 - prod_k (1 - alpha_kj)

Validated against the run's own `.cs` file: reproduces SuSiEx's CS_PIP to ~1e-5, which is
the printing precision of the LogBF columns. The reconstruction also reproduces every
published headline count exactly -- 193 loci with a credible set, 215 credible sets, 136
loci with a variant above PIP 0.1 (= ST27's locus count), 11 above 0.95.

What it still cannot recover: the 57 loci where no credible set survived emit no LogBF
column at all, and a component dropped by the purity gate leaves no trace in the file. So
this is the full posterior *over surviving components*, not over the unfiltered fit.

Usage:
  python3 reconstruct_published_full_posterior.py --src ../inputs/202.EUR_ASN.SuSiEx_spec \
      --allele-source ../aou_ld/task_a_out/pip_vectors_full250_qcexcl.csv.gz \
      --out pip_vectors_susiex_published_full.csv.gz
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

METHOD = "SuSiEx-published-full"

# Yuan et al. trimmed PGC's chr15 window to the observed variant extent (343 bp off the
# start, 668 off the end). Verified 2026-09-11: our 384 rows for this locus span exactly
# Kai's boundaries, so the two ids name the same locus with the same variants. Aliased
# rather than renamed -- renaming would touch nine arms' committed vectors for no gain.
LOCUS_ALIAS = {"chr15_34646330_34663925": "chr15_34645987_34664593"}


def locus_id_from_summary(path):
    first = open(path).readline().strip()      # "# chr7:1827867-2344867"
    m = re.match(r"#\s*chr(\w+):(\d+)-(\d+)", first)
    if not m:
        sys.exit(f"[error] cannot parse region from {path}: {first!r}")
    lid = f"chr{m.group(1)}_{m.group(2)}_{m.group(3)}"
    return LOCUS_ALIAS.get(lid, lid)


def posterior_from_snp(path):
    """Full per-variant PIP for one locus. Empty frame if no credible set survived."""
    d = pd.read_csv(path, sep="\t")
    cs_ids = sorted({int(m.group(1)) for c in d.columns
                     if (m := re.match(r"LogBF\(CS(\d+),Pop\d+\)", c))})
    if not cs_ids:
        return pd.DataFrame({"rsid": d["SNP"], "bp": d["Pos"], "pip": 0.0})
    alphas = []
    for k in cs_ids:
        pops = [c for c in d.columns if c.startswith(f"LogBF(CS{k},Pop")]
        tot = d[pops].sum(axis=1).to_numpy(dtype=float)
        a = np.exp(tot - tot.max())
        alphas.append(a / a.sum())
    pip = 1.0 - np.prod(1.0 - np.vstack(alphas), axis=0)
    return pd.DataFrame({"rsid": d["SNP"], "bp": d["Pos"], "pip": pip})


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", required=True, help="directory of SuSiEx.*.cs99.{snp,cs,summary}")
    ap.add_argument("--allele-source", required=True,
                    help="our own PIP vectors, to backfill a1/a2 by (locus, rsid)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    src = os.path.expanduser(args.src)
    snps = sorted(glob.glob(os.path.join(src, "*.cs99.snp")))
    if not snps:
        sys.exit(f"[error] no *.cs99.snp under {src}")

    frames = []
    for snp_path in snps:
        stem = snp_path[:-4]
        lid = locus_id_from_summary(stem + ".summary")
        g = posterior_from_snp(snp_path)
        g.insert(0, "locus", lid)
        frames.append(g)
    d = pd.concat(frames, ignore_index=True)
    print(f"[info] {len(snps)} loci, {len(d):,} variants reconstructed")

    # Alleles: the .snp carries none, and the merged-LD arms key on (bp, allele pair),
    # so without them the join silently degrades to rsid-only for half the comparison.
    src_pips = pd.read_csv(os.path.expanduser(args.allele_source),
                           usecols=["locus", "rsid", "a1", "a2"]).dropna(subset=["a1", "a2"])
    src_pips = src_pips.drop_duplicates(["locus", "rsid"])
    d = d.merge(src_pips, on=["locus", "rsid"], how="left")
    filled = d["a1"].notna().sum()
    print(f"[info] alleles backfilled for {filled:,}/{len(d):,} ({100*filled/len(d):.1f}%)")

    d["method"] = METHOD
    d = d[["method", "locus", "rsid", "bp", "a1", "a2", "pip"]]
    d.to_csv(args.out, index=False, compression="gzip")

    # Reproduce the published headline counts as a self-check, not a claim.
    per_locus = d.groupby("locus")["pip"].max()
    print(f"[check] loci with a variant PIP>0.10 : {(per_locus>0.10).sum():3d}  (published: 136)")
    print(f"[check] loci with a variant PIP>0.95 : {(per_locus>0.95).sum():3d}  (published:  11)")
    print(f"[done] wrote {args.out} ({len(d):,} rows)")


if __name__ == "__main__":
    main()
