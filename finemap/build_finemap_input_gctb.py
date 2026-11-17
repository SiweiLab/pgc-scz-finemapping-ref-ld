#!/usr/bin/env python3
"""
Build FINEMAP (Benner) inputs for one locus from the GCTB-UKB-LD pipeline outputs
(build_ukb_ld_matrix.py + prep_pgc_locus.py) — i.e. the SAME z and LD that
run_susie_gctb.R uses, just reformatted, so FINEMAP runs on identical data.

Inputs (per locus):
  <ld_dir>/{locus}_ukb_R.npz        key "R": signed-r, aligned to sumstats Allele1,
                                    row/col order = rows of {locus}_ukb_snps.txt.gz
  <ld_dir>/{locus}_ukb_snps.txt.gz  cols incl SNP_ID, MarkerName  (defines R order)
  <sumstats_dir>/{locus}_sumstats.tsv.gz  cols incl MarkerName, Zscore, Weight, Freq1

Writes (in --out), same format as build_finemap_input.py:
  {locus}.z       header: rsid chromosome position allele1 allele2 maf beta se
                  (rsid = SNP_ID; beta = z, se = 1, so FINEMAP z = beta/se = z)
  {locus}.ld      space-delimited signed-r matrix, same SNP order as .z
  {locus}.master  z;ld;snp;config;cred;log;n_samples  (+ one data row)
  {locus}_snps.txt.gz  only with --drop-ambiguous: the ukb_snps rows actually written, in
                  .z order, for run_susie_gctb_aou.R's snps_file (B34 §2b)

--drop-ambiguous (B34) drops A/T and C/G variants (alleles from SNP_ID, the same ones
written to .z) after the sumstats alignment, from the one index that builds .z, .ld and the
snps table. Without it the output is unchanged.

Run: finemap --sss --in-files {locus}.master --n-causal-snps 5 --log
     (k_max=5 matches the PGC3-SCZ paper; we previously used 10. See
      multi_ancestry/archive/FINEMAP_replication_notes.md.)

Usage:
  python build_finemap_input_gctb.py <locus> --ld-dir <eur_ukb_ld> \\
      --sumstats-dir <eur_sumstats> --out <finemap_in> [--n <override>]
"""
import os
import argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("locus")
    ap.add_argument("--ld-dir", required=True)
    ap.add_argument("--sumstats-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=None,
                    help="n_samples override; default = round(median Weight of matched SNPs)")
    ap.add_argument("--drop-ambiguous", action="store_true",
                    help="drop A/T and C/G variants and write the matching {locus}_snps.txt.gz")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    R = np.load(os.path.join(args.ld_dir, f"{args.locus}_ukb_R.npz"))["R"].astype(float)
    ld_snps = pd.read_csv(os.path.join(args.ld_dir, f"{args.locus}_ukb_snps.txt.gz"), sep="\t")
    ss = pd.read_csv(os.path.join(args.sumstats_dir, f"{args.locus}_sumstats.tsv.gz"), sep="\t")

    assert R.shape[0] == len(ld_snps), \
        f"R {R.shape} vs {len(ld_snps)} snps in ukb_snps file"

    # Align sumstats to the LD SNP order (LD order defines R rows/cols).
    ss = ss.set_index("MarkerName")
    keep_mask = ld_snps["MarkerName"].isin(ss.index).values
    if not keep_mask.all():
        idx = np.where(keep_mask)[0]
        R = R[np.ix_(idx, idx)]
        ld_snps = ld_snps.iloc[idx].reset_index(drop=True)
    if args.drop_ambiguous:
        PAL = {frozenset(("A", "T")), frozenset(("C", "G"))}
        al = ld_snps["SNP_ID"].astype(str).str.split(":")
        keep = np.array([frozenset((x[2].upper(), x[3].upper())) not in PAL for x in al])
        idx = np.where(keep)[0]
        R = R[np.ix_(idx, idx)]
        ld_snps = ld_snps.iloc[idx].reset_index(drop=True)
    sub = ss.loc[ld_snps["MarkerName"].values]

    z = sub["Zscore"].astype(float).values
    weight = sub["Weight"].astype(float).values
    freq1 = pd.to_numeric(sub.get("Freq1"), errors="coerce").values \
        if "Freq1" in sub.columns else np.full(len(z), np.nan)
    snp_ids = ld_snps["SNP_ID"].astype(str).values

    n_samples = args.n if args.n is not None else int(round(np.nanmedian(weight)))
    P = len(snp_ids)

    base = os.path.join(args.out, args.locus)
    zf, ldf = f"{base}.z", f"{base}.ld"
    with open(zf, "w") as f:
        f.write("rsid chromosome position allele1 allele2 maf beta se\n")
        for i, sid in enumerate(snp_ids):
            c, pos, a1, a2 = sid.split(":")
            f1 = freq1[i]
            maf = min(f1, 1 - f1) if np.isfinite(f1) else 0.0
            f.write(f"{sid} {c.replace('chr','')} {pos} {a1} {a2} "
                    f"{maf:.6g} {z[i]:.6g} 1\n")
    np.savetxt(ldf, R, fmt="%.6f")
    if args.drop_ambiguous:
        assert R.shape == (P, P) and len(ld_snps) == P
        ld_snps.to_csv(f"{base}_snps.txt.gz", sep="\t", index=False, compression="gzip")
    with open(f"{base}.master", "w") as f:
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(f"{zf};{ldf};{base}.snp;{base}.config;{base}.cred;{base}.log;{n_samples}\n")
    print(f"{args.locus}: {P} SNPs, n_samples={n_samples} -> {zf}, {ldf}, {base}.master")


if __name__ == "__main__":
    main()
