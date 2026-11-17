#!/usr/bin/env python3
"""Pipeline v2, step 6: one harmoniser for every arm (PIPELINE.md §4 and §6).

Replaces harmonize_loci.py, harmonize_loci_v2_union.py, merge_eur_eas_ld.py, the alignment in
build_finemap_input_gctb.py, and subset_noambig_b34.py. It reads only through v2_io, so every
z and R arrives indexed by vkey and signed to the master's ref_a1. **This file flips no
allele.** Each arm is a set expression over the loaded sets, and the merged LD is computed here
at run time (step 4).

Outputs go into a run's scratch tree, never the bucket, and use v1's file names, so
run_primary_2026.py only needs its paths pointed here:

  OUT/union/{l}_susiex/{EUR,EAS}.{sst.txt,ld.bin} {EUR,EAS}_ref.bim {EUR,EAS}_frq.frq   SuSiEx
  OUT/union/{l}_harmonized.npz            MultiSuSiE: both-LD domain, missing z = NaN, R real
  OUT/union_censored/{l}_harmonized.npz   MultiSuSiE-censored: FULL union, and per ancestry a
                                          variant without z or LD has z = 0, R row/col = 0
  OUT/union_censored_bothld/...npz        sensitivity: censored on the both-LD domain, so
                                          (vs union/) censoring alone, (vs union_censored/) the
                                          extra variants alone (option C, 2026-09-30)
  OUT/intersect/{l}.{EUR,EAS}.{ma,ld.gz} {l}.snps.tsv                                  MESuSiE
  OUT/merged/{l}.z {l}.ld                 SuSiE-merged, FINEMAP-merged, RSparsePro, CARMA
  OUT/eur/{l}.z {l}.ld {l}.master {l}_snps.txt.gz                   EUR-SuSiE-RSS, EUR-FINEMAP
  OUT/_logs/{l}.json                      domains, merge weights, dropped counts

Domains ("valid" = has sumstats and LD in that ancestry; --drop-ambiguous removes palindromic
vkeys from everything first):
  SuSiEx        each ancestry's own valid set (ragged)
  MultiSuSiE    (EUR valid | EAS valid) & EUR LD & EAS LD      (and censored_bothld)
  MultiSuSiE-censored   EUR valid | EAS valid
  MESuSiE       EUR valid & EAS valid
  merged arms   EUR LD & EAS LD & core sumstats
  EUR arms      EUR valid
Merge weights: the median sumstats Weight (NEFF) of each ancestry over EUR LD & EAS LD, taken
before the core and palindromic filters, as merge_eur_eas_ld.py did.

The tool-level variant ID is PGC3's rsid from the master (unique within every locus). A1/A2 are
always the master's ref_a1/ref_a2, in every file and every ancestry.

Usage:
  python harmonize_v2.py --loci chr12_110129434_111489434 --out ~/work/v2_inputs_run [--drop-ambiguous]
      [--arms susiex,multisusie,mesusie,merged,eur] [--eur-ld eur_gctb] [--ld-root D] [--ss-root D]
"""
import argparse
import gzip
import json
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_params as P  # noqa: E402
import v2_io  # noqa: E402

ARMS = ("susiex", "multisusie", "mesusie", "merged", "eur")
ANC = {"EUR": "eur", "EAS": "eas"}


def ordered(master, keys):
    """keys in master order."""
    keys = set(keys)
    return [v for v in master.index if v in keys]


def ss_block(ss, vkeys):
    """beta, se, z, freq, n arrays for vkeys (all present in ss), already signed by v2_io.
    Without an Effect column, beta = z * se, as v1's harmonisers did."""
    s = ss.loc[vkeys]
    se = s["se"].to_numpy(float) if "se" in s else np.ones(len(s))
    z = s["z"].to_numpy(float)
    beta = s["beta"].to_numpy(float) if "beta" in s else z * se
    freq = s["freq"].to_numpy(float) if "freq" in s else np.full(len(s), np.nan)
    n = s["n"].to_numpy(float) if "n" in s else np.full(len(s), np.nan)
    return beta, se, z, freq, n


def write_zld(path_base, master, vkeys, z, freq, R):
    """FINEMAP-format .z (rsid chromosome position allele1 allele2 maf beta se, beta = z,
    se = 1) and a %.6f .ld, as build_finemap_input_gctb.py wrote them."""
    m = master.loc[vkeys]
    with open(f"{path_base}.z", "w") as f:
        f.write("rsid chromosome position allele1 allele2 maf beta se\n")
        for rs, c, bp, a1, a2, zi, fi in zip(m["rsid"], m["chr"], m["bp"], m["ref_a1"], m["ref_a2"], z, freq):
            maf = min(fi, 1 - fi) if np.isfinite(fi) else 0.0
            f.write(f"{rs} {c} {bp} {a1} {a2} {maf:.6g} {zi:.6g} 1\n")
    np.savetxt(f"{path_base}.ld", R, fmt="%.6f")


def write_susiex(prefix, master, vkeys, ss, R):
    m = master.loc[vkeys]
    beta, se, _, freq, n = ss_block(ss, vkeys)
    stat = beta / se
    pval = [math.erfc(abs(x) / math.sqrt(2)) for x in stat]
    chrom, bp = m["chr"].to_numpy(), m["bp"].to_numpy()
    snp, a1, a2 = m["rsid"].to_numpy(), m["ref_a1"].to_numpy(), m["ref_a2"].to_numpy()
    pd.DataFrame({"chr": chrom, "snp": snp, "cm": 0, "bp": bp, "a1": a1, "a2": a2}).to_csv(
        f"{prefix}_ref.bim", sep="\t", header=False, index=False)
    pd.DataFrame({"CHR": chrom, "SNP": snp, "A1": a1, "A2": a2, "MAF": freq,
                  "NCHROBS": int(round(2 * np.nanmedian(n)))}).to_csv(f"{prefix}_frq.frq", sep=" ", index=False)
    R.astype(np.float32).tofile(f"{prefix}.ld.bin")
    pd.DataFrame({"chr": chrom, "snp": snp, "bp": bp, "A1": a1, "A2": a2, "beta": beta, "se": se,
                  "stat": stat, "p": pval}).to_csv(f"{prefix}.sst.txt", sep="\t", index=False)


def harmonize(locus, out, arms, drop_ambiguous, eur_ld, ld_root, ss_root):
    master = v2_io.load_master(locus)
    allowed = set(master.index[~master["palindromic"]]) if drop_ambiguous else set(master.index)
    ss = {a: v2_io.load_ss(ANC[a], locus, ss_root=ss_root) for a in ANC}
    ld = {"EUR": v2_io.load_ld(eur_ld, locus, ld_root=ld_root),
          "EAS": v2_io.load_ld("eas_aou", locus, ld_root=ld_root)}
    ld_keys = {a: set(ld[a][1].index) for a in ANC}
    valid = {a: set(ss[a].index) & ld_keys[a] & allowed for a in ANC}
    log = dict(locus=locus, drop_ambiguous=drop_ambiguous, eur_ld=eur_ld,
               n_ss={a: len(ss[a]) for a in ANC}, n_ld={a: len(ld_keys[a]) for a in ANC},
               n_valid={a: len(valid[a]) for a in ANC}, domains={})

    def ld_sub(a, vkeys):
        return v2_io.subset(*ld[a], vkeys)[0]

    if "susiex" in arms:
        d = os.path.join(out, "union", f"{locus}_susiex")
        os.makedirs(d, exist_ok=True)
        for a in ANC:
            vk = ordered(master, valid[a])
            write_susiex(os.path.join(d, a), master, vk, ss[a], ld_sub(a, vk))
            log["domains"][f"susiex_{a}"] = len(vk)

    if "multisusie" in arms:
        n_med = {a: (np.nanmedian(ss_block(ss[a], own)[4]) if (own := ordered(master, valid[a])) else np.nan)
                 for a in ANC}

        def write_npz(tree, vk, censor):
            """Uncensored: missing z = NaN, R real everywhere (needs LD in both ancestries).
            Censored (MultiSuSiE's z-path contract, susiepy_ss.py:98-106): for each ancestry, a
            variant that isn't valid there (no z, or no LD) gets z = 0 and an all-zero R row and
            column, so its LD in that ancestry is never needed."""
            K, Pn = len(ANC), len(vk)
            pos = {v: i for i, v in enumerate(vk)}
            arr = {key: np.full((K, Pn), np.nan) for key in ("beta", "se", "z", "freq")}
            R = np.zeros((K, Pn, Pn), dtype=np.float32)
            miss = {}
            for k, a in enumerate(ANC):
                present = [v for v in vk if v in valid[a]]
                idx = np.array([pos[v] for v in present], dtype=int)
                if present:
                    b, s_, z, f, _ = ss_block(ss[a], present)
                    for key, val in zip(("beta", "se", "z", "freq"), (b, s_, z, f)):
                        arr[key][k, idx] = val
                if censor:
                    if present:
                        R[k][np.ix_(idx, idx)] = ld_sub(a, present)
                    for key in ("z", "beta"):
                        arr[key][k] = np.where(np.isfinite(arr[key][k]), arr[key][k], 0.0)
                else:
                    R[k] = ld_sub(a, vk)
                miss[a] = Pn - len(present)
            m = master.loc[vk]
            np.savez_compressed(
                os.path.join(out, tree, f"{locus}_harmonized.npz"), R=R, **arr,
                ancestries=np.array(list(ANC)),
                snp_ids=np.array([f"{c}:{p}:{x}:{y}" for c, p, x, y in
                                  zip(m["chr"], m["bp"], m["ref_a1"], m["ref_a2"])]),
                rsids=m["rsid"].to_numpy(), chr=m["chr"].to_numpy(), pos=m["bp"].to_numpy(),
                effect_allele=m["ref_a1"].to_numpy(), other_allele=m["ref_a2"].to_numpy(),
                n=np.array([n_med[a] for a in ANC]), varY=np.ones(K))
            log["domains"][tree] = Pn
            log["domains"][f"{tree}_missing"] = miss

        both_ld = ordered(master, (valid["EUR"] | valid["EAS"]) & ld_keys["EUR"] & ld_keys["EAS"] & allowed)
        full_union = ordered(master, valid["EUR"] | valid["EAS"])
        for tree, vk, censor in (("union", both_ld, False), ("union_censored_bothld", both_ld, True),
                                 ("union_censored", full_union, True)):
            os.makedirs(os.path.join(out, tree), exist_ok=True)
            write_npz(tree, vk, censor)
        log["domains"]["multisusie"] = log["domains"]["union"]

    if "mesusie" in arms:
        vk = ordered(master, valid["EUR"] & valid["EAS"])
        d = os.path.join(out, "intersect")
        os.makedirs(d, exist_ok=True)
        m = master.loc[vk]
        pd.DataFrame({"rsid": m["rsid"], "chr": m["chr"], "pos": m["bp"], "effect_allele": m["ref_a1"],
                      "other_allele": m["ref_a2"]}).to_csv(os.path.join(d, f"{locus}.snps.tsv"),
                                                           sep="\t", index=False)
        for a in ANC:
            b, s_, _, f, nn = ss_block(ss[a], vk)
            pd.DataFrame({"SNP": m["rsid"].to_numpy(), "A1": m["ref_a1"].to_numpy(), "A2": m["ref_a2"].to_numpy(),
                          "freq": f, "b": b, "se": s_, "N": int(np.nanmedian(nn))}).to_csv(
                os.path.join(d, f"{locus}.{a}.ma"), sep="\t", index=False)
            with gzip.open(os.path.join(d, f"{locus}.{a}.ld.gz"), "wt") as fh:
                np.savetxt(fh, ld_sub(a, vk), fmt="%.8g", delimiter="\t")
        log["domains"]["mesusie"] = len(vk)

    if "merged" in arms:
        both = ordered(master, ld_keys["EUR"] & ld_keys["EAS"])
        w = {}
        for a in ANC:
            have = [v for v in both if v in ss[a].index]
            w[a] = float(np.nanmedian(ss_block(ss[a], have)[4])) if have else float("nan")
        ss_core = v2_io.load_ss("core", locus, ss_root=ss_root)
        vk = [v for v in both if v in ss_core.index and v in allowed]
        Rm = (w["EUR"] * ld_sub("EUR", vk).astype(np.float64)
              + w["EAS"] * ld_sub("EAS", vk).astype(np.float64)) / (w["EUR"] + w["EAS"])
        np.fill_diagonal(Rm, 1.0)
        d = os.path.join(out, "merged")
        os.makedirs(d, exist_ok=True)
        _, _, z, f, _ = ss_block(ss_core, vk)
        write_zld(os.path.join(d, locus), master, vk, z, f, Rm)
        log["domains"]["merged"] = len(vk)
        log["merge_weights"] = w
        log["merged_ld_both"] = len(both)
        log["n_ss"]["core"] = len(ss_core)

    if "eur" in arms:
        vk = ordered(master, valid["EUR"])
        d = os.path.join(out, "eur")
        os.makedirs(d, exist_ok=True)
        _, _, z, f, _ = ss_block(ss["EUR"], vk)
        base = os.path.join(d, locus)
        write_zld(base, master, vk, z, f, ld_sub("EUR", vk))
        rs = master.loc[vk, "rsid"].to_numpy()
        pd.DataFrame({"SNP_ID": rs, "MarkerName": rs, "vkey": vk}).to_csv(
            f"{base}_snps.txt.gz", sep="\t", index=False, compression="gzip")
        with open(f"{base}.master", "w") as fh:
            fh.write("z;ld;snp;config;cred;log;n_samples\n")
            fh.write(f"{base}.z;{base}.ld;{base}.snp;{base}.config;{base}.cred;{base}.log;"
                     f"{P.N_GWAS_PUBLISHED_EUR}\n")
        log["domains"]["eur"] = len(vk)

    log["dropped_ss"] = {a: ss[a].attrs.get("dropped") for a in ANC}
    os.makedirs(os.path.join(out, "_logs"), exist_ok=True)
    lp = os.path.join(out, "_logs", f"{locus}.json")
    if os.path.exists(lp):  # --arms subset: keep the other arms' domains and the merge weights (B54)
        with open(lp) as fh:
            old = json.load(fh)
        old_dom = old.get("domains", {})
        old.update({k: v for k, v in log.items() if k != "domains"})
        old["domains"] = {**old_dom, **log["domains"]}
        log = old
    with open(lp, "w") as fh:
        json.dump(log, fh, indent=1, default=str)
    return log


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", help="comma-separated")
    ap.add_argument("--loci-file")
    ap.add_argument("--out", required=True)
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--drop-ambiguous", action="store_true")
    ap.add_argument("--eur-ld", default="eur_gctb", choices=["eur_gctb", "eur_aou", "eur_polyfun"])
    ap.add_argument("--ld-root")
    ap.add_argument("--ss-root")
    a = ap.parse_args()
    loci = a.loci.split(",") if a.loci else [l.strip() for l in open(a.loci_file) if l.strip()]
    arms = a.arms.split(",")
    bad = set(arms) - set(ARMS)
    if bad:
        sys.exit(f"unknown arms {sorted(bad)}")
    for L in loci:
        log = harmonize(L, os.path.expanduser(a.out), arms, a.drop_ambiguous, a.eur_ld, a.ld_root, a.ss_root)
        print(L, json.dumps(log["domains"]), json.dumps(log.get("merge_weights", {})), flush=True)


if __name__ == "__main__":
    main()
