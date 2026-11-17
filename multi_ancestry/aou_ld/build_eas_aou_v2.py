#!/usr/bin/env python3
"""Pipeline v2 step 3c (B41): align the AoU EAS extraction to the master. Replaces
align_eas_ld_to_sumstats.py.

Input per locus: build_eas_ld_npz.py's {locus}_eas_R.npz + {locus}_eas_snps.tsv (rows of the
extraction's .z, NaN-diagonal rows already dropped) and the v2 target file.
- Key: the .z `rsid` column is the vkey (Hail carried the target's row ID). No position or
  rsID join.
- src_a1/src_a2 = the .z allele1/allele2 = ALT/REF. R is ALT-dosage correlation, so it is
  coded to src_a1 (confirmed against v1's hail_A1 at chr12 by the B41 checks).
- sign: (src_a1, src_a2) vs the target's (allele1, allele2), which are the master's
  ref_a1/ref_a2, complemented on minus-strand rows, i.e. on Hail's strand. Exact same = +1,
  swapped = -1, anything else is not stored and is listed. A vkey on more than one extracted
  row is ambiguous and not stored either. There is no unresolved state.
- R is NOT flipped. eas_af and call_rate are carried from the .z when present (B45), else NA.

Usage: python build_eas_aou_v2.py --loci L1,L2 --extract-dir DIR --out-dir DIR [--anc eas]
  DIR/{locus}/ holds the npz outputs. Writes OUT/ld_v2/{anc}_aou/{locus}/{R.npz,snps.tsv.gz}
  and OUT/{anc}_aou_not_stored.tsv. --anc (B48) names the input files, the output tree and the
  AF column ({anc}_af); the default keeps B45's EAS behaviour.
"""
import argparse
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def align(locus, extract_dir, master, anc="eas"):
    d = os.path.join(extract_dir, locus)
    R = np.load(os.path.join(d, f"{locus}_{anc}_R.npz"))["R"]
    z = pd.read_csv(os.path.join(d, f"{locus}_{anc}_snps.tsv"), sep=" ", dtype=str)
    assert len(z) == R.shape[0]
    tg = pd.read_csv(os.path.join(HERE, "v2_inputs", "hg38_targets", f"{locus}_hg38_targets.tsv"),
                     sep="\t", dtype=str).set_index("rsid")
    mst = master[master["locus"] == locus].set_index("vkey")
    z = z.rename(columns={"rsid": "vkey"})
    z["src_a1"], z["src_a2"] = z["allele1"].str.upper(), z["allele2"].str.upper()
    dup = z["vkey"].duplicated(keep=False)
    unknown = ~z["vkey"].isin(tg.index)
    t1 = z["vkey"].map(tg["allele1"]).str.upper()
    t2 = z["vkey"].map(tg["allele2"]).str.upper()
    same = (z["src_a1"] == t1) & (z["src_a2"] == t2)
    swap = (z["src_a1"] == t2) & (z["src_a2"] == t1)
    z["sign"] = np.select([same, swap], [1, -1], 0)
    reason = np.select([unknown, dup, z["sign"] == 0], ["not_a_target", "duplicate_vkey", "allele_mismatch"], "")
    keep = reason == ""
    bad = z.loc[~keep, ["vkey", "chromosome", "position", "src_a1", "src_a2"]].assign(
        locus=locus, reason=reason[~keep], tgt_a1=t1[~keep], tgt_a2=t2[~keep])
    k = np.where(keep)[0]
    zs = z.iloc[k].reset_index(drop=True)
    Rs = R[np.ix_(k, k)]
    out = pd.DataFrame({
        "vkey": zs["vkey"],
        "rsid": zs["vkey"].map(mst["rsid"]),
        "src_id": "chr" + zs["chromosome"] + ":" + zs["position"] + ":" + zs["src_a2"] + ":" + zs["src_a1"],
        "src_a1": zs["src_a1"], "src_a2": zs["src_a2"], "sign": zs["sign"],
        "strand": zs["vkey"].map(tg["strand"])})
    # B45: the extraction's EAS-subset ALT frequency and call rate, when the .z carries them.
    for c in (f"{anc}_af", "call_rate"):
        out[c] = zs[c].values if c in zs else np.nan
    assert out["rsid"].notna().all() and out["sign"].isin([1, -1]).all()
    return Rs, out, bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci", required=True)
    ap.add_argument("--extract-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--anc", default="eas")
    a = ap.parse_args()
    master = pd.read_csv(os.path.join(HERE, "v2_inputs", "master_variants_250.tsv.gz"), sep="\t",
                         dtype=str, keep_default_na=False)
    bads = []
    for L in a.loci.split(","):
        Rs, out, bad = align(L, os.path.expanduser(a.extract_dir), master, a.anc)
        root = os.path.join(os.path.expanduser(a.out_dir), "ld_v2", f"{a.anc}_aou", L)
        os.makedirs(root, exist_ok=True)
        np.savez_compressed(os.path.join(root, "R.npz"), R=Rs)
        out.to_csv(os.path.join(root, "snps.tsv.gz"), sep="\t", index=False, compression="gzip")
        bads.append(bad)
        print(f"{L}: stored {len(out):,} (sign +1 {int((out['sign'] == 1).sum()):,}, -1 "
              f"{int((out['sign'] == -1).sum()):,}); not stored {len(bad)} {bad['reason'].value_counts().to_dict()}")
    pd.concat(bads).to_csv(os.path.join(os.path.expanduser(a.out_dir), f"{a.anc}_aou_not_stored.tsv"),
                           sep="\t", index=False)


if __name__ == "__main__":
    main()
