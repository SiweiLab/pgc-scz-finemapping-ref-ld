#!/usr/bin/env python3
"""Task B39 (pipeline v2, step 1): the master variant table.

One row per variant present in any of the EUR, EAS or core per-locus sumstats inside the
locus window, keyed on vkey = chr:bp:A:B (GRCh37, alleles upper-cased and sorted). rsid and
ref_a1/ref_a2 come from EUR, else EAS, else core (README.md, Pipeline design). Read-only on its inputs.

Usage: python build_master_b39.py SS_ROOT LOCI_TXT OUT_DIR
  SS_ROOT holds {eur,eas,core}_sumstats/{locus}_sumstats.tsv.gz; OUT_DIR gets
  {locus}_variants.tsv.gz, master_variants_250.tsv.gz and build.log.
"""
import os
import re
import sys
import pandas as pd

ss_root, loci_txt, out_dir = sys.argv[1:4]
os.makedirs(out_dir, exist_ok=True)
log = open(os.path.join(out_dir, "build.log"), "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    log.write(s + "\n")
    log.flush()


ANC = ["eur", "eas", "core"]
COLS = ["locus", "vkey", "rsid", "chr", "bp", "ref_a1", "ref_a2", "palindromic",
        "in_eur_ss", "in_eas_ss", "in_core_ss"]
PAL = {frozenset("AT"), frozenset("CG")}
loci = [l.strip() for l in open(loci_txt) if l.strip()]
assert len(loci) == len(set(loci)) == 250, len(loci)

per_anc = {a: [] for a in ANC}
outside = {a: 0 for a in ANC}
for L in loci:
    c, s, e = L.split("_")
    chrom, start, end = int(c.replace("chr", "")), int(s), int(e)
    for a in ANC:
        d = pd.read_csv(os.path.join(ss_root, f"{a}_sumstats", f"{L}_sumstats.tsv.gz"),
                        sep="\t", dtype=str, keep_default_na=False)
        bad = ~d["POS"].str.fullmatch(r"\d+")
        assert not bad.any(), f"{a} {L}: non-integer POS {d.loc[bad, 'POS'].head().tolist()}"
        d["bp"] = d["POS"].astype("int64")
        assert d["bp"].dtype == "int64"
        inwin = (d["bp"] >= start) & (d["bp"] <= end)
        outside[a] += int((~inwin).sum())
        d = d[inwin].copy()
        d["locus"], d["chr"] = L, chrom
        d["a1"], d["a2"] = d["Allele1"].str.upper(), d["Allele2"].str.upper()
        lo = d[["a1", "a2"]].min(axis=1)
        hi = d[["a1", "a2"]].max(axis=1)
        d["vkey"] = d["chr"].astype(str) + ":" + d["bp"].astype(str) + ":" + lo + ":" + hi
        per_anc[a].append(d[["locus", "vkey", "chr", "bp", "MarkerName", "a1", "a2", "SNP_ID"]])

S = {a: pd.concat(per_anc[a], ignore_index=True) for a in ANC}
for a in ANC:
    d = S[a]
    say(f"{a}: {len(d):,} rows in window ({outside[a]} outside window dropped); "
        f"dup (locus,vkey): {d.duplicated(['locus', 'vkey']).sum()}; "
        f"SNP_ID chr:pos agrees with chr/POS: "
        f"{(d['SNP_ID'].str.split(':').str[1].astype('int64') == d['bp']).sum():,}/{len(d):,}")
    assert not d.duplicated(["locus", "vkey"]).any(), f"duplicate vkey within {a}"

# Union, with EUR > EAS > core precedence for rsid and ref alleles.
keys = pd.concat([S[a][["locus", "vkey", "chr", "bp"]] for a in ANC]).drop_duplicates(["locus", "vkey"])
m = keys.copy()
for a in ANC:
    x = S[a].set_index(["locus", "vkey"])[["MarkerName", "a1", "a2"]].add_prefix(f"{a}_")
    m = m.join(x, on=["locus", "vkey"])
    m[f"in_{a}_ss"] = m[f"{a}_MarkerName"].notna()
m["rsid"] = m["eur_MarkerName"].fillna(m["eas_MarkerName"]).fillna(m["core_MarkerName"])
m["ref_a1"] = m["eur_a1"].fillna(m["eas_a1"]).fillna(m["core_a1"])
m["ref_a2"] = m["eur_a2"].fillna(m["eas_a2"]).fillna(m["core_a2"])
m["palindromic"] = [frozenset((x, y)) in PAL for x, y in zip(m["ref_a1"], m["ref_a2"])]
# Numeric (chr, locus start, bp, vkey): the laptop reference's order.
m["_start"] = m["locus"].str.split("_").str[1].astype("int64")
m = m.sort_values(["chr", "_start", "bp", "vkey"]).drop(columns="_start").reset_index(drop=True)
m.to_pickle(os.path.join(out_dir, "_master_full.pkl"))  # with per-ancestry columns, for checks

out = m[COLS]
out.to_csv(os.path.join(out_dir, "master_variants_250.tsv.gz"), sep="\t", index=False, compression="gzip")
pl = os.path.join(out_dir, "per_locus")
os.makedirs(pl, exist_ok=True)
for L, g in out.groupby("locus"):
    g.to_csv(os.path.join(pl, f"{L}_variants.tsv.gz"), sep="\t", index=False, compression="gzip")
say(f"master: {len(out):,} rows over {out['locus'].nunique()} loci; palindromic {out['palindromic'].sum():,}")
