#!/usr/bin/env python3
"""Derive the manuscript's 175/122/53 locus strata from Yuan et al. Suppl. Table 26.

The strata are defined in hierarchy_v3_gws_confident.py against files that live only
inside AoU. This reconstructs them from the published workbook alone, so the
variant-level metric can be scored on exactly the manuscript's own locus set.

Validated: reproduces 175 single-CS loci, 122 confident (max PIP>=0.1) and 53 weak
exactly, so the GWS filter the manuscript also applies does not change these counts
for the published SuSiEx arm -- Yuan et al. had already applied their own.

Usage:
  python3 derive_confident_stratum.py --xlsx <SuSiEx supp .xlsx> --out-dir .
"""
import argparse

import openpyxl
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True, help="../inputs/Yuan2024_SupplementaryTables_MOESM4.xlsx (41588_2024_1870_MOESM4_ESM.xlsx)")
    ap.add_argument("--sheet", default="Suppl Table 26")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--our-loci", required=True,
                    help="a table carrying our own 250 locus ids in a 'locus' column "
                         "(e.g. task_a_out/pip_vectors_full250_qcexcl.csv.gz)")
    args = ap.parse_args()

    w = openpyxl.load_workbook(args.xlsx, read_only=True, data_only=True)
    rows = list(w[args.sheet].iter_rows(min_row=2, values_only=True))
    hdr = [str(h).replace("\n", " ").strip() if h else "" for h in rows[0]]
    data = [r for r in rows[1:] if r and r[0]]
    d = pd.DataFrame(data, columns=hdr[:len(data[0])])
    d["max_pip"] = pd.to_numeric(d["Max PIP"], errors="coerce")
    d["locus"] = d["Loci"].astype(str).str.replace(":", "_").str.replace("-", "_")

    # Locus windows are PGC3-SCZ's own (loci_250_pgc3scz.tsv), which is what
    # PGC-FINEMAP used and what our pipeline adopted -- PGC-FINEMAP's published
    # variants land on all 250 of our ids with nothing left over. Yuan et al.'s
    # SuSiEx paper re-derived one window slightly narrower (ST26 calls chr15
    # 34646330-34663925 where PGC has 34645987-34664593, trimmed 343 bp at the
    # start and 668 at the end; reason unknown, and it is the only one of their 193
    # loci that differs). Task A2's export already harmonised the benchmark's
    # variants onto the PGC ids, so a stratum list built from raw ST26 names
    # silently loses that locus when used as a --loci filter. Map by overlap.
    ours = pd.read_csv(args.our_loci, usecols=["locus"])["locus"].drop_duplicates()
    ref = pd.DataFrame({"our": ours})
    parts = ref["our"].str.extract(r"^chr([0-9XY]+)_(\d+)_(\d+)$")
    ref["chrom"], ref["start"], ref["stop"] = parts[0], parts[1].astype(int), parts[2].astype(int)

    theirs = d["locus"].drop_duplicates()
    tp = theirs.str.extract(r"^chr([0-9XY]+)_(\d+)_(\d+)$")
    tp.columns = ["chrom", "start", "stop"]
    tp["start"] = tp["start"].astype(int)
    tp["stop"] = tp["stop"].astype(int)
    tp["theirs"] = theirs.values

    remap, unmatched = {}, []
    for r in tp.itertuples():
        c = ref[(ref["chrom"] == r.chrom) & (ref["start"] < r.stop) & (ref["stop"] > r.start)]
        if len(c) == 1:
            remap[r.theirs] = c["our"].iloc[0]
        elif len(c) == 0:
            unmatched.append(r.theirs)
        else:  # overlaps several of ours -- take the one covering the most of it
            ov = (c["stop"].clip(upper=r.stop) - c["start"].clip(lower=r.start))
            remap[r.theirs] = c.loc[ov.idxmax(), "our"]
    renamed = sum(1 for k, v in remap.items() if k != v)
    print(f"locus ids: {len(remap)} of {len(tp)} mapped onto ours "
          f"({renamed} renamed, {len(unmatched)} unmatched)")
    for u in unmatched:
        print(f"  unmatched (kept as-is): {u}")
    d["locus"] = d["locus"].map(lambda x: remap.get(x, x))

    n_cs = d.groupby("locus").size()
    single = set(n_cs[n_cs == 1].index)
    mx = d[d["locus"].isin(single)].groupby("locus")["max_pip"].max()
    confident = set(mx[mx >= 0.1].index)
    weak = single - confident

    for name, s, expect in (("single175", single, 175),
                            ("confident122", confident, 122),
                            ("weak53", weak, 53)):
        path = f"{args.out_dir}/loci_{name}.txt"
        with open(path, "w") as f:
            f.write("\n".join(sorted(s)) + "\n")
        flag = "OK" if len(s) == expect else f"MISMATCH, expected {expect}"
        print(f"{path}: {len(s)} loci [{flag}]")


if __name__ == "__main__":
    main()
