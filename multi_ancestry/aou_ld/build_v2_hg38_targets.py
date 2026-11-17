#!/usr/bin/env python3
"""Pipeline v2, step 3a: lift every master variant hg19 -> GRCh38 for the AoU EAS extraction.

Replaces lift_eur_z_to_hg38_targets.py (Category A: EUR SNPs only) and the chr13 pilot's
build_union_hg38_targets.py. The input is the whole master table (RUN_ORDER.md step 1), so the
targets are the union of EUR, EAS and core sumstats: Category B for all 250 loci.

Two departures from the v1 scripts, both deliberate:

1. STRAND. v1 took the first mapping and ignored its strand. A variant that maps to the minus
   strand of GRCh38 has its alleles on the other strand there. Its alleles are therefore
   complemented here. Otherwise Hail's exact (locus, alleles) match either misses it, or, for
   an A/T or C/G variant, matches it with the wrong orientation.
2. UNIQUENESS. A variant with more than one mapping, or landing off its own chromosome, is
   dropped. Two master variants landing on the same GRCh38 (position, allele pair) are both
   dropped, since Hail could not tell them apart. Every drop is written out with its reason.

The `rsid` column holds the **vkey**, not an rsID. build_eas_ld_matrix.py carries that column
through Hail as the row ID, so its output rows come back keyed on vkey with no code change.
PGC3's rsID is `pgc_rsid`. allele1/allele2 are the master's ref_a1/ref_a2 (complemented on a
minus-strand mapping). The loader tries both orderings, so REF/ALT order doesn't matter.

Uses the same pyliftover call as v1: convert_coordinate('chr'+chrom, pos-1), then +1.

Usage:
  python3 build_v2_hg38_targets.py --master v2_inputs/master_expected_laptop.tsv.gz \
      --out v2_inputs/hg38_targets
"""
import argparse
import os
from collections import Counter

import pandas as pd
from pyliftover import LiftOver

COMP = str.maketrans("ACGT", "TGCA")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pad", type=int, default=5000,
                    help="bp beyond the lifted min/max for the Hail query interval (v1's 5000)")
    args = ap.parse_args()

    m = pd.read_csv(args.master, sep="\t", dtype={"chr": int, "bp": int})
    lo = LiftOver("hg19", "hg38")
    rows, why = [], Counter()
    for r in m.itertuples(index=False):
        res = lo.convert_coordinate(f"chr{r.chr}", int(r.bp) - 1)
        rec = dict(locus=r.locus, vkey=r.vkey, pgc_rsid=r.rsid, chr=r.chr, bp=r.bp,
                   ref_a1=r.ref_a1, ref_a2=r.ref_a2, palindromic=r.palindromic)
        if not res:
            rec["status"] = "unmapped"
        elif len(res) > 1:
            rec["status"] = "multi_mapped"
        elif res[0][0] != f"chr{r.chr}":
            rec["status"] = "off_chromosome"
        else:
            _, pos0, strand, _ = res[0]
            a1, a2 = r.ref_a1, r.ref_a2
            if strand == "-":
                a1, a2 = a1.translate(COMP), a2.translate(COMP)
            rec.update(status="ok", hg38_pos=pos0 + 1, strand=strand, allele1=a1, allele2=a2)
        why[rec["status"]] += 1
        rows.append(rec)
    d = pd.DataFrame(rows)

    ok = d.status == "ok"
    # two master variants on one GRCh38 site and unordered allele pair: Hail can't separate them
    site = (d.chr.astype(str) + ":" + d.hg38_pos.astype("Int64").astype(str) + ":"
            + [":".join(sorted((a, b))) if isinstance(a, str) else "" for a, b in
               zip(d.allele1, d.allele2)])
    coll = ok & site.duplicated(keep=False)
    d.loc[coll, "status"] = "hg38_collision"
    why["hg38_collision"] += int(coll.sum())
    why["ok"] -= int(coll.sum())

    os.makedirs(args.out, exist_ok=True)
    d.to_csv(os.path.join(args.out, "liftover_all.tsv.gz"), sep="\t", index=False)
    keep = d[d.status == "ok"]
    ivals = []
    for locus, g in keep.groupby("locus", sort=True):
        g = g.sort_values("hg38_pos")
        t = pd.DataFrame({"rsid": g.vkey, "chromosome": g.chr, "position": g.hg38_pos.astype(int),
                          "allele1": g.allele1, "allele2": g.allele2, "pgc_rsid": g.pgc_rsid,
                          "strand": g.strand})
        t.to_csv(os.path.join(args.out, f"{locus}_hg38_targets.tsv"), sep="\t", index=False)
        lo_bp, hi_bp = int(g.hg38_pos.min()) - args.pad, int(g.hg38_pos.max()) + args.pad
        ivals.append(dict(locus=locus, n_master=int((d.locus == locus).sum()), n_targets=len(g),
                          interval=f"chr{g.chr.iloc[0]}:{max(lo_bp, 1)}-{hi_bp}"))
    pd.DataFrame(ivals).to_csv(os.path.join(args.out, "intervals.tsv"), sep="\t", index=False)

    n = len(d)
    print(f"master {n:,}: " + ", ".join(f"{k} {v:,}" for k, v in why.most_common()))
    print(f"minus-strand mappings (alleles complemented): {int((keep.strand == '-').sum()):,}")
    print(f"targets written: {len(keep):,} over {keep.locus.nunique()} loci -> {args.out}")


if __name__ == "__main__":
    main()
