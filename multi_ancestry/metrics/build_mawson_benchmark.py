#!/usr/bin/env python3
"""Build Mawson et al.'s credible-set tables into our PIP-vector schema.

Closes a real reproducibility gap: every Mawson comparison so far has read two CSVs that
were made by hand and left in /tmp, which do not survive a reboot and which no script
produced. Task B22 could not run its Mawson comparison at all for this reason.

Their supplementary workbooks repeat each variant once per annotated transcript, so a
variant appears many times with an identical PIP. We deduplicate to one row per
(locus, CS, rsid, BP), which is what "dedup" in the old filenames meant.

**The full SuSiE table is not among the workbooks we hold.** `media-2` is the full FINEMAP
table and rebuilds exactly (21,598 variants / 255 loci, matching ESTABLISHED_FACTS sec 17),
but `media-5` is only "Full SuSiE results of **95% prioritised genes**" -- 2,498 variants
over 104 loci after dedup, against the 11,463 / 232 the full table has. So the SuSiE side
is read from a committed dedup CSV (`mawson_susie_dedup.csv`, same schema, verified to
reproduce sec 17's counts) rather than rebuilt, and `--susie` accepts either. **Re-download
the full SuSiE supplementary table when convenient** so this side is rebuildable too.

Locus IDs are rebuilt as `chr{chromosome}_{start}_{stop}` to match our own convention, and
the known chr15 alias is applied (see reconstruct_published_full_posterior.py).

Alleles are NOT in their tables, so the output carries empty a1/a2 and the metric layer's
`variant_key` falls back to the rsID for these rows. That is the documented fallback, but
it means a Mawson join is rsID-keyed rather than position+allele-keyed -- state it wherever
these numbers are quoted.

Usage:
  python3 build_mawson_benchmark.py --susie <Mawson_media-5.xlsx> \
      --finemap ../inputs/Mawson2026_media-2_FINEMAP_full.xlsx --out-dir .
"""
import argparse
import os

import pandas as pd

LOCUS_ALIAS = {"chr15_34646330_34663925": "chr15_34645987_34664593",
               "chr15_34645987_34659593": "chr15_34645987_34664593"}
KEEP = ["locus", "CS", "chromosome", "start", "stop", "k", "rsid", "BP", "PIP"]


def load(path, label):
    path = os.path.expanduser(path)
    if path.endswith(".csv"):
        d = pd.read_csv(path)
        missing = set(KEEP) - set(d.columns)
        if missing:
            raise SystemExit(f"[error] {path} missing {sorted(missing)}")
        d = d[KEEP].drop_duplicates(subset=["locus", "CS", "rsid", "BP"])
        d["locus_id"] = ("chr" + d.chromosome.astype(str) + "_"
                         + d.start.astype(int).astype(str) + "_"
                         + d.stop.astype(int).astype(str)).replace(LOCUS_ALIAS)
        print(f"  (pre-built CSV) -> {label}: {len(d):,} rows, {d.locus_id.nunique()} loci, "
              f"{d.groupby(['locus_id','BP']).ngroups:,} unique variants")
        return d
    x = pd.ExcelFile(path)
    sheets = [s for s in x.sheet_names if s != "README"]
    frames = []
    for sh in sheets:
        d = pd.read_excel(x, sheet_name=sh, header=0)
        missing = set(KEEP) - set(d.columns)
        if missing:
            raise SystemExit(f"[error] {path}:{sh} missing {sorted(missing)}; "
                             f"found {list(d.columns)[:12]}")
        frames.append(d[KEEP])
        print(f"  {sh}: {len(d):,} rows (annotation-expanded)")
    d = pd.concat(frames, ignore_index=True)
    # one row per variant per credible set: the workbooks repeat a variant per transcript
    d = d.drop_duplicates(subset=["locus", "CS", "rsid", "BP"])
    d["locus_id"] = ("chr" + d.chromosome.astype(str) + "_"
                     + d.start.astype(int).astype(str) + "_"
                     + d.stop.astype(int).astype(str)).replace(LOCUS_ALIAS)
    print(f"  -> {label}: {len(d):,} unique variant-CS rows, "
          f"{d.locus_id.nunique()} loci, {d.groupby(['locus_id','BP']).ngroups:,} unique variants")
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--susie", default="mawson_susie_dedup.csv",
                    help="the full SuSiE table: a dedup CSV, or an .xlsx if one is obtained")
    ap.add_argument("--finemap", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "inputs", "Mawson2026_media-2_FINEMAP_full.xlsx"))
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    for path, label, method in ((args.susie, "SuSiE", "Mawson-SuSiE"),
                                (args.finemap, "FINEMAP", "Mawson-FINEMAP")):
        print(f"{label} <- {path}")
        d = load(path, label)
        d.to_csv(os.path.join(args.out_dir, f"mawson_{label.lower()}_dedup.csv"), index=False)
        # PIP-vector schema: one row per variant, its best PIP across that locus's sets
        v = (d.groupby(["locus_id", "BP"])
               .agg(pip=("PIP", "max"), rsid=("rsid", "first")).reset_index()
               .rename(columns={"locus_id": "locus", "BP": "bp"}))
        v["method"], v["a1"], v["a2"] = method, pd.NA, pd.NA
        out = os.path.join(args.out_dir, f"pip_vectors_{method.lower()}.csv.gz")
        v[["method", "locus", "rsid", "bp", "a1", "a2", "pip"]].to_csv(out, index=False)
        print(f"  wrote {out}\n")


if __name__ == "__main__":
    main()
