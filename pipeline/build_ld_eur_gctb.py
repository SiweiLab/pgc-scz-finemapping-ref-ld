#!/usr/bin/env python3
"""Pipeline v2 step 2: EUR LD from the GCTB UKB 13M FullLDM, joined to the master by
position.

For every master row: join to snp.info on (Chrom, PhysPos, unordered {A1, A2}); record
sign = +1 if (GCTB A1) == (ref_a1, ref_a2) as exact strings, -1 if swapped. Read each
needed block once with ld/read_ukb_ld.py::read_gctb_ldm_block (v1's reader), cut every
locus's rows out of it, and assemble per-locus R in master order. Pairs in different blocks
are 0 (as GCTB stores them, as v1). R is NOT flipped: it stays in GCTB's allele coding.
The diagonal is set to 1.0, as v1's builder does (GCTB stores 0.99985-0.99995).

Usage:
  python build_ld_eur_gctb.py MASTER_TSV GCTB_DIR OUT_DIR [--workers N]
Writes OUT_DIR/ld_v2/eur_gctb/{locus}/{R.npz,snps.tsv.gz}, OUT_DIR/gctb_presence.tsv.gz,
OUT_DIR/gctb_multiallelic.tsv, OUT_DIR/build_gctb.log.
"""
import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from read_ukb_ld import read_gctb_ldm_block  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("master")
ap.add_argument("gctb_dir")
ap.add_argument("out_dir")
ap.add_argument("--workers", type=int, default=4)
ap.add_argument("--big-n", type=int, default=20000,
                help="blocks with more SNPs than this are read one at a time (v1's reader "
                     "peaks at ~2.5x the dense block: ~22 GB for the 47k-SNP block)")
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)
LOG = open(os.path.join(args.out_dir, "build_gctb.log"), "a")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")
    LOG.flush()


def join_master(master, si):
    """Position + unordered-allele join. Returns master with GCTB columns and match_note."""
    m = master.copy()
    m["aset"] = [frozenset(p) for p in zip(m["ref_a1"], m["ref_a2"])]
    cand = si.merge(m[["chr", "bp"]].drop_duplicates(), left_on=["Chrom", "PhysPos"],
                    right_on=["chr", "bp"])
    cand["aset"] = [frozenset(p) for p in zip(cand["A1"], cand["A2"])]
    hit = m[["locus", "vkey", "chr", "bp", "aset"]].merge(
        cand[["chr", "bp", "aset", "ID", "A1", "A2", "A1Freq", "Block", "row"]], on=["chr", "bp", "aset"])
    dup = hit.duplicated("vkey", keep=False)
    assert not dup.any(), f"{dup.sum()} master rows match >1 snp.info row"
    m = m.merge(hit.drop(columns=["locus", "chr", "bp", "aset"]), on="vkey", how="left")
    m["in_gctb"] = m["ID"].notna()
    same = (m["A1"] == m["ref_a1"]) & (m["A2"] == m["ref_a2"])
    swap = (m["A1"] == m["ref_a2"]) & (m["A2"] == m["ref_a1"])
    m["sign"] = np.select([m["in_gctb"] & same, m["in_gctb"] & swap], [1, -1], 0)
    assert (m.loc[m["in_gctb"], "sign"] != 0).all()  # unordered match => same or swapped
    # Positions GCTB carries with a different allele pair: record, don't store.
    pos_alleles = cand.groupby(["chr", "bp"]).apply(
        lambda d: ";".join(sorted(f"{a}/{b}" for a, b in zip(d["A1"], d["A2"]))), include_groups=False)
    m["gctb_pairs_at_pos"] = [pos_alleles.get((c, b), "") for c, b in zip(m["chr"], m["bp"])]
    note = np.where(m["in_gctb"], "matched", np.where(m["gctb_pairs_at_pos"] != "",
                    "position_in_gctb_other_alleles", "absent"))
    rs_diff = m["in_gctb"] & (m["ID"] != m["rsid"])
    note = np.where(rs_diff, note + ";rsid_differs", note)
    m["match_note"] = note
    m["eur_maf"] = np.minimum(m["A1Freq"], 1 - m["A1Freq"])
    return m.drop(columns=["aset"])


def cut_block(task):
    """Worker: read one block once, return {locus: (rows_in_block_order_positions, sub)}."""
    bid, per_locus, gctb_dir, si_sub = task
    ld, _ = read_gctb_ldm_block(bid, gctb_dir, os.path.join(gctb_dir, "ldm.info"), si_sub)
    out = {}
    for locus, rows in per_locus.items():
        rows = np.asarray(rows)
        out[locus] = ld[np.ix_(rows, rows)].copy()
    del ld
    return bid, out


def main():
    t0 = time.time()
    master = pd.read_csv(args.master, sep="\t", dtype=str, keep_default_na=False)
    master["chr"], master["bp"] = master["chr"].astype(int), master["bp"].astype(int)
    for c in ["palindromic", "in_eur_ss", "in_eas_ss", "in_core_ss"]:
        assert master[c].isin(["True", "False"]).all(), c
        master[c] = master[c] == "True"
    si = pd.read_csv(os.path.join(args.gctb_dir, "snp.info"), sep=r"\s+",
                     usecols=["Chrom", "ID", "PhysPos", "A1", "A2", "A1Freq", "Block"],
                     dtype={"Chrom": int, "PhysPos": int, "ID": str, "A1": str, "A2": str})
    si["row"] = si.groupby("Block").cumcount()  # as build_ukb_ld_matrix.load_snp_info
    si["A1"], si["A2"] = si["A1"].str.upper(), si["A2"].str.upper()
    say(f"snp.info {len(si):,} rows; master {len(master):,} rows")

    m = join_master(master, si)
    say(f"in GCTB: {m['in_gctb'].sum():,}; notes: {m['match_note'].value_counts().to_dict()}")
    pres = m[["locus", "vkey", "in_gctb", "ID", "sign", "Block", "eur_maf", "match_note",
              "gctb_pairs_at_pos"]].rename(columns={"ID": "src_id", "Block": "block"})
    pres["block"] = pres["block"].astype("Int64")
    pres.to_csv(os.path.join(args.out_dir, "gctb_presence.tsv.gz"), sep="\t", index=False,
                compression="gzip")

    # Multi-allelic positions in the master: how each resolved.
    mp = m.groupby(["locus", "bp"])["vkey"].transform("size") > 1
    m.loc[mp, ["locus", "vkey", "rsid", "ref_a1", "ref_a2", "in_eur_ss", "in_eas_ss", "in_core_ss",
               "in_gctb", "ID", "sign", "match_note", "gctb_pairs_at_pos"]].to_csv(
        os.path.join(args.out_dir, "gctb_multiallelic.tsv"), sep="\t", index=False)

    g = m[m["in_gctb"]].copy()
    g["Block"] = g["Block"].astype(int)
    g["row"] = g["row"].astype(int)
    blocks = sorted(g["Block"].unique())
    missing = [b for b in blocks if not os.path.exists(os.path.join(args.gctb_dir, f"block{b}.ldm.bin"))]
    assert not missing, f"blocks not on disk: {missing}"
    si_sub = os.path.join(args.out_dir, "_snp_info_needed_blocks.txt")
    # read_gctb_ldm_block uses snp.info only to count the block's rows; a copy cut to the
    # needed blocks gives the same count and saves a 1.4 GB scan per block.
    pd.read_csv(os.path.join(args.gctb_dir, "snp.info"), sep=r"\s+").pipe(
        lambda d: d[d["Block"].isin(blocks)]).to_csv(si_sub, sep=" ", index=False)
    say(f"{len(blocks)} blocks needed over {g['locus'].nunique()} loci")

    tasks = []
    for b in blocks:
        gb = g[g["Block"] == b]
        tasks.append((b, {L: d["row"].tolist() for L, d in gb.groupby("locus")}, args.gctb_dir, si_sub))
    nsnp = pd.read_csv(os.path.join(args.gctb_dir, "ldm.info"), sep=r"\s+").set_index("Block")["NumSnps"]
    big = [t for t in tasks if nsnp[t[0]] > args.big_n]
    small = [t for t in tasks if nsnp[t[0]] <= args.big_n]
    say(f"{len(big)} blocks > {args.big_n} SNPs read serially; {len(small)} with {args.workers} workers")
    cuts = {}
    i = 0
    for pool_n, ts in ((1, big), (args.workers, small)):
        with ProcessPoolExecutor(pool_n) as ex:
            for bid, out in ex.map(cut_block, ts):
                i += 1
                for L, sub in out.items():
                    cuts[(L, bid)] = sub
                if i % 20 == 0 or i == len(tasks) or pool_n == 1:
                    say(f"  blocks read {i}/{len(tasks)} (block {bid}, {time.time()-t0:.0f}s)")

    root = os.path.join(args.out_dir, "ld_v2", "eur_gctb")
    for L, d in g.groupby("locus", sort=False):
        d = d.reset_index(drop=True)  # master order (m keeps master order)
        n = len(d)
        R = np.zeros((n, n), dtype=np.float32)
        for b, db in d.groupby("Block"):
            pos = db.index.to_numpy()
            R[np.ix_(pos, pos)] = cuts.pop((L, b))
        np.fill_diagonal(R, 1.0)
        os.makedirs(os.path.join(root, L), exist_ok=True)
        np.savez_compressed(os.path.join(root, L, "R.npz"), R=R)
        d.assign(src_id=d["ID"], src_a1=d["A1"], src_a2=d["A2"], block=d["Block"])[
            ["vkey", "rsid", "src_id", "src_a1", "src_a2", "sign", "block", "eur_maf"]].to_csv(
            os.path.join(root, L, "snps.tsv.gz"), sep="\t", index=False, compression="gzip")
    assert not cuts
    say(f"wrote {g['locus'].nunique()} loci, {len(g):,} variants in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
