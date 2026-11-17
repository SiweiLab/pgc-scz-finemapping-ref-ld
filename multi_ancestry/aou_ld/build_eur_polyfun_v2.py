#!/usr/bin/env python3
"""Pipeline v2 step 2b (B40): EUR LD from Broad/PolyFun UKBB-337k, joined to the master by
position.

Windows: task_b19_out/all250_best_window_v2.tsv (one per locus); chr3_48182866_53553354 is
excluded (wider than any 3 Mb window, as in B22), so 249 loci over 234 windows.

Streaming: download one window's .npz/.gz into a scratch dir, cut every locus it serves,
write the cuts, delete the window. The reader is B22's, imported unchanged:
task_b22_out/build_polyfun_arm_b22.py::load_window_raw (64 MiB sequential reads, raw COO
triplets, no n x n build) and ::submatrix_from_coo (lower triangle + half diagonal, filtered
to the locus's indices before anything is built).

Join on (chr, position, unordered {allele1, allele2}); sign = +1 if (allele1, allele2) ==
(ref_a1, ref_a2) as exact strings, -1 if swapped. R is NOT flipped. The .gz carries no allele
frequency, so eur_maf is NA here. `block` holds the window name.

Usage:
  python build_eur_polyfun_v2.py MASTER_TSV OUT_DIR SCRATCH_DIR [--workers N] [--only W,...]
Writes OUT_DIR/ld_v2/eur_polyfun/{locus}/{R.npz,snps.tsv.gz}, OUT_DIR/polyfun_presence/
{window}.tsv.gz (merged at the end into polyfun_presence.tsv.gz), build_polyfun.log.
Resumable: a window whose loci all have R.npz is skipped.
"""
import argparse
import os
import shutil
import sys
import time
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "task_b22_out"))
import build_polyfun_arm_b22 as b22  # noqa: E402

URL = "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/{}"
EXCLUDED_LOCUS = "chr3_48182866_53553354"


def download(url, dst):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(dst + ".part", "wb") as fh:
                size = int(r.headers["Content-Length"])
                while True:
                    c = r.read(64 << 20)
                    if not c:
                        break
                    fh.write(c)
            assert os.path.getsize(dst + ".part") == size, "short download"
            os.replace(dst + ".part", dst)
            return size
        except Exception as e:  # noqa: BLE001
            print(f"  retry {attempt + 1} {url}: {e!r}", flush=True)
            time.sleep(15)
    raise IOError(f"download failed: {url}")


def join_window(mst, meta, window):
    """Master rows of one locus x window meta, on (chr, position, unordered alleles)."""
    meta = meta.assign(i=np.arange(len(meta)), A1=meta["allele1"].str.upper(),
                       A2=meta["allele2"].str.upper())
    meta["aset"] = [frozenset(p) for p in zip(meta["A1"], meta["A2"])]
    mm = mst.assign(aset=[frozenset(p) for p in zip(mst["ref_a1"], mst["ref_a2"])])
    hit = mm[["vkey", "bp", "aset"]].merge(meta[["position", "aset", "i", "rsid", "A1", "A2"]],
                                           left_on=["bp", "aset"], right_on=["position", "aset"])
    ndup = int(hit.duplicated("vkey", keep=False).sum())
    hit = hit.drop_duplicates("vkey", keep=False)  # >1 window row for one vkey: ambiguous, dropped
    out = mst.merge(hit[["vkey", "i", "rsid", "A1", "A2"]].rename(columns={"rsid": "src_id"}),
                    on="vkey", how="left")
    out["in_polyfun"] = out["i"].notna()
    same = (out["A1"] == out["ref_a1"]) & (out["A2"] == out["ref_a2"])
    swap = (out["A1"] == out["ref_a2"]) & (out["A2"] == out["ref_a1"])
    out["sign"] = np.select([out["in_polyfun"] & same, out["in_polyfun"] & swap], [1, -1], 0)
    posset = meta.groupby("position").apply(
        lambda d: ";".join(sorted(f"{a}/{b}" for a, b in zip(d["A1"], d["A2"]))), include_groups=False)
    out["pf_pairs_at_pos"] = out["bp"].map(posset).fillna("")
    note = np.where(out["in_polyfun"], "matched",
                    np.where(out["pf_pairs_at_pos"] != "", "position_in_polyfun_other_alleles", "absent"))
    note = np.where(out["in_polyfun"] & (out["src_id"] != out["rsid"]), note + ";rsid_differs", note)
    out["match_note"] = note
    out["block"] = window
    return out, ndup


def process_window(window, loci, master_path, out_dir, scratch):
    t0 = time.time()
    root = os.path.join(out_dir, "ld_v2", "eur_polyfun")
    wdir = os.path.join(scratch, window)
    os.makedirs(wdir, exist_ok=True)
    nb = download(URL.format(f"{window}.npz"), os.path.join(wdir, f"{window}.npz"))
    download(URL.format(f"{window}.gz"), os.path.join(wdir, f"{window}.gz"))
    t_dl = time.time() - t0
    b22.RAW_DIR = wdir
    meta, data, row, col, n = b22.load_window_raw(window)
    master = pd.read_csv(master_path, sep="\t", dtype={"chr": int, "bp": int}, keep_default_na=False)
    res, pres = [], []
    for L in loci:
        mst = master[master["locus"] == L].reset_index(drop=True)
        j, ndup = join_window(mst, meta, window)
        pres.append(j[["locus", "vkey", "in_polyfun", "src_id", "sign", "block", "match_note",
                       "pf_pairs_at_pos"]])
        k = j[j["in_polyfun"]].reset_index(drop=True)
        R = b22.submatrix_from_coo(data, row, col, n, k["i"].astype(int).to_numpy()).astype(np.float32)
        os.makedirs(os.path.join(root, L), exist_ok=True)
        k.assign(src_a1=k["A1"], src_a2=k["A2"], eur_maf=np.nan)[
            ["vkey", "rsid", "src_id", "src_a1", "src_a2", "sign", "block", "eur_maf"]].to_csv(
            os.path.join(root, L, "snps.tsv.gz"), sep="\t", index=False, compression="gzip")
        np.savez_compressed(os.path.join(root, L, "R.npz"), R=R)
        res.append(dict(locus=L, window=window, n_master=len(mst), n_polyfun=len(k), dup_dropped=ndup,
                        diag_min=float(np.diag(R).min()) if len(k) else np.nan,
                        diag_max=float(np.diag(R).max()) if len(k) else np.nan))
    b22._window_cache.pop(window, None)
    del meta, data, row, col
    shutil.rmtree(wdir)
    pd.concat(pres).to_csv(os.path.join(out_dir, "polyfun_presence", f"{window}.tsv.gz"), sep="\t",
                           index=False, compression="gzip")
    return window, res, nb, t_dl, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("master")
    ap.add_argument("out_dir")
    ap.add_argument("scratch")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out_dir, "polyfun_presence"), exist_ok=True)
    log = open(os.path.join(a.out_dir, "build_polyfun.log"), "a")

    def say(s):
        print(s, flush=True)
        log.write(s + "\n")
        log.flush()

    loci = sorted(pd.read_csv(a.master, sep="\t", usecols=["locus"])["locus"].unique())
    wm = pd.read_csv(os.path.join(HERE, "task_b19_out", "all250_best_window_v2.tsv"), sep="\t")
    wm = wm[wm["locus"].isin(loci) & (wm["locus"] != EXCLUDED_LOCUS)]
    by_w = wm.groupby("window")["locus"].apply(list).to_dict()
    if a.only:
        by_w = {w: by_w[w] for w in a.only.split(",")}
    root = os.path.join(a.out_dir, "ld_v2", "eur_polyfun")
    todo = {w: ls for w, ls in by_w.items()
            if not (all(os.path.exists(os.path.join(root, L, "R.npz")) for L in ls)
                    and os.path.exists(os.path.join(a.out_dir, "polyfun_presence", f"{w}.tsv.gz")))}
    say(f"{len(by_w)} windows over {sum(map(len, by_w.values()))} loci; {len(todo)} to do; "
        f"{a.workers} workers")
    t0 = time.time()
    tot = 0
    with ProcessPoolExecutor(a.workers) as ex:
        futs = {ex.submit(process_window, w, ls, a.master, a.out_dir, a.scratch): w for w, ls in todo.items()}
        for i, f in enumerate(as_completed(futs), 1):
            w = futs[f]
            try:
                w, res, nb, t_dl, t_all = f.result()
            except Exception as e:  # noqa: BLE001
                say(f"[FAIL] {w}: {e!r}")
                continue
            tot += nb
            say(f"[{i}/{len(todo)}] {w} {nb / 1e9:.2f} GB dl {t_dl:.0f}s total {t_all:.0f}s | "
                + "; ".join(f"{r['locus']} {r['n_polyfun']}/{r['n_master']} dup{r['dup_dropped']} "
                            f"diag[{r['diag_min']:.3f},{r['diag_max']:.3f}]" for r in res)
                + f" | cum {tot / 1e9:.0f} GB, {(time.time() - t0) / 60:.0f} min")
    parts = [os.path.join(a.out_dir, "polyfun_presence", f"{w}.tsv.gz") for w in by_w]
    if all(os.path.exists(p) for p in parts):
        pd.concat(pd.read_csv(p, sep="\t", keep_default_na=False) for p in parts).to_csv(
            os.path.join(a.out_dir, "polyfun_presence.tsv.gz"), sep="\t", index=False, compression="gzip")
        say("wrote polyfun_presence.tsv.gz")


if __name__ == "__main__":
    main()
