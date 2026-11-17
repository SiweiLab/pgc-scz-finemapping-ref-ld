#!/usr/bin/env python3
"""Task B22: build EUR-SuSiE-RSS_polyfun337k's .z/.ld/.master for 249 loci (all 250
except chr3_48182866_53553354, excluded per the brief -- 20.3% of its variant pairs
are >3 Mb apart, outside any single 3 Mb Broad window, so no stitching recovers a real
correlation for them).

Windows: task_b19_out/all250_best_window_v2.tsv (B19's corrected, no-liftover picker --
native hg19 locus coordinates matched directly against the Broad bucket's own window
boundaries). Matching/flip logic: identical to build_polyfun_v2.py's match_and_flip
(rsid + allele-set match, sign-flip on swap, drop on genuine allele mismatch) -- reused,
not reimplemented, scaled from 5 pilot loci to all 249. Same EUR_SUMSTATS_DIR source as
the GCTB arm (B23), so the underlying sumstats are identical by construction; only the
LD panel and its own variant coverage differ.

Writes the same .z/.ld/.master trio as build_finemap_input_gctb.py, so the same R
fitting script fits both arms unchanged. n_samples fixed to N_GWAS_PUBLISHED_EUR.
"""
import concurrent.futures
import gzip
import io
import os
import sys

import numpy as np
import pandas as pd

AOU_LD_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MULTI_ANCESTRY_DIR = os.path.dirname(AOU_LD_DIR)
sys.path.insert(0, AOU_LD_DIR)
sys.path.insert(0, MULTI_ANCESTRY_DIR)
from run_params import EUR_SUMSTATS_DIR, N_GWAS_PUBLISHED_EUR  # noqa: E402
from ld_diagnostic import estimate_s_scaled_legacy  # noqa: E402

RAW_DIR = os.path.expanduser("~/work/task_b19_out/ukbb_broad_raw_v2")
OUT_DIR = os.path.expanduser("~/work/task_b22_out/finemap_eur_polyfun337k")
EUR_SS_DIR = os.path.expanduser(EUR_SUMSTATS_DIR)
WINDOW_MAP = os.path.join(AOU_LD_DIR, "task_b19_out", "all250_best_window_v2.tsv")
LOCI_FILE = os.path.join(AOU_LD_DIR, "task_b10_out", "loci250_b10.txt")
EXCLUDED_LOCUS = "chr3_48182866_53553354"

os.makedirs(OUT_DIR, exist_ok=True)
_window_cache = {}


def load_eur_sumstats(locus):
    with gzip.open(f"{EUR_SS_DIR}/{locus}_sumstats.tsv.gz", "rt") as f:
        return pd.read_csv(f, sep="\t")


def load_window_raw(window):
    """Returns (meta, data, row, col, n) -- the raw COO triplets, no matrix built at
    all. Building the full n x n sparse matrix (coo_matrix().tocsr(), then M + M.T)
    before knowing which few hundred to few thousand indices a locus actually needs
    is what caused a near-OOM on this cluster's densest windows: one window is
    48783 x 48783 at ~50% density -- 1.18 BILLION nonzeros -- and scipy's own
    COO->CSR conversion plus the M + M.T transpose-and-add each make roughly a full
    additional copy, ballooning VmSize past 70 GB (confirmed: killed three times
    before this fix, each kill needed to protect the wider session). Returning the
    raw triplets and filtering them down to the needed submatrix BEFORE building
    anything sparse (see submatrix_from_coo) avoids ever materializing the full
    structure -- the unavoidable floor is just loading data/row/col themselves
    (~25-28 GB for the 1.18B-nonzero window), not the multiplied-up copies on top."""
    if window in _window_cache:
        return _window_cache[window]
    meta = pd.read_csv(f"{RAW_DIR}/{window}.gz", sep="\t")
    # np.load() on a path opens the .npz (a zip archive) and seeks around it in small
    # reads to find each array -- on this cluster's network-attached disk that seek
    # pattern measured ~500 KB/s effective throughput (50+ min for a 1.3 GB file,
    # confirmed via /proc/<pid>/fdinfo read offset while a job was "stuck"). Reading
    # the whole file into memory first with one large sequential read (~29 MB/s,
    # matching a plain `dd` benchmark) and parsing from that buffer instead took 64s
    # total for a 1.5 GB file -- a ~50x speedup, confirmed head-to-head before this
    # fix went in.
    #
    # For files > 2 GiB, a plain f.read() hung indefinitely -- confirmed via
    # /proc/<pid>/fdinfo's read offset sitting motionless at exactly 2147479552
    # bytes (0x7ffff000, Linux's MAX_RW_COUNT). A first fix chunked reads at 1 GiB
    # -- and ALL THREE workers hung again, every one motionless at exactly
    # 1,073,741,824 bytes (2^30), a full chunk in. So this is not specifically
    # about MAX_RW_COUNT -- something about this cluster's fs/kernel/glibc combo
    # stalls after any single very-large read() completes, independent of the
    # exact size chosen, as long as it's still "very large". Chunking much
    # smaller (64 MiB, comfortably far from any GiB-scale boundary this cluster
    # has misbehaved at) is the safer fix -- more read() calls, but each one
    # small enough that whatever this issue is stays clear of it.
    size = os.path.getsize(f"{RAW_DIR}/{window}.npz")
    chunk = 64 << 20  # 64 MiB
    parts = []
    with open(f"{RAW_DIR}/{window}.npz", "rb") as f:
        remaining = size
        while remaining > 0:
            piece = f.read(min(chunk, remaining))
            if not piece:
                raise IOError(f"short read on {window}.npz: {remaining} bytes remaining")
            parts.append(piece)
            remaining -= len(piece)
    buf = b"".join(parts)
    d = np.load(io.BytesIO(buf), allow_pickle=True)
    n = int(d["shape"][0])
    result = (meta, d["data"], d["row"], d["col"], n)
    _window_cache[window] = result
    return result


def submatrix_from_coo(data, row, col, n, idx):
    """Build the small symmetric R submatrix for exactly `idx` (window-space row
    indices, any order) straight from the raw lower-triangle-plus-half-diagonal COO
    triplets -- never touches an n x n matrix. `idx`'s order becomes the output
    matrix's row/col order.

    The raw convention (build_polyfun_v2.py, confirmed against the untouched old
    pilot's matrix): strictly lower triangle at true values, diagonal at HALF value
    (0.5). The original M + M.T recipe recovers a true diagonal of 1.0 because the
    diagonal entry gets added to itself once from M and once from M.T. Scattering
    each kept entry to BOTH (r,c) and (c,r) via np.add.at reproduces that exactly,
    including on the diagonal (r2==c2 adds 0.5 twice into the same cell -> 1.0)."""
    keep = np.zeros(n, dtype=bool)
    keep[idx] = True
    mask = keep[row] & keep[col]
    r, c, v = row[mask], col[mask], data[mask]
    remap = np.full(n, -1, dtype=np.int64)
    remap[idx] = np.arange(len(idx))
    r2, c2 = remap[r], remap[c]
    m = len(idx)
    R = np.zeros((m, m), dtype=np.float64)
    np.add.at(R, (r2, c2), v)
    np.add.at(R, (c2, r2), v)
    return R


def match_and_flip(sumstats, meta):
    """rsid + allele-set match against the window; flip sign on swap, drop on genuine
    allele mismatch. Identical logic to build_polyfun_v2.py's match_and_flip."""
    meta_pos = {row.rsid: (i, row.allele1.upper(), row.allele2.upper())
                for i, row in enumerate(meta.itertuples(index=False))}
    idx, z, rsids, freq1 = [], [], [], []
    for row in sumstats.itertuples(index=False):
        rsid = row.MarkerName
        if rsid not in meta_pos:
            continue
        i, ma1, ma2 = meta_pos[rsid]
        a1, a2 = str(row.Allele1).upper(), str(row.Allele2).upper()
        if a1 == ma1 and a2 == ma2:
            sign = 1.0
        elif a1 == ma2 and a2 == ma1:
            sign = -1.0
        else:
            continue
        idx.append(i)
        z.append(row.Zscore * sign)
        rsids.append(rsid)
        f1 = getattr(row, "Freq1", np.nan)
        freq1.append(f1 if sign > 0 else (1 - f1 if pd.notna(f1) else np.nan))
    return np.array(idx), np.array(z), rsids, freq1


def write_finemap_trio(locus, snp_ids, chrom, pos, a1, a2, maf, z, R, n_samples):
    base = f"{OUT_DIR}/{locus}"
    with open(f"{base}.z", "w") as f:
        f.write("rsid chromosome position allele1 allele2 maf beta se\n")
        for i in range(len(snp_ids)):
            m = maf[i] if np.isfinite(maf[i]) else 0.0
            m = min(m, 1 - m)
            f.write(f"{snp_ids[i]} {chrom} {pos[i]} {a1[i]} {a2[i]} {m:.6g} {z[i]:.6g} 1\n")
    np.savetxt(f"{base}.ld", R, fmt="%.6f")
    with open(f"{base}.master", "w") as f:
        f.write("z;ld;snp;config;cred;log;n_samples\n")
        f.write(f"{base}.z;{base}.ld;{base}.snp;{base}.config;{base}.cred;{base}.log;{n_samples}\n")


def process_locus(locus, window):
    """One locus against its already-cached (in this worker process) window."""
    if os.path.exists(f"{OUT_DIR}/{locus}.master"):
        z_df = pd.read_table(f"{OUT_DIR}/{locus}.z", sep=" ")
        R = np.loadtxt(f"{OUT_DIR}/{locus}.ld")
        s, _, _ = estimate_s_scaled_legacy(z_df["beta"].values.astype(float), R)
        n_gwas = len(load_eur_sumstats(locus))
        return dict(locus=locus, ok=True, window=window, n_sumstats=n_gwas,
                    n_matched=len(z_df), pct_matched=round(100 * len(z_df) / n_gwas, 1),
                    s=round(float(s), 4))
    sumstats = load_eur_sumstats(locus)
    n_gwas = len(sumstats)
    meta, data, row, col, n = load_window_raw(window)
    idx, z, rsids, freq1 = match_and_flip(sumstats, meta)
    if len(idx) == 0:
        return dict(locus=locus, ok=False, msg="zero matched variants")
    # filter the raw COO triplets straight to the needed submatrix -- never builds
    # the full n x n matrix (up to 1.18 BILLION nonzeros on this dataset's densest
    # window; see load_window_raw's docstring)
    R = submatrix_from_coo(data, row, col, n, idx)
    n_poly = len(idx)

    meta_sub = meta.iloc[idx].reset_index(drop=True)
    chrom = locus.split("_")[0].replace("chr", "")
    snp_ids = [f"{locus.split('_')[0]}:{p}:{a1}:{a2}"
               for p, a1, a2 in zip(meta_sub["position"], meta_sub["allele1"],
                                     meta_sub["allele2"])]
    write_finemap_trio(locus, snp_ids, chrom, meta_sub["position"].tolist(),
                        meta_sub["allele1"].str.upper().tolist(),
                        meta_sub["allele2"].str.upper().tolist(),
                        np.array(freq1, dtype=float), z, R, N_GWAS_PUBLISHED_EUR)

    s_poly, _, _ = estimate_s_scaled_legacy(z, R)
    return dict(locus=locus, ok=True, window=window, n_sumstats=n_gwas,
                n_matched=n_poly, pct_matched=round(100 * n_poly / n_gwas, 1),
                s=round(float(s_poly), 4))


def process_window(window, loci_for_window):
    """Runs in a worker process: load `window` ONCE, process every locus that maps to
    it (11 windows are shared by >1 locus, per B8 -- this avoids reloading the same
    2 GB+ npz once per locus sharing it)."""
    out = []
    for locus in loci_for_window:
        try:
            out.append(process_locus(locus, window))
        except Exception as e:
            out.append(dict(locus=locus, ok=False, msg=str(e)))
    # ProcessPoolExecutor reuses worker processes across many submitted window-tasks,
    # and _window_cache is module-level -- left uncleared, each worker accumulates
    # every window (each up to several GB as a CSR matrix) it has ever processed for
    # its entire lifetime, not just the current one. Confirmed as the cause of a slow
    # multi-GB memory climb (32 GB used and rising) that nearly repeated the earlier
    # OOM. Clear after this window's loci are done -- nothing later needs it again.
    _window_cache.pop(window, None)
    return out


def main():
    loci = [l.strip() for l in open(LOCI_FILE) if l.strip()]
    win_map = pd.read_csv(WINDOW_MAP, sep="\t").set_index("locus")["window"].to_dict()

    by_window = {}
    for locus in loci:
        if locus == EXCLUDED_LOCUS:
            print(f"[EXCLUDED] {locus}: no window contains it (20.3% of pairs >3 Mb apart)",
                  flush=True)
            continue
        window = win_map.get(locus)
        if window is None:
            print(f"[FAIL] {locus}: no window in {WINDOW_MAP}", flush=True)
            continue
        if not (os.path.exists(f"{RAW_DIR}/{window}.npz") and os.path.exists(f"{RAW_DIR}/{window}.gz")):
            print(f"[NOT YET DOWNLOADED] {locus}: window {window} not on disk, skipping "
                  f"this pass -- rerun once the download finishes", flush=True)
            continue
        by_window.setdefault(window, []).append(locus)

    n_workers = int(os.environ.get("N_CONCURRENT", 8))
    print(f"[info] {len(by_window)} distinct windows covering "
          f"{sum(len(v) for v in by_window.values())} loci, {n_workers} worker processes",
          flush=True)

    overlap_rows, s_rows = [], []
    n_ok = n_fail = 0
    # max_tasks_per_child=1 was tried here to force a fresh process per window-task
    # (glibc's allocator doesn't reliably hand freed multi-GB numpy/scipy arrays back
    # to the OS, so RSS crept up even with the cache cleared). Reverted: with it set,
    # the pool stopped replacing exited workers after a few cycles (zero worker
    # processes alive, no forward progress, for 25+ min) -- a worse failure than the
    # memory growth it was meant to prevent. The _window_cache.pop() fix alone kept
    # memory roughly flat in practice; rely on that plus periodic manual restarts if
    # RSS climbs again, not this executor option.
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = {ex.submit(process_window, w, ls): w for w, ls in by_window.items()}
        for fut in concurrent.futures.as_completed(futures):
            window = futures[fut]
            try:
                results = fut.result()
            except Exception as e:
                print(f"[WINDOW FAIL] {window}: {e}", flush=True)
                continue
            for r in results:
                if r["ok"]:
                    n_ok += 1
                    overlap_rows.append(dict(locus=r["locus"], window=r["window"],
                                              n_sumstats=r["n_sumstats"], n_matched=r["n_matched"],
                                              pct_matched=r["pct_matched"]))
                    s_rows.append(dict(locus=r["locus"], n_variants=r["n_matched"], s=r["s"]))
                    print(f"[OK] {r['locus']}: window={r['window']} "
                          f"matched={r['n_matched']}/{r['n_sumstats']} ({r['pct_matched']}%) "
                          f"s={r['s']}", flush=True)
                else:
                    n_fail += 1
                    print(f"[FAIL] {r['locus']}: {r['msg']}", flush=True)

    print(f"=== DONE: {n_ok} ok, {n_fail} failed, 1 excluded ({EXCLUDED_LOCUS}) ===", flush=True)
    HERE = os.path.dirname(os.path.abspath(__file__))
    pd.DataFrame(overlap_rows).to_csv(f"{HERE}/polyfun337k_overlap.csv", index=False)
    pd.DataFrame(s_rows).to_csv(f"{HERE}/polyfun337k_s_stats.csv", index=False)
    print(f"wrote {HERE}/polyfun337k_{{overlap,s_stats}}.csv", flush=True)


if __name__ == "__main__":
    main()
