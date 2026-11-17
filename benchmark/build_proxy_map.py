#!/usr/bin/env python3
"""export the r2>0.95 proxy map that unblocks LD-aware sensitivity/FDR
scoring .

Side a = the 389 published SuSiEx benchmark variants (results/benchmark/published_pip_vectors.csv.gz,
method == 'SuSiEx-published'); side b = every variant in that locus's EUR reference-LD
panel with r2 > 0.95 to it. Primary r2 is computed in the EUR panel (per_anc_ld/EUR/{locus}_R.npz
+ _snps.txt.gz -- restored from the workspace backup of the per-ancestry EUR LD,
the same panel the per-ancestry methods' *_r2_to_published columns use). An r2_eas column from the EAS
panel is added alongside (bonus, per the brief) -- EAS panel rows are matched to the same
(rsid_a, rsid_b) pair by MarkerName, since the EAS snps file carries no bp column of its own.

Emits each pair once (a,b symmetric; the metric's load_proxy_map() indexes both directions).
"""
import os
import sys

import numpy as np
import pandas as pd

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "benchmark")
EUR_DIR = os.path.expanduser("~/work/per_anc_ld/EUR")
EAS_DIR = os.path.expanduser("~/work/per_anc_ld/EAS")
PUB_PATH = os.path.join(RESULTS, "published_pip_vectors.csv.gz")
OUT_PATH = os.path.join(RESULTS, "proxy_map_r095.csv.gz")
R2_THRESHOLD = 0.95


def load_panel(locus, d):
    r_path = os.path.join(d, f"{locus}_R.npz")
    s_path = os.path.join(d, f"{locus}_snps.txt.gz")
    if not (os.path.exists(r_path) and os.path.exists(s_path)):
        return None
    R = np.load(r_path)["R"]
    snps = pd.read_csv(s_path, sep="\t")
    if len(snps) != R.shape[0]:
        sys.exit(f"[error] {locus}: snps table ({len(snps)}) / R matrix ({R.shape[0]}) size mismatch")
    return R, snps


def eur_bp(snps):
    # SNP_ID is chr:pos:a1:a2 (hg38 tag) in the EUR panel -- see merge_eur_eas_ld.py.
    return snps["SNP_ID"].str.split(":").str[1].astype(int)


def main():
    pub = pd.read_csv(PUB_PATH)
    side_a = pub[pub["method"] == "SuSiEx-published"].reset_index(drop=True)
    print(f"[info] {len(side_a)} published SuSiEx variants across {side_a['locus'].nunique()} loci")

    rows = []
    n_loci_missing_eur = 0
    n_a_missing_eur = 0
    n_loci_missing_eas = 0
    n_a_missing_eas = 0

    for locus, grp in side_a.groupby("locus"):
        eur = load_panel(locus, EUR_DIR)
        if eur is None:
            n_loci_missing_eur += 1
            n_a_missing_eur += len(grp)
            continue
        R_eur, snps_eur = eur
        snps_eur = snps_eur.copy()
        snps_eur["bp"] = eur_bp(snps_eur)
        rsid_to_idx_eur = {}
        for i, rs in enumerate(snps_eur["MarkerName"]):
            rsid_to_idx_eur.setdefault(rs, i)

        eas = load_panel(locus, EAS_DIR)
        if eas is None:
            n_loci_missing_eas += 1
            rsid_to_idx_eas = {}
            R_eas = None
        else:
            R_eas, snps_eas = eas
            rsid_to_idx_eas = {}
            for i, rs in enumerate(snps_eas["MarkerName"]):
                rsid_to_idx_eas.setdefault(rs, i)

        for _, a in grp.iterrows():
            ia = rsid_to_idx_eur.get(a["rsid"])
            if ia is None:
                # fallback: match by bp (published bp is hg38, same as SNP_ID's pos)
                hit = snps_eur.index[snps_eur["bp"] == int(a["bp"])]
                ia = int(hit[0]) if len(hit) else None
            if ia is None:
                n_a_missing_eur += 1
                continue

            r2_col = R_eur[ia, :] ** 2
            hit_idx = np.where(r2_col > R2_THRESHOLD)[0]
            hit_idx = hit_idx[hit_idx != ia]  # drop self
            if len(hit_idx) == 0:
                continue

            ia_eas = rsid_to_idx_eas.get(a["rsid"]) if rsid_to_idx_eas else None
            if ia_eas is None and rsid_to_idx_eas:
                n_a_missing_eas += 1

            for jb in hit_idx:
                b = snps_eur.iloc[jb]
                r2_eas = np.nan
                if ia_eas is not None:
                    jb_eas = rsid_to_idx_eas.get(b["MarkerName"])
                    if jb_eas is not None and R_eas is not None:
                        r2_eas = float(R_eas[ia_eas, jb_eas]) ** 2
                rows.append({
                    "locus": locus,
                    "rsid_a": a["rsid"], "bp_a": int(a["bp"]), "a1_a": a["a1"], "a2_a": a["a2"],
                    "rsid_b": b["MarkerName"], "bp_b": int(b["bp"]), "a1_b": b["A1"], "a2_b": b["A2"],
                    "r2": float(r2_col[jb]),
                    "r2_eas": r2_eas,
                    "ld_panel": "EUR (per_anc_ld/EUR, GCTB/UKB reference)",
                })

    out = pd.DataFrame(rows)
    # symmetric pairs can arise twice if two different published variants both
    # proxy to the same third variant AND to each other -- dedupe on the
    # unordered (locus, rsid_a, rsid_b) pair, keeping the first.
    key = out.apply(lambda r: (r["locus"],) + tuple(sorted([r["rsid_a"], r["rsid_b"]])), axis=1)
    before = len(out)
    out = out.loc[~key.duplicated()].reset_index(drop=True)
    print(f"[info] {before} raw pairs -> {len(out)} after symmetric-dedup")

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    out.to_csv(OUT_PATH, index=False, compression="gzip")
    print(f"[info] wrote {OUT_PATH}: {len(out)} rows, {out['locus'].nunique()} loci")
    print(f"[info] side-a EUR-panel misses: {n_a_missing_eur} variants "
          f"({n_loci_missing_eur} whole loci with no EUR panel file)")
    print(f"[info] side-a EAS-panel misses: {n_a_missing_eas} variants "
          f"({n_loci_missing_eas} whole loci with no EAS panel file) -- r2_eas is NaN for those")
    print(f"[info] r2 (non-null): {out['r2'].notna().sum()}, "
          f"r2_eas (non-null): {out['r2_eas'].notna().sum()}")


if __name__ == "__main__":
    main()
