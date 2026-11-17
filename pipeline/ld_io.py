#!/usr/bin/env python3
"""Pipeline v2, step 5: the one loader, and the only code in the pipeline that flips alleles.

Everything downstream (harmonize.py, the merged LD, every arm's inputs) reads LD and sumstats
through this module and gets them indexed by vkey and signed to the master's ref_a1. Nothing
else may re-align an allele (README.md, Pipeline design; gate 6 greps for it).

Stored layout (README.md, Pipeline design; as promoted, RUN_ORDER.md step 3):
  LD_ROOT/ld_per_anc_{source}/{locus}/R.npz        key "R", n x n, in the SOURCE's own coding
  LD_ROOT/ld_per_anc_{source}/{locus}/snps.tsv.gz  vkey, rsid, src_a1, src_a2, sign, ...
  SS_ROOT/{anc}_sumstats/{locus}_sumstats.tsv.gz   raw PGC3 per-locus sumstats (never signed on disk)
  MASTER                                           ../inputs/master_variants_250.tsv.gz

sources: eur_gctb, eas_aou, eur_aou, eur_polyfun.  ancestries: eur, eas, core.

sign = +1 if (src_a1, src_a2) == (ref_a1, ref_a2) as exact strings, -1 if swapped. It is stored at
build time; load_ld recomputes it from the alleles and refuses a file where they disagree.

Roots come from arguments or the environment (LD_ROOT, SS_ROOT, MASTER_TABLE), and are local
paths: sync the bucket to local disk first.
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES = ("eur_gctb", "eas_aou", "eur_aou", "eur_polyfun")
ANCESTRIES = ("eur", "eas", "core")
DEFAULT_MASTER = os.path.join(HERE, "..", "inputs", "master_variants_250.tsv.gz")

# Optional sumstats columns carried through. z and beta flip with the allele, freq becomes 1 - freq.
SS_OPTIONAL = {"Effect": "beta", "StdErr": "se", "Freq1": "freq", "Weight": "n", "NEFF": "n"}


class V2IOError(ValueError):
    pass


def _root(arg, env):
    r = arg or os.environ.get(env)
    if not r:
        raise V2IOError(f"no root given: pass it or set ${env}")
    return r


def make_vkey(chrom, bp, a, b):
    """chr:bp:A:B, GRCh37, alleles upper-cased and sorted. Vectorised over pandas Series."""
    ab = pd.DataFrame({"a": a.str.upper(), "b": b.str.upper()})
    lo, hi = ab.min(axis=1), ab.max(axis=1)
    return chrom.astype(str).str.replace("chr", "", regex=False) + ":" + bp.astype(str) + ":" + lo + ":" + hi


_COMP = str.maketrans("ACGT", "TGCA")


def complement(alleles):
    return alleles.str.upper().str.translate(_COMP)


def allele_sign(a1, a2, ref_a1, ref_a2):
    """+1 same orientation, -1 swapped, 0 neither. Exact upper-case string match. The only
    orientation rule in the pipeline; strand is handled by the caller complementing a1/a2."""
    a1, a2, r1, r2 = (s.str.upper().to_numpy() for s in (a1, a2, ref_a1, ref_a2))
    return np.select([(a1 == r1) & (a2 == r2), (a1 == r2) & (a2 == r1)], [1, -1], 0)


_master_cache = {}


def load_master(locus, master=None):
    """The master rows for one locus, indexed by vkey, in master order (chr, bp, vkey)."""
    path = master or os.environ.get("MASTER_TABLE") or DEFAULT_MASTER
    if path not in _master_cache:
        m = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        for c in ("palindromic", "in_eur_ss", "in_eas_ss", "in_core_ss"):
            m[c] = m[c] == "True"
        m["bp"] = m["bp"].astype("int64")
        _master_cache[path] = m
    m = _master_cache[path]
    m = m[m["locus"] == locus]
    if m.empty:
        raise V2IOError(f"{locus}: not in the master table")
    if m["vkey"].duplicated().any():
        raise V2IOError(f"{locus}: duplicate vkey in the master")
    return m.set_index("vkey")


def load_ld(source, locus, ld_root=None, master=None, dtype=np.float32):
    """R signed to ref_a1 (R * s s^T) and its rows, indexed by vkey in R's order.

    Returns (R, snps). snps carries rsid (PGC3's, from the master), sign and palindromic.
    A variant the source doesn't carry is absent, never a zero row."""
    if source not in SOURCES:
        raise V2IOError(f"unknown LD source {source!r}; expected one of {SOURCES}")
    d = os.path.join(_root(ld_root, "LD_ROOT"), f"ld_per_anc_{source}", locus)
    R = np.load(os.path.join(d, "R.npz"))["R"]
    s = pd.read_csv(os.path.join(d, "snps.tsv.gz"), sep="\t", dtype=str, keep_default_na=False)
    if R.ndim != 2 or R.shape[0] != R.shape[1] or R.shape[0] != len(s):
        raise V2IOError(f"{source}/{locus}: R {R.shape} vs {len(s)} snps rows")
    if s["vkey"].duplicated().any():
        raise V2IOError(f"{source}/{locus}: duplicate vkey in snps.tsv.gz")
    m = load_master(locus, master)
    missing = ~s["vkey"].isin(m.index)
    if missing.any():
        raise V2IOError(f"{source}/{locus}: {int(missing.sum())} vkeys not in the master, "
                        f"e.g. {s.loc[missing, 'vkey'].head(3).tolist()}")
    mm = m.loc[s["vkey"]]
    sign = s["sign"].astype(int).to_numpy()
    a1, a2 = s["src_a1"], s["src_a2"]
    if "strand" in s.columns:  # AoU (step 3): src alleles on +, signed through the complement on -
        minus = (s["strand"] == "-").to_numpy()
        a1, a2 = a1.where(~minus, complement(a1)), a2.where(~minus, complement(a2))
    recomputed = allele_sign(a1, a2, mm["ref_a1"], mm["ref_a2"])
    if not np.array_equal(sign, recomputed):
        bad = np.flatnonzero(sign != recomputed)
        raise V2IOError(f"{source}/{locus}: stored sign disagrees with the alleles at {len(bad)} rows, "
                        f"e.g. {s['vkey'].iloc[bad[:3]].tolist()}")
    R = R.astype(dtype, copy=True)
    R *= np.outer(sign, sign).astype(dtype)
    snps = pd.DataFrame({"rsid": mm["rsid"].to_numpy(), "sign": sign,
                         "palindromic": mm["palindromic"].to_numpy()}, index=pd.Index(s["vkey"], name="vkey"))
    return R, snps


def load_ss(anc, locus, ss_root=None, master=None):
    """Sumstats signed to ref_a1, indexed by vkey, in master order.

    Columns: rsid (PGC3's, from the master), z, and whichever of beta, se, freq, n the file
    carries; plus sign and palindromic. Rows outside the window, or matching the master's
    alleles in neither orientation, are dropped. The count is in .attrs["dropped"]."""
    if anc not in ANCESTRIES:
        raise V2IOError(f"unknown ancestry {anc!r}; expected one of {ANCESTRIES}")
    f = os.path.join(_root(ss_root, "SS_ROOT"), f"{anc}_sumstats", f"{locus}_sumstats.tsv.gz")
    d = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
    chrom = locus.split("_")[0].replace("chr", "")
    d["vkey"] = make_vkey(pd.Series(chrom, index=d.index), d["POS"], d["Allele1"], d["Allele2"])
    m = load_master(locus, master)
    n_in = len(d)
    d = d[d["vkey"].isin(m.index)]
    if d["vkey"].duplicated().any():
        raise V2IOError(f"{anc}/{locus}: duplicate vkey in the sumstats, "
                        f"e.g. {d.loc[d['vkey'].duplicated(), 'vkey'].head(3).tolist()}")
    mm = m.loc[d["vkey"]]
    sign = allele_sign(d["Allele1"], d["Allele2"], mm["ref_a1"], mm["ref_a2"])
    keep = sign != 0
    d, mm, sign = d[keep], mm[keep], sign[keep]
    out = pd.DataFrame({"rsid": mm["rsid"].to_numpy(), "z": d["Zscore"].astype(float).to_numpy() * sign},
                       index=pd.Index(d["vkey"].to_numpy(), name="vkey"))
    for col, name in SS_OPTIONAL.items():
        if col in d.columns and name not in out.columns:
            v = pd.to_numeric(d[col], errors="coerce").to_numpy()
            if name == "beta":
                v = v * sign
            elif name == "freq":
                v = np.where(sign == 1, v, 1.0 - v)
            out[name] = v
    out["sign"] = sign
    out["palindromic"] = mm["palindromic"].to_numpy()
    out = out.loc[[v for v in m.index if v in out.index]]  # master order
    out.attrs["dropped"] = {"outside_master": int(n_in - len(keep)), "allele_mismatch": int((~keep).sum())}
    return out


def subset(R, snps, vkeys):
    """R and snps restricted and reordered to vkeys. Every vkey must be present: callers take set
    intersections first, so a missing one is a bug, not a variant to zero-fill."""
    pos = pd.Series(np.arange(len(snps)), index=snps.index)
    missing = [v for v in vkeys if v not in pos.index]
    if missing:
        raise V2IOError(f"{len(missing)} requested vkeys absent from this LD, e.g. {missing[:3]}")
    i = pos.loc[list(vkeys)].to_numpy()
    return R[np.ix_(i, i)], snps.iloc[i]
