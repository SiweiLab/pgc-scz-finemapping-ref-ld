#!/usr/bin/env python3
"""B53: EUR reference-panel comparison, three panels (supplement; extends v1's B22).

Same EUR z, three LD panels: GCTB-10k (`eur_gctb`), AoU-267k (`eur_aou`), PolyFun/UKBB-337k
(`eur_polyfun`). Settings:
  - **Common variant set per locus:** EUR sumstats & GCTB & AoU & PolyFun, non-palindromic. z is
    identical across the three fits by construction, so any difference in s or PIP comes from LD
    alone. Per-panel coverage is reported separately (coverage.csv).
  - **Primary EUR-SuSiE-RSS parameters:** run_susie_gctb_aou.R with L_SUSIE_RSS, COVERAGE,
    MIN_PURITY, and n = N_GWAS_PUBLISHED_EUR from run_params. That's the same runner and settings
    as the primary arm; only the domain differs.
  - 249 loci: PolyFun has no chr3_48182866_53553354.
Everything reads through v2_io, and nothing here flips an allele.

  build   OUT/{panel}/{l}.z .ld .master _snps.txt.gz; OUT/s_stats.csv (s, n, min eigenvalue per
          locus x panel; s from ld_diagnostic.estimate_s_scaled_legacy, then named estimate_s_and_kriging, as B22 and gate 2b);
          OUT/coverage.csv
  fit     susie_rss per panel x locus -> OUT/fit/{panel}/{l}/eur_susie_gctb_aou_pip_all.tsv
  export  OUT/export/pip_{panel}.csv.gz (locus, rsid, vkey, bp, pip, in_cs), pip_agreement.csv
          (pairwise r, max |dPIP|, crossings), top_variant.csv

Usage: python eur_panel_compare_b53.py {build,fit,export} --loci-file F --out OUT [--jobs N]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):  # before numpy: build runs --jobs processes (B53)
    os.environ.setdefault(_v, "2")
import argparse
import concurrent.futures as cf
import itertools
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MA_DIR = os.path.dirname(HERE)
FM_DIR = os.path.dirname(MA_DIR)
sys.path.insert(0, HERE)
sys.path.insert(0, MA_DIR)
import run_params as P  # noqa: E402
import v2_io  # noqa: E402
from harmonize_v2 import ordered, write_zld  # noqa: E402
from ld_diagnostic import estimate_s_scaled_legacy  # noqa: E402

PANELS = ("eur_gctb", "eur_aou", "eur_polyfun")
NO_POLYFUN = "chr3_48182866_53553354"
RSCRIPT = os.path.expanduser(os.environ.get("ILAE4_ENV", "~/envs/ilae4_finemap")) + "/bin/Rscript"
SUSIE_R = os.path.join(FM_DIR, "susie", "run_susie_gctb_aou.R")


def build_one(L, out):
    m = v2_io.load_master(L)
    ss = v2_io.load_ss("eur", L)
    ld = {p: v2_io.load_ld(p, L) for p in PANELS}
    nonpal = set(m.index[~m["palindromic"]])
    base = set(ss.index) & nonpal
    common = ordered(m, base.intersection(*[set(ld[p][1].index) for p in PANELS]))
    z = ss.loc[common, "z"].to_numpy(float)
    freq = ss.loc[common, "freq"].to_numpy(float) if "freq" in ss else np.full(len(common), np.nan)
    cov = dict(locus=L, eur_ss_nonpal=len(base), common=len(common),
               **{f"in_{p}": len(base & set(ld[p][1].index)) for p in PANELS})
    srows = []
    rs = m.loc[common, "rsid"].to_numpy()
    for p in PANELS:
        R = v2_io.subset(*ld[p], common)[0].astype(np.float64)
        d = os.path.join(out, p)
        os.makedirs(d, exist_ok=True)
        b = os.path.join(d, L)
        write_zld(b, m, common, z, freq, R)
        pd.DataFrame({"SNP_ID": rs, "MarkerName": rs, "vkey": common}).to_csv(
            f"{b}_snps.txt.gz", sep="\t", index=False, compression="gzip")
        with open(f"{b}.master", "w") as fh:
            fh.write("z;ld;snp;config;cred;log;n_samples\n")
            fh.write(f"{b}.z;{b}.ld;{b}.snp;{b}.config;{b}.cred;{b}.log;{P.N_GWAS_PUBLISHED_EUR}\n")
        s = estimate_s_scaled_legacy(z, R)[0] if len(common) > 1 else np.nan
        srows.append(dict(locus=L, panel=p, n=len(common), s=float(s),
                          min_eig=float(np.linalg.eigvalsh(R).min()) if len(common) else np.nan))
    return cov, srows


def fit_one(p, L, out):
    o = os.path.join(out, "fit", p, L)
    done = os.path.join(o, "eur_susie_gctb_aou_pip_all.tsv")
    if os.path.exists(done):
        return p, L, "skip", ""
    os.makedirs(o, exist_ok=True)
    d = os.path.join(out, p)
    argv = [RSCRIPT, "--vanilla", SUSIE_R, "eur", L, d, f"{d}/{L}_snps.txt.gz", o,
            str(P.L_SUSIE_RSS), P.COVERAGE, P.MIN_PURITY]
    r = subprocess.run(argv, capture_output=True, text=True, env={**os.environ, **P.THREAD_LIMIT_ENV})
    with open(os.path.join(o, "run.log"), "w") as fh:
        fh.write(" ".join(argv) + "\n" + r.stdout + r.stderr)
    ok = r.returncode == 0 and os.path.exists(done)
    return p, L, "ok" if ok else "FAIL", (r.stderr or "")[-300:]


def export(loci, out):
    e = os.path.join(out, "export")
    os.makedirs(e, exist_ok=True)
    pip = {}
    for p in PANELS:
        rows = []
        for L in loci:
            f = os.path.join(out, "fit", p, L, "eur_susie_gctb_aou_pip_all.tsv")
            if not os.path.exists(f):
                continue
            d = pd.read_csv(f, sep="\t", dtype={"MarkerName": str})
            snps = pd.read_csv(os.path.join(out, p, f"{L}_snps.txt.gz"), sep="\t", dtype=str)
            d = d.merge(snps[["MarkerName", "vkey"]], on="MarkerName", how="left")
            rows.append(pd.DataFrame({"locus": L, "rsid": d["MarkerName"], "vkey": d["vkey"],
                                      "bp": d["position"], "pip": d["PIP"], "in_cs": d["in_cs"]}))
        pip[p] = pd.concat(rows, ignore_index=True)
        pip[p].to_csv(os.path.join(e, f"pip_{p}.csv.gz"), index=False)
    agree, top = [], []
    for a, b in itertools.combinations(PANELS, 2):
        j = pip[a].merge(pip[b], on=["locus", "vkey"], suffixes=("_a", "_b"))
        dd = (j["pip_a"] - j["pip_b"]).abs()
        k = dd.idxmax()
        agree.append(dict(pair=f"{a}~{b}", loci=j["locus"].nunique(), n=len(j),
                          pearson_r=round(float(np.corrcoef(j["pip_a"], j["pip_b"])[0, 1]), 4),
                          max_abs_dpip=round(float(dd.max()), 4), max_at=f"{j.at[k, 'locus']}|{j.at[k, 'rsid_a']}",
                          **{f"crossing_{t}": int(((j["pip_a"] > t) != (j["pip_b"] > t)).sum()) for t in (0.1, 0.5, 0.95)}))
    for L in loci:
        r = {"locus": L}
        for p in PANELS:
            d = pip[p][pip[p]["locus"] == L]
            if len(d):
                i = d["pip"].idxmax()
                r[f"top_{p}"], r[f"pip_{p}"] = d.at[i, "rsid"], round(float(d.at[i, "pip"]), 4)
        top.append(r)
    pd.DataFrame(agree).to_csv(os.path.join(e, "pip_agreement.csv"), index=False)
    t = pd.DataFrame(top)
    t["all_agree"] = t[[f"top_{p}" for p in PANELS]].nunique(axis=1) == 1
    t.to_csv(os.path.join(e, "top_variant.csv"), index=False)
    print(pd.DataFrame(agree).to_string(index=False))
    print(f"top variant identical across all three panels at {int(t['all_agree'].sum())}/{len(t)} loci")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["build", "fit", "export"])
    ap.add_argument("--loci-file", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    out = os.path.expanduser(a.out)
    loci = [l.strip() for l in open(a.loci_file) if l.strip() and l.strip() != NO_POLYFUN]
    if a.step == "build":
        covs, srows = [], []
        with cf.ProcessPoolExecutor(a.jobs) as ex:
            for cov, s in ex.map(build_one, loci, itertools.repeat(out)):
                covs.append(cov)
                srows += s
                print(cov["locus"], cov["common"], [round(x["s"], 4) for x in s], flush=True)
        pd.DataFrame(covs).to_csv(os.path.join(out, "coverage.csv"), index=False)
        S = pd.DataFrame(srows)
        S.to_csv(os.path.join(out, "s_stats.csv"), index=False)
        print(S.groupby("panel")["s"].describe().to_string())
    elif a.step == "fit":
        jobs = [(p, L) for p in PANELS for L in loci]
        with cf.ThreadPoolExecutor(a.jobs) as ex:
            for p, L, st, msg in ex.map(lambda x: fit_one(x[0], x[1], out), jobs):
                if st == "FAIL":
                    print("FAIL", p, L, msg, flush=True)
        print("fit done")
    else:
        export(loci, out)


if __name__ == "__main__":
    main()
