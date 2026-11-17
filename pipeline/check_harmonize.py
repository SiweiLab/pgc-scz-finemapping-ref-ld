#!/usr/bin/env python3
"""Gates 4 and 6 of pipeline v2, over harmonize.py's output.

Gate 4 (merged LD, a check, nothing built): over v1's merged domain, the run-time merge of the
promoted EUR and EAS LD must reproduce v1 `ld_eur_eas_merged/{l}.ld`. v1's weights aren't stored,
so the check fits the one free parameter, alpha = w_EUR / (w_EUR + w_EAS), by least squares over
the off-diagonal. It reports alpha_fit against v2's run-time alpha, the max |residual| at alpha_fit
(expect <= 5e-7, v1's %.6f text), and every row whose residual exceeds 1e-5. Those rows are
v2's EUR or EAS side changing, which is expected only at a few known rows.

Gate 6 (per-arm inputs):
  a. Every domain equals its set expression, recomputed independently of harmonize.py: from
     the master's in_{anc}_ss flags and the promoted snps.tsv.gz vkeys (valid = flag & in LD).
  b. Each tool file's row count equals its logged domain.
  c. Against v1's exports (the original repository's git history; skipped where absent): per arm and locus, the v1 variants
     v2 lacks, and the v2 variants v1 lacked. The domain of a v1 arm is its PIP export's variant set.
  d. Only ld_io.py flips an allele: grep harmonize.py for sign arithmetic.

Usage:
  python check_harmonize.py --run-dir RUN --v1-merged DIR --ld-root D --loci-file F --out OUT [--noambig]
"""
import argparse
import gzip
import json
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ld_io  # noqa: E402

# v1 exports, read from the original working repository's git history (not part of this release; skipped when absent)
V1_EXPORT = {False: ("b59d67a", "task_b30_out/primary_{arm}_pip.csv.gz", ""),
             True: ("b59d67a", "task_b34_out/primary_{arm}-noambig_pip.csv.gz", "-noambig")}
ARM_DOMAIN = {"MultiSuSiE": "multisusie", "MESuSiE": "mesusie", "SuSiE-merged": "merged",
              "FINEMAP-merged": "merged", "CARMA": "merged", "RSparsePro": "merged",
              "EUR-SuSiE-RSS": "eur", "EUR-FINEMAP": "eur"}


def vkey_from(chrom, bp, a1, a2):
    lo, hi = sorted([str(a1).upper(), str(a2).upper()])
    return f"{str(chrom).replace('chr', '')}:{bp}:{lo}:{hi}"


def v1_export(arm, noambig):
    ref, tpl, _ = V1_EXPORT[noambig]
    path = "Post-hoc_analyses/fine-mapping/multi_ancestry/aou_ld/" + tpl.format(arm=arm)
    root = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True,
                          cwd=HERE).stdout.strip()
    blob = subprocess.run(["git", "show", f"{ref}:{path}"], capture_output=True, cwd=root)
    if blob.returncode:
        return None
    import io
    d = pd.read_csv(io.BytesIO(gzip.decompress(blob.stdout)), dtype=str)
    d["vkey"] = [vkey_from(L.split("_")[0], bp, x, y) for L, bp, x, y in zip(d["locus"], d["bp"], d["a1"], d["a2"])]
    return d.groupby("locus")["vkey"].apply(set).to_dict()


def gate4(L, run, v1_dir, ld_root, log):
    zf, lf = os.path.join(v1_dir, f"{L}.z"), os.path.join(v1_dir, f"{L}.ld")
    if not (os.path.exists(zf) and os.path.exists(lf)):
        return dict(locus=L, status="v1_missing")
    z1 = pd.read_csv(zf, sep=" ", dtype=str)
    R1 = np.loadtxt(lf)
    m = ld_io.load_master(L)
    vk = [vkey_from(c, p, x, y) for c, p, x, y in zip(z1["chromosome"], z1["position"], z1["allele1"], z1["allele2"])]
    ok = np.array([v in m.index for v in vk])
    s1 = np.array([(1 if (a == m.at[v, "ref_a1"]) else -1) if o else 0
                   for v, a, o in zip(vk, z1["allele1"].str.upper(), ok)])
    Re, se = ld_io.load_ld("eur_gctb", L, ld_root=ld_root)
    Ra, sa = ld_io.load_ld("eas_aou", L, ld_root=ld_root)
    keep = np.array([o and v in se.index and v in sa.index for v, o in zip(vk, ok)])
    idx = np.flatnonzero(keep)
    shared = [vk[i] for i in idx]
    A = (R1 * np.outer(s1, s1))[np.ix_(idx, idx)]
    E = ld_io.subset(Re, se, shared)[0].astype(np.float64)
    S = ld_io.subset(Ra, sa, shared)[0].astype(np.float64)
    off = ~np.eye(len(idx), dtype=bool)
    dE, dY = (E - S)[off], (A - S)[off]
    alpha = float(dE @ dY / (dE @ dE)) if dE @ dE > 0 else float("nan")
    res = np.abs(A - (alpha * E + (1 - alpha) * S))
    np.fill_diagonal(res, 0)
    w = log.get("merge_weights", {})
    alpha_v2 = w["EUR"] / (w["EUR"] + w["EAS"]) if w else float("nan")
    res_v2 = np.abs(A - (alpha_v2 * E + (1 - alpha_v2) * S))
    np.fill_diagonal(res_v2, 0)
    bad_rows = np.flatnonzero(res.max(axis=1) > 1e-5)
    return dict(locus=L, status="ok", v1_rows=len(vk), compared=len(idx), v1_not_in_v2=int((~keep).sum()),
                alpha_fit=alpha, alpha_v2=alpha_v2, max_res_fit=float(res.max()) if res.size else 0.0,
                max_res_v2_alpha=float(res_v2.max()) if res_v2.size else 0.0,
                rows_res_gt_1e5=";".join(m.at[shared[i], "rsid"] for i in bad_rows[:20]), n_rows_res_gt_1e5=len(bad_rows))


def gate6(L, run, ld_root, log, noambig, v1_sets):
    m = ld_io.load_master(L)
    allow = ~m["palindromic"] if noambig else pd.Series(True, index=m.index)
    ldk = {}
    for src, a in (("eur_gctb", "EUR"), ("eas_aou", "EAS")):
        s = pd.read_csv(os.path.join(ld_root, f"ld_per_anc_{src}", L, "snps.tsv.gz"), sep="\t", dtype=str)
        ldk[a] = set(s["vkey"])
    ok = set(m.index[allow])
    flag = {a: set(m.index[m[f"in_{a.lower()}_ss"]]) for a in ("EUR", "EAS")}
    core = set(m.index[m["in_core_ss"]])
    valid = {a: flag[a] & ldk[a] & ok for a in ("EUR", "EAS")}
    expect = {"susiex_EUR": len(valid["EUR"]), "susiex_EAS": len(valid["EAS"]),
              "multisusie": len((valid["EUR"] | valid["EAS"]) & ldk["EUR"] & ldk["EAS"] & ok),
              "union_censored_bothld": len((valid["EUR"] | valid["EAS"]) & ldk["EUR"] & ldk["EAS"] & ok),
              "union_censored": len(valid["EUR"] | valid["EAS"]),
              "mesusie": len(valid["EUR"] & valid["EAS"]),
              "merged": len(ldk["EUR"] & ldk["EAS"] & core & ok), "eur": len(valid["EUR"])}
    got = log["domains"]
    r = dict(locus=L)
    for k, e in expect.items():
        r[f"{k}_expect"], r[f"{k}_got"] = e, got.get(k)
    # b. file row counts
    def nlines(p, header=True):
        op = gzip.open if p.endswith(".gz") else open
        with op(p, "rt") as f:
            return sum(1 for _ in f) - (1 if header else 0)
    r["file_merged_z"] = nlines(os.path.join(run, "merged", f"{L}.z"))
    r["file_eur_z"] = nlines(os.path.join(run, "eur", f"{L}.z"))
    r["file_mesusie_ma"] = nlines(os.path.join(run, "intersect", f"{L}.EUR.ma"))
    r["file_susiex_eas"] = nlines(os.path.join(run, "union", f"{L}_susiex", "EAS.sst.txt"))
    r["file_npz"] = len(np.load(os.path.join(run, "union", f"{L}_harmonized.npz"), allow_pickle=True)["rsids"])
    r["file_npz_cens"] = len(np.load(os.path.join(run, "union_censored", f"{L}_harmonized.npz"), allow_pickle=True)["rsids"])
    # c. against v1 exports
    dom_sets = {"merged": set(pd.read_csv(os.path.join(run, "merged", f"{L}.z"), sep=" ", dtype=str)["rsid"]),
                "eur": set(pd.read_csv(os.path.join(run, "eur", f"{L}.z"), sep=" ", dtype=str)["rsid"]),
                "mesusie": set(pd.read_csv(os.path.join(run, "intersect", f"{L}.EUR.ma"), sep="\t", dtype=str)["SNP"]),
                "multisusie": set(np.load(os.path.join(run, "union", f"{L}_harmonized.npz"), allow_pickle=True)["rsids"])}
    rs2vk = m.reset_index().set_index("rsid")["vkey"]
    for arm, dom in ARM_DOMAIN.items():
        v1 = v1_sets.get(arm, {}).get(L)
        if v1 is None:
            continue
        v2 = {rs2vk[x] for x in dom_sets[dom] if x in rs2vk.index}
        r[f"{arm}_v1_only"], r[f"{arm}_v2_only"] = len(v1 - v2), len(v2 - v1)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--v1-merged", required=True, help="local copy of v1 ld_eur_eas_merged/ (or _noambig)")
    ap.add_argument("--ld-root", required=True)
    ap.add_argument("--loci-file", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--noambig", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    loci = [l.strip() for l in open(a.loci_file) if l.strip()]
    v1_sets = {arm: s for arm in ARM_DOMAIN if (s := v1_export(arm, a.noambig)) is not None}
    print("v1 exports found:", sorted(v1_sets), flush=True)
    g4, g6 = [], []
    for L in loci:
        with open(os.path.join(a.run_dir, "_logs", f"{L}.json")) as f:
            log = json.load(f)
        g4.append(gate4(L, a.run_dir, a.v1_merged, a.ld_root, log))
        g6.append(gate6(L, a.run_dir, a.ld_root, log, a.noambig, v1_sets))
        print(L, g4[-1].get("max_res_fit"), flush=True)
    g4, g6 = pd.DataFrame(g4), pd.DataFrame(g6)
    tag = "noambig" if a.noambig else "primary"
    g4.to_csv(os.path.join(a.out, f"gate4_{tag}.csv"), index=False)
    g6.to_csv(os.path.join(a.out, f"gate6_{tag}.csv"), index=False)
    # d. grep: nothing but ld_io flips an allele
    src = open(os.path.join(HERE, "harmonize.py")).read()
    flips = [l for l in src.splitlines() if re.search(r"\*\s*-1|\*=\s*-1|np\.outer\(sign|1\.0\s*-\s*\w*freq", l)
             and not l.strip().startswith("#")]
    mism = [k for k in ("susiex_EUR", "susiex_EAS", "multisusie", "union_censored_bothld", "union_censored",
                        "mesusie", "merged", "eur")
            if not (g6[f"{k}_expect"] == g6[f"{k}_got"]).all()]
    fmis = [c for c, k in (("file_merged_z", "merged_got"), ("file_eur_z", "eur_got"), ("file_mesusie_ma", "mesusie_got"),
                           ("file_susiex_eas", "susiex_EAS_got"), ("file_npz", "multisusie_got"),
                           ("file_npz_cens", "union_censored_got")) if not (g6[c] == g6[k]).all()]
    ok4 = g4[g4["status"] == "ok"]
    print(f"GATE4 {tag}: loci {len(ok4)}/{len(g4)}; max residual at fitted alpha {ok4['max_res_fit'].max():.3g}; "
          f"alpha_fit {ok4['alpha_fit'].min():.4f}-{ok4['alpha_fit'].max():.4f} vs v2 "
          f"{ok4['alpha_v2'].min():.4f}-{ok4['alpha_v2'].max():.4f}; loci with rows > 1e-5: "
          f"{int((ok4['n_rows_res_gt_1e5'] > 0).sum())}")
    print(f"GATE6 {tag}: set-expression mismatches {mism or 'none'}; file-count mismatches {fmis or 'none'}; "
          f"grep flips in harmonize.py: {flips or 'none'}")
    for k in ("susiex_EUR", "susiex_EAS", "multisusie", "union_censored_bothld", "union_censored",
              "mesusie", "merged", "eur"):
        print(f"  domain {k}: {int(g6[f'{k}_got'].sum()):,}")
    for arm in v1_sets:
        if f"{arm}_v1_only" in g6:
            print(f"  vs v1 {arm}: v1-only {int(g6[f'{arm}_v1_only'].sum()):,}, v2-only {int(g6[f'{arm}_v2_only'].sum()):,}")


if __name__ == "__main__":
    main()
