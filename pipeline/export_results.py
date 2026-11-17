#!/usr/bin/env python3
"""Pipeline v2, step 7: export a run_finemap.py v2 run in the metric layer's schema.

It reuses export_readers.py's per-arm readers unchanged, with three v2 differences:
- **Alleles come from the master.** Every v2 tool file carries PGC3's rsid, which is unique within
  a locus, and bp/a1/a2 = master bp/ref_a1/ref_a2. v1's merged and EUR arms instead parsed a
  chr:pos:a1:a2 tag, so here export_readers.parse_tag is pointed at the master.
- **common_purity reads the promoted panels through ld_io** (ld_per_anc_eur_gctb +
  ld_per_anc_eas_aou). v1's ld_per_anc_eur/eas are deleted. The composition is export_cs_purity's max-of-mins:
  an ancestry counts only if every member resolves in it.
- **9 arms, no MACHINE.** The MultiSuSiE label is the censored both-LD fit (run_finemap.MULTISUSIE_NPZ).

Outputs, in OUT (default results/finemap/main/ or with_ambiguous/), in the v1 schemas:
  primary_<label>_pip.csv.gz, primary_<label>_cs_purity.csv.gz, primary_status.csv.gz
  <label> = the v1 name, with "-noambig" appended for the noambig tree, so the metric layer can load both.

Usage (environment as for the run):
  HARMONIZED_INPUTS=~/work/v2_harmonized/noambig LD_ROOT=... python export_results.py
L sweep: --L n --arms CARMA exports WORK/CARMA_L{n} as primary_CARMA[-noambig]_L{n}_{pip,cs_purity}.csv.gz, plus
primary_status_L{n}.csv.gz, into OUT (default results/finemap/L_sweep/).
"""
import argparse
import csv
import gzip
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_finemap as D  # noqa: E402
import export_readers as E  # noqa: E402
import ld_io  # noqa: E402
from export_cs_purity import max_of_mins  # noqa: E402

CUR = {}  # rsid -> (bp, ref_a1, ref_a2) for the locus being exported


def parse_tag_v2(tag):
    return CUR.get(tag, ("", "", ""))


def panel(src, L, m):
    """(R, {rsid: index}) for one promoted panel, or (None, None) if the locus has none."""
    try:
        R, s = ld_io.load_ld(src, L)
    except FileNotFoundError:
        return None, None
    rs = m.loc[s.index, "rsid"].to_numpy()
    return R, {r: i for i, r in enumerate(rs)}


def results_dir(variant, L=None, coverage=None):
    """results/finemap/<main | with_ambiguous | L_sweep | coverage_95>[_with_ambiguous]/"""
    if coverage is not None or L is not None:
        sub = ("coverage_95" if coverage is not None else "L_sweep") + ("" if variant == "noambig" else "_with_ambiguous")
    else:
        sub = {"noambig": "main", "primary": "with_ambiguous"}[variant]
    return os.path.join(HERE, "..", "results", "finemap", sub)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=None)
    ap.add_argument("--loci-file", default=os.path.join(HERE, "..", "inputs", "loci_250.txt"))
    ap.add_argument("--L", type=int, default=None, help="L-sweep point: read WORK/<arm>_L{n} (run_finemap --L)")
    ap.add_argument("--coverage", default=None, help="read WORK/<arm>[_L{n}]_cov<NNN> (run_finemap --coverage); "
                    "output goes to results/finemap/coverage_95/ with a _cov095 label suffix")
    ap.add_argument("--arms", default=None, help="comma-separated export labels, e.g. CARMA (default: all nine)")
    ap.add_argument("--eur-ld", default="eur_gctb", choices=["eur_gctb", "eur_aou", "eur_polyfun"],
                    help="EUR panel for common_purity; must match the tree's harmonize --eur-ld (eur_aou)")
    a = ap.parse_args()
    if not D.V2:
        sys.exit("set HARMONIZED_INPUTS to the harmonize.py tree this run used")
    noambig = D.V2_VARIANT == "noambig"
    suffix = ("-noambig" if noambig else "") + (f"_L{a.L}" if a.L is not None else "")
    D.L_OVERRIDE = a.L  # out_dir() then points at WORK/<arm>_L{n}
    D.COV_OVERRIDE = a.coverage
    if a.coverage is not None:
        suffix += D.cov_tag()
        if not a.arms or not set(a.arms.split(",")) <= set(D.COVERAGE_ARMS):
            sys.exit(f"--coverage needs --arms from {list(D.COVERAGE_ARMS)}")
    if a.arms:
        want = set(a.arms.split(","))
        E.ARMS = [(m, f) for m, f in E.ARMS if m in want]
        if not E.ARMS:
            sys.exit(f"no export arm matches {sorted(want)}")
    if a.L is not None and a.coverage is None and any(D.primary_L(E.DRIVER_ARM.get(m, m)) == a.L for m, _ in E.ARMS):
        sys.exit(f"L={a.L} is the primary run for these arms; export the primary run instead")
    out = a.out or results_dir(D.V2_VARIANT, a.L, a.coverage)
    os.makedirs(out, exist_ok=True)
    E.parse_tag = parse_tag_v2
    E.OUT = out  # load_carma_cs writes its temporary dump here
    E.CARMA_CS = E.load_carma_cs() if any(m == "CARMA" for m, _ in E.ARMS) else {}
    loci = [l.strip() for l in open(a.loci_file) if l.strip()]

    def name(method):  # the sweep runs SuSiEx mult-step off, a different fit from SuSiEx-ours (export_readers.label)
        return "SuSiEx-nomultistep" if a.L is not None and method == E.SUSIEX else method

    pw, cw, fh = {}, {}, []
    for method, _ in E.ARMS:
        lab = name(method) + suffix
        f1 = gzip.open(f"{out}/primary_{lab}_pip.csv.gz", "wt", newline="")
        f2 = gzip.open(f"{out}/primary_{lab}_cs_purity.csv.gz", "wt", newline="")
        fh += [f1, f2]
        pw[method], cw[method] = csv.writer(f1), csv.writer(f2)
        pw[method].writerow(["method", "locus", "rsid", "bp", "a1", "a2", "pip"])
        cw[method].writerow(["method", "locus", "rsid", "bp", "a1", "a2", "pip", "cs_ids_all", "n_cs_member",
                             "native_purity", "native_purity_source", "common_purity"])

    status, n_rows, miss, empty, cp_miss, cp_total = [], Counter(), Counter(), [], Counter(), Counter()
    for i, l in enumerate(loci, 1):
        m = ld_io.load_master(l)
        CUR.clear()
        CUR.update({r: (str(bp), a1, a2) for r, bp, a1, a2 in zip(m["rsid"], m["bp"], m["ref_a1"], m["ref_a2"])})
        R_eur, eur_idx = panel(a.eur_ld, l, m)
        R_eas, eas_idx = panel("eas_aou", l, m)
        for method, fn in E.ARMS:
            lab = name(method) + suffix
            rows, cs, native, src, st = fn(l, CUR)
            if not rows:
                empty.append((lab, l))
            member, common = {}, {}
            for k, ms in cs.items():
                for x in ms:
                    member.setdefault(x, []).append(k)
                common[k] = max_of_mins(ms, eur_idx, eas_idx, R_eur, R_eas)[0]
                cp_total[lab] += 1
                cp_miss[lab] += common[k] is None
            for rsid, bp, a1, a2, pip in rows:
                if rsid in CUR:  # the master is the one source of bp/alleles in v2
                    bp, a1, a2 = CUR[rsid]
                else:
                    miss[lab] += 1
                pw[method].writerow([lab, l, rsid, bp, a1, a2, pip])
                ids = member.get(rsid, [])
                nat = [float(native[k]) for k in ids if native.get(k) not in (None, "")]
                com = [common[k] for k in ids if common[k] is not None]
                cw[method].writerow([lab, l, rsid, bp, a1, a2, pip, ";".join(ids), len(ids),
                                     min(nat) if nat else "", src if ids else "", min(com) if com else ""])
            n_rows[lab] += len(rows)
            status.append([lab, l, st.get("n_cs", ""), st.get("converged", ""), st.get("lambda", ""), len(rows),
                           st.get("n_cs_ungated", "")])
        if i % 25 == 0:
            print(f"  {i}/{len(loci)}", flush=True)
    for f in fh:
        f.close()
    with gzip.open(f"{out}/primary_status{f'_L{a.L}' if a.L is not None else ''}{D.cov_tag()}.csv.gz", "wt", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "locus", "n_cs", "converged", "lambda", "n_variants", "n_cs_ungated"])
        w.writerows(status)
    print("rows per arm:", dict(n_rows))
    print("rows whose rsid is not in the master:", dict(miss) or "none")
    print("empty (arm, locus):", empty or "none")
    print("CS with no common_purity / total:", {k: f"{cp_miss[k]}/{cp_total[k]}" for k in cp_total})


if __name__ == "__main__":
    main()
