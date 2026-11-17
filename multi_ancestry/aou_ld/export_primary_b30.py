#!/usr/bin/env python3
"""Task B30 Step 3: export the 2026-09 primary run (run_primary_2026.py, 10 arms x 250 loci)
in the metric layer's schema.

Outputs (task_b30_out/), one pair per arm:
  primary_<arm>_pip.csv.gz        Task-A schema: method, locus, rsid, bp, a1, a2, pip.
                                  load_pips() reads it as-is.
  primary_<arm>_cs_purity.csv.gz  Task-B12 schema: ... pip, cs_ids_all, n_cs_member,
                                  native_purity, native_purity_source, common_purity. Every
                                  arm that defines credible sets (all but MACHINE).
  primary_status.csv.gz           method, locus, n_cs, converged, lambda, n_variants,
                                  n_cs_ungated (SuSiEx only)

L sweep (B31): `--L n` exports the 7 sweep arms (run_primary_2026.SWEEP_L) at that L into
task_b31_out/ as primary_<label>_L<n>_{pip,cs_purity}.csv.gz and primary_status_L<n>.csv.gz,
with `method` = <label>_L<n>.

No-ambiguous rerun (B33): `--noambig`, with PRIMARY_WORK=~/work/primary_2026_noambig and
PRIMARY_UNION=<inputs>/ld_eur_eas_union_noambig in the environment, exports SuSiEx alone into
task_b33_out/ as primary_SuSiEx-noambig_{pip,cs_purity}.csv.gz and primary_status.csv.gz,
`method` = SuSiEx-noambig. Loci come from the SuSiEx outputs. B34: `--noambig --arms
<driver arms> --out-dir task_b34_out` exports the other arms of the same WORK as
<arm>-noambig; the allele backfill is that WORK's MultiSuSiE table when it exists (as in the
primary export), else empty (B33 ran SuSiEx alone; susiex() never uses it). A separate directory on purpose: the laptop's loaders glob
primary_*_cs_purity.csv.gz, so sweep files next to the B30 export would enter the primary
figures as extra arms. Where the primary run already IS the point (run_primary_2026.primary_L),
it is read from WORK/<arm> -- re-exported, not recomputed -- and the status file's `source`
column says so. SuSiEx's sweep ran mult-step off, a different setting from the primary
SuSiEx-ours arm, so it is labelled SuSiEx-nomultistep and never reads the primary output.

The two exceptions the brief names:
  * MACHINE returns fit@CL as an M x K matrix, one column per ancestry. Its pip file carries the
    standard columns with `pip` left EMPTY, plus CL_EUR and CL_EAS: the runner's pooled `CL`
    column is a row max, and how to aggregate is a judgement call. MACHINE's fit@CS was never
    written to disk (run_machine_diag_b29.R keeps only its length), so it has no cs_purity file.
  * CARMA: membership from `Credible set` (rho-based), not `Credible model`
    (PRIMARY_RUN_PLAN 4.2). Indices are 1-based rows of the arm's own z file, which is also
    the row order of {l}_carma_pip.tsv.

SuSiEx ran at --min_purity 0 (PRIMARY_RUN_PLAN 4.3). PIP is the full posterior -- the max over
every PIP(CS*) column, the same reduction export_full_pip_vectors.py used -- and is NOT gated.
Credible sets ARE gated: a component is kept only if its CS_PURITY >= run_params.MIN_PURITY.
n_cs_ungated records how many components existed before the gate.

native_purity is filled only where the tool itself wrote one to disk (SuSiEx CS_PURITY,
MultiSuSiE min_abs_corr, FINEMAP #min(|ld|)); purity on FINEMAP x2 / CARMA / RSparsePro is
accepted as unalignable (PRIMARY_RUN_PLAN, signed off 2026-09-23). common_purity is B12's:
max over ancestries of the min pairwise |r| among CS members, on one shared pair of panels for
every arm -- ld_per_anc_eur + ld_per_anc_eas, the same bytes B12 read from b6b7 -- with B12's
rule that an ancestry contributes only if EVERY member resolves in it; blank if neither does.
Members are resolved to the panels' MarkerName through (bp, unordered alleles) against the
SuSiEx sst files, which list every union variant with its rsid.

Identifier space per arm is each arm's own (variant_key in the metric layer joins across them):
  SuSiEx, MultiSuSiE, MESuSiE, MACHINE: dbSNP rsid (union / intersection trees). MESuSiE and
    MACHINE take bp/a1/a2 from the MultiSuSiE pip table of the same union input; SuSiEx sees
    the wider union of its EUR/EAS sst files, so it takes A1/A2 from those, its own inputs.
  SuSiE-merged, FINEMAP-merged, RSparsePro, CARMA, EUR-FINEMAP: chr:pos:a1:a2 tag.
  EUR-SuSiE-RSS: MarkerName (dbSNP rsid), bp/a1/a2 parsed from its own SNP_ID tag.
"""
import csv
import gzip
import os
import re
import argparse
import subprocess
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_params import MIN_PURITY  # noqa: E402
import run_primary_2026 as D  # noqa: E402
from run_primary_2026 import INPUTS, RSCRIPT, UNION, WORK, out_dir  # noqa: E402
from export_cs_purity_b12 import max_of_mins  # noqa: E402

MIN_PURITY_F = float(MIN_PURITY)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "task_b30_out")
OUT_SWEEP = os.path.join(HERE, "task_b31_out")
OUT_NOAMBIG = os.path.join(HERE, "task_b33_out")
NOAMBIG = False
SUFFIX = "-noambig"
MACHINE_ALT = False
TAG_RE = re.compile(r"^chr[\dXYM]+:(\d+):([ACGT]+):([ACGT]+)$", re.IGNORECASE)

EUR_LD = f"{INPUTS}/ld_per_anc_eur"
EAS_LD = f"{INPUTS}/ld_per_anc_eas"

SUSIEX = "SuSiEx-ours"   # method label used by Task A / B12 for our SuSiEx runs


def read_tsv(path, delim="\t"):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        if delim == " ":
            return list(csv.DictReader(f, delimiter=" ", skipinitialspace=True))
        return list(csv.DictReader(f, delimiter=delim))


def parse_tag(tag):
    m = TAG_RE.match(tag)
    return (m.group(1), m.group(2), m.group(3)) if m else ("", "", "")


def loci_list(driver_arms=("SuSiEx",)):
    if NOAMBIG:
        # loci every requested arm has an output for (B37 reruns two loci for some arms)
        cand = sorted(n[:-len("_harmonized.npz")] for n in os.listdir(UNION)
                      if n.endswith("_harmonized.npz"))
        return [l for l in cand if all(os.path.exists(D.ARMS[a][3].format(o=out_dir(a), l=l))
                                       for a in driver_arms)]
    d = out_dir("MultiSuSiE")
    return sorted(n[:-len("_multisusie_pip.tsv")] for n in os.listdir(d)
                  if n.endswith("_multisusie_pip.tsv"))


# ---------------------------------------------------------------------------
# One function per arm. Each returns (rows, cs_members, native, native_source, status):
#   rows        [(rsid, bp, a1, a2, pip)] in the tool's own order
#   cs_members  {cs_id: [rsid, ...]}
#   native      {cs_id: purity}   (only where the tool wrote one)
#   status      dict of status-file fields
# ---------------------------------------------------------------------------

def multisusie(l, ref):
    rows = [(r["rsid"], r["pos"], r["effect_allele"], r["other_allele"], r["pip"])
            for r in read_tsv(f"{out_dir('MultiSuSiE')}/{l}_multisusie_pip.tsv")]
    cs, native = {}, {}
    for r in read_tsv(f"{out_dir('MultiSuSiE')}/{l}_multisusie_cs.tsv"):
        cs.setdefault(r["cs_id"], []).append(r["rsid"])
        native[r["cs_id"]] = r["min_abs_corr"]
    return rows, cs, native, "on_disk", {"n_cs": len(cs)}


def mesusie(l, ref):
    rows = [(r["rsid"], *ref.get(r["rsid"], ("", "", "")), r["pip"])
            for r in read_tsv(f"{out_dir('MESuSiE')}/{l}_mesusie_pip.tsv")]
    cs = {}
    for r in read_tsv(f"{out_dir('MESuSiE')}/{l}_mesusie_config.tsv"):
        cs.setdefault(r["cs_id"], []).append(r["rsid"])
    return rows, cs, {}, "", {"n_cs": len(cs)}


def susiex(l, ref):
    # SuSiEx sees the union of its EUR and EAS sst files, wider than the union npz the
    # MultiSuSiE ref covers, so alleles come from its own inputs (EUR first, then EAS).
    own = {}
    for anc in ("EAS", "EUR"):
        for r in read_tsv(f"{UNION}/{l}_susiex/{anc}.sst.txt"):
            own[r["snp"]] = (r["A1"], r["A2"])
    rows = []
    for r in read_tsv(f"{out_dir('SuSiEx')}/{l}_susiex.snp"):
        pip = max((float(r[c]) for c in r if c.startswith("PIP(CS")), default=0.0)
        a1, a2 = own.get(r["SNP"], ("", ""))
        rows.append((r["SNP"], r["BP"], a1, a2, pip))
    purity = {}
    with open(f"{out_dir('SuSiEx')}/{l}_susiex.summary") as f:
        lines = f.readlines()
    if len(lines) > 1 and lines[1].strip() != "NULL":
        for r in csv.DictReader(lines[1:], delimiter="\t"):
            purity[r["CS_ID"]] = float(r["CS_PURITY"])
    keep = {k for k, v in purity.items() if v >= MIN_PURITY_F}
    cs = {}
    for r in read_tsv(f"{out_dir('SuSiEx')}/{l}_susiex.cs"):
        if r["CS_ID"] in keep:
            cs.setdefault(r["CS_ID"], []).append(r["SNP"])
    native = {k: purity[k] for k in keep}
    return rows, cs, native, "on_disk", {"n_cs": len(keep), "n_cs_ungated": len(purity)}


def _in_cs_arm(pip_path, id_col, pip_col, tag_col):
    rows, cs = [], {}
    for r in read_tsv(pip_path):
        bp, a1, a2 = parse_tag(r[tag_col])
        rows.append((r[id_col], bp or r["position"], a1, a2, r[pip_col]))
        if r["in_cs"] != ".":
            cs.setdefault(r["in_cs"], []).append(r[id_col])
    return rows, cs


def _summary(path):
    s = read_tsv(path)
    return ({"n_cs": s[0]["n_cs"], "converged": s[0]["converged"], "lambda": s[0]["lambda"]}
            if s else {})


def susie_merged(l, ref):
    o = out_dir("SuSiE-merged")
    rows, cs = _in_cs_arm(f"{o}/{l}_susie_pip.tsv", "rsid", "PIP", "rsid")
    return rows, cs, {}, "", _summary(f"{o}/{l}_susie_summary.tsv")


def eur_susie(l, ref):
    o = f"{out_dir('EUR-SuSiE-RSS')}/{l}"
    rows, cs = _in_cs_arm(f"{o}/eur_susie_gctb_aou_pip_all.tsv", "MarkerName", "PIP", "SNP_ID")
    return rows, cs, {}, "", _summary(f"{o}/eur_susie_gctb_aou_summary.tsv")


def finemap_kstar(log_path):
    """argmax_k Post-Pr(# causal = k), k >= 1 -- export_full_pip_vectors.finemap_n_cs's rule."""
    if not os.path.exists(log_path):
        return None
    text = open(log_path).read()
    num = r"[\d.]+(?:[eE][+-]?\d+)?"
    m = re.search(r"Post-Pr\(# of causal SNPs is k\)\s*:\s*\n((?:\s*\(?\d+\)?\s*->\s*"
                  + num + r"\)?\s*\n?)+)", text)
    if not m:
        return None
    k = {}
    for line in m.group(1).strip().splitlines():
        mm = re.match(r"\s*\(?(\d+)\)?\s*->\s*(" + num + r")", line)
        if mm and int(mm.group(1)) > 0:
            k[int(mm.group(1))] = float(mm.group(2))
    return max(k, key=lambda x: (k[x], -x)) if k else None


def finemap(arm):
    def f(l, ref):
        o = out_dir(arm)
        rows = [(r["rsid"], r["position"], r["allele1"], r["allele2"], r["prob"])
                for r in read_tsv(f"{o}/{l}.snp", " ")]
        k = finemap_kstar(f"{o}/{l}.log_sss")
        cs, native = {}, {}
        if k:
            with open(f"{o}/{l}.cred{k}") as fh:
                lines = fh.readlines()
            minld = next((ln.split()[1:] for ln in lines if ln.startswith("#min(|ld|)")), [])
            body = [ln for ln in lines if not ln.startswith("#")]
            for r in csv.DictReader(body, delimiter=" ", skipinitialspace=True):
                for j in range(1, k + 1):
                    v = r.get(f"cred{j}")
                    if v and v != "NA":
                        cs.setdefault(str(j), []).append(v)
            for j in range(1, k + 1):          # "#min(|ld|) v1 NA v2 NA": value per cred column
                if 2 * (j - 1) < len(minld) and minld[2 * (j - 1)] != "NA":
                    native[str(j)] = minld[2 * (j - 1)]
        return rows, cs, native, "on_disk", {"n_cs": k if k else ""}
    return f


def rsparsepro(l, ref):
    rows, cs = [], {}
    for r in read_tsv(f"{out_dir('RSparsePro')}/{l}.rsparsepro.txt"):
        bp, a1, a2 = parse_tag(r["RSID"])
        rows.append((r["RSID"], bp, a1, a2, r["PIP"]))
        if r["cs"] != "0":
            cs.setdefault(r["cs"], []).append(r["RSID"])
    return rows, cs, {}, "", {"n_cs": len(cs)}


CARMA_CS = None


def load_carma_cs():
    """One Rscript call over every locus: 1-based row indices of each `Credible set` member."""
    dump = os.path.join(OUT, f"_carma_credible_set_{os.getpid()}.tsv")
    r = f'''
    d <- "{out_dir("CARMA")}"
    fs <- list.files(d, pattern = "_carma[.]rds$")
    out <- list()
    for (f in fs) {{
      locus <- sub("_carma[.]rds$", "", f)
      cs <- readRDS(file.path(d, f))[[1]][["Credible set"]][[2]]
      if (length(cs)) for (i in seq_along(cs))
        out[[length(out) + 1]] <- data.frame(locus = locus, cs = i, idx = cs[[i]])
    }}
    res <- if (length(out)) do.call(rbind, out) else data.frame(locus = character(), cs = integer(), idx = integer())
    write.table(res, "{dump}", sep = "\\t", row.names = FALSE, quote = FALSE)
    cat(length(fs), "rds read\\n")
    '''
    p = subprocess.run([RSCRIPT, "-e", r], capture_output=True, text=True)
    if p.returncode:
        sys.exit(f"CARMA credible-set dump failed:\n{p.stderr}")
    print("CARMA:", p.stdout.strip())
    cs = {}
    for row in read_tsv(dump):
        cs.setdefault(row["locus"], {}).setdefault(row["cs"], []).append(int(row["idx"]))
    os.remove(dump)
    return cs


def carma(l, ref):
    rows = []
    for r in read_tsv(f"{out_dir('CARMA')}/{l}_carma_pip.tsv"):
        _, a1, a2 = parse_tag(r["rsid"])
        rows.append((r["rsid"], r["bp"], a1, a2, r["pip"]))
    cs = {k: [rows[i - 1][0] for i in idx] for k, idx in CARMA_CS.get(l, {}).items()}
    return rows, cs, {}, "", {"n_cs": len(cs)}


ARMS = [(SUSIEX, susiex), ("MultiSuSiE", multisusie), ("MESuSiE", mesusie),
        ("SuSiE-merged", susie_merged), ("FINEMAP-merged", finemap("FINEMAP-merged")),
        ("RSparsePro", rsparsepro), ("CARMA", carma),
        ("EUR-SuSiE-RSS", eur_susie), ("EUR-FINEMAP", finemap("EUR-FINEMAP"))]


def unordered_key(bp, a1, a2):
    lo, hi = sorted((a1.upper(), a2.upper()))
    return f"{int(float(bp))}:{lo}:{hi}"


def load_panel(ld_dir, l):
    """(R, {MarkerName: index}) for one ancestry's per-ancestry panel, or (None, None)."""
    r_path, s_path = f"{ld_dir}/{l}_R.npz", f"{ld_dir}/{l}_snps.txt.gz"
    if not (os.path.exists(r_path) and os.path.exists(s_path)):
        return None, None
    snps = pd.read_csv(s_path, sep="\t", usecols=["MarkerName"])
    return np.load(r_path)["R"], {m: i for i, m in enumerate(snps.MarkerName)}


def common_purity_fn(l):
    """cs members (any arm's own ids) -> max-of-mins on the shared per-ancestry panels."""
    R_eur, eur_idx = load_panel(EUR_LD, l)
    R_eas, eas_idx = load_panel(EAS_LD, l)
    to_rsid = {}
    for anc in ("EAS", "EUR"):
        for r in read_tsv(f"{UNION}/{l}_susiex/{anc}.sst.txt"):
            to_rsid[unordered_key(r["bp"], r["A1"], r["A2"])] = r["snp"]
            to_rsid[r["snp"]] = r["snp"]

    def resolve(m, bp, a1, a2):
        return to_rsid.get(m) or (to_rsid.get(unordered_key(bp, a1, a2)) if bp and a1 and a2 else None)

    def purity(members):
        if any(m is None for m in members):
            return None
        return max_of_mins(members, eur_idx, eas_idx, R_eur, R_eas)[0]
    return resolve, purity


DRIVER_ARM = {SUSIEX: "SuSiEx"}          # export label -> run_primary_2026 arm name


def use_dir(method, L):
    """Point out_dir() at this (arm, L)'s outputs: the sweep dir, or the primary dir when the
    primary run already is that point. Returns 'run' or 'reused-primary' for provenance."""
    arm = DRIVER_ARM.get(method, method)
    if L is None or D.primary_L(arm) == L:
        D.L_OVERRIDE = None
        return "primary" if L is None else "reused-primary"
    D.L_OVERRIDE = L
    return "run"


def label(method, L):
    if NOAMBIG:
        return f"{DRIVER_ARM.get(method, method)}{SUFFIX}"
    if L is None:
        return method
    return f"{'SuSiEx-nomultistep' if method == SUSIEX else method}_L{L}"


def main():
    global CARMA_CS, OUT, NOAMBIG, SUFFIX, MACHINE_ALT
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, default=None, help="export one L-sweep point (B31)")
    ap.add_argument("--arms", default=None, help="with --L: comma-separated driver arm names")
    ap.add_argument("--noambig", action="store_true",
                    help="B33/B34: <arm>-noambig labels; SuSiEx alone into task_b33_out/ "
                         "unless --arms / --out-dir say otherwise")
    ap.add_argument("--out-dir", default=None, help="with --noambig: output directory")
    ap.add_argument("--suffix", default="-noambig",
                    help="with --noambig: label suffix (B37: -orientfix, -noambig-orientfix)")
    a_ = ap.parse_args()
    L = a_.L
    NOAMBIG = a_.noambig
    SUFFIX = a_.suffix
    if NOAMBIG:
        if L is not None:
            sys.exit("--noambig does not take --L")
        if a_.arms and a_.arms != "SuSiEx" and not a_.out_dir:
            sys.exit("--noambig with arms other than SuSiEx needs --out-dir (not task_b33_out)")
        if not os.path.basename(os.path.normpath(UNION)).startswith("ld_eur_eas_union_"):
            sys.exit(f"--noambig needs PRIMARY_UNION at an ld_eur_eas_union_* rebuild; got {UNION}")
        if os.path.normpath(WORK) == os.path.normpath(os.path.expanduser("~/work/primary_2026")):
            sys.exit("--noambig needs PRIMARY_WORK pointed at the no-ambiguous run")
    arms = ARMS if L is None else [(m, f) for m, f in ARMS if DRIVER_ARM.get(m, m) in D.SWEEP_L]
    if a_.arms and not NOAMBIG:
        want = set(a_.arms.split(","))
        arms = [(m, f) for m, f in arms if DRIVER_ARM.get(m, m) in want]
    if NOAMBIG:
        want = set((a_.arms or "SuSiEx").split(","))
        MACHINE_ALT = "MACHINE" in want
        arms = [(m, f) for m, f in ARMS if DRIVER_ARM.get(m, m) in want]
        if len(arms) + MACHINE_ALT != len(want):
            sys.exit(f"unknown arm(s) in {want}")
        OUT = os.path.abspath(a_.out_dir) if a_.out_dir else OUT_NOAMBIG
    if L is not None:
        OUT = OUT_SWEEP
    os.makedirs(OUT, exist_ok=True)
    D.L_OVERRIDE = None
    loci = loci_list(sorted(want)) if NOAMBIG else loci_list()
    print(f"{len(loci)} loci under {WORK}; L = {L}", flush=True)
    source = {m: use_dir(m, L) for m, _ in arms}
    if any(m == "CARMA" for m, _ in arms):
        use_dir("CARMA", L)
        CARMA_CS = load_carma_cs()
    for m, _ in arms:
        use_dir(m, L)
        print(f"  {label(m, L):<26} <- {out_dir(DRIVER_ARM.get(m, m))}  [{source[m]}]")

    header = ["method", "locus", "rsid", "bp", "a1", "a2", "pip"]
    fh, pw, cw = [], {}, {}
    for method, _ in arms:
        lab = label(method, L)
        f = gzip.open(f"{OUT}/primary_{lab}_pip.csv.gz", "wt", newline="")
        g = gzip.open(f"{OUT}/primary_{lab}_cs_purity.csv.gz", "wt", newline="")
        fh += [f, g]
        pw[method], cw[method] = csv.writer(f), csv.writer(g)
        pw[method].writerow(header)
        cw[method].writerow(header + ["cs_ids_all", "n_cs_member", "native_purity",
                                      "native_purity_source", "common_purity"])
    MLAB = f"MACHINE{SUFFIX}" if NOAMBIG else "MACHINE"
    if L is None and (not NOAMBIG or MACHINE_ALT):
        f = gzip.open(f"{OUT}/primary_{MLAB}_pip.csv.gz", "wt", newline="")
        fh.append(f)
        mw = csv.writer(f)
        mw.writerow(header + ["CL_EUR", "CL_EAS"])

    status, n_rows, n_loci, miss, empty = [], Counter(), Counter(), Counter(), []
    cp_miss, cp_total = Counter(), Counter()
    for i, l in enumerate(loci, 1):
        D.L_OVERRIDE = None     # allele backfill: the primary MultiSuSiE table, same union input
        ms = f"{out_dir('MultiSuSiE')}/{l}_multisusie_pip.tsv"
        if NOAMBIG and not os.path.exists(ms):
            # no MultiSuSiE in this run: take the same fields from the union npz it copies
            z = np.load(f"{UNION}/{l}_harmonized.npz", allow_pickle=True)
            ref = {str(r): (str(p), str(e), str(o)) for r, p, e, o in
                   zip(z["rsids"], z["pos"], z["effect_allele"], z["other_allele"])}
        else:
            ref = {r["rsid"]: (r["pos"], r["effect_allele"], r["other_allele"])
                   for r in read_tsv(ms)}
        resolve, purity = common_purity_fn(l)
        for method, fn in arms:
            use_dir(method, L)
            rows, cs, native, src, st = fn(l, ref)
            lab = label(method, L)
            if not rows:
                empty.append((lab, l))
            info = {r[0]: r for r in rows}
            member, common = {}, {}
            for k, ms in cs.items():
                for m in ms:
                    member.setdefault(m, []).append(k)
                common[k] = purity([resolve(m, *info.get(m, (m, "", "", ""))[1:4]) for m in ms])
                cp_total[lab] += 1
                cp_miss[lab] += common[k] is None
            for rsid, bp, a1, a2, pip in rows:
                if not (bp and a1 and a2):
                    miss[lab] += 1
                pw[method].writerow([lab, l, rsid, bp, a1, a2, pip])
                ids = member.get(rsid, [])
                nat = [float(native[k]) for k in ids if native.get(k) not in (None, "")]
                com = [common[k] for k in ids if common[k] is not None]
                cw[method].writerow([lab, l, rsid, bp, a1, a2, pip, ";".join(ids), len(ids),
                                     min(nat) if nat else "", src if ids else "",
                                     min(com) if com else ""])
            n_rows[lab] += len(rows)
            n_loci[lab] += bool(rows)
            status.append([lab, l, st.get("n_cs", ""), st.get("converged", ""),
                           st.get("lambda", ""), len(rows), st.get("n_cs_ungated", "")]
                          + ([source[method]] if L is not None else []))
        if L is not None or (NOAMBIG and not MACHINE_ALT):
            if i % 25 == 0:
                print(f"  {i}/{len(loci)}", flush=True)
            continue

        # MACHINE: the full M x K fit@CL matrix, no pooled pip
        mrows = read_tsv(f"{out_dir('MACHINE')}/{l}_machine_b29_cl_bNULL.tsv")
        if not mrows:
            empty.append((MLAB, l))
        for r in mrows:
            bp, a1, a2 = ref.get(r["rsid"], ("", "", ""))
            if not (bp and a1 and a2):
                miss[MLAB] += 1
            mw.writerow([MLAB, l, r["rsid"], bp, a1, a2, "", r["CL_EUR"], r["CL_EAS"]])
        n_rows[MLAB] += len(mrows)
        n_loci[MLAB] += bool(mrows)
        # n_cs left blank: the diag file's n_cs is length(fit@CS), which is 2 at all 250 loci --
        # the slot's shape, not a count of credible sets.
        status.append([MLAB, l, "", "", "", len(mrows), ""])
        if i % 25 == 0:
            print(f"  {i}/{len(loci)}", flush=True)

    for f in fh:
        f.close()
    name = ("primary_status.csv.gz" if L is None else
            f"primary_status_L{L}.csv.gz" if not a_.arms else
            f"primary_status_L{L}_{'+'.join(sorted(DRIVER_ARM.get(m, m) for m, _ in arms))}.csv.gz")
    with gzip.open(f"{OUT}/{name}", "wt", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "locus", "n_cs", "converged", "lambda", "n_variants", "n_cs_ungated"]
                   + (["source"] if L is not None else []))
        w.writerows(status)

    print("rows per arm:", dict(n_rows))
    print("loci per arm:", dict(n_loci))
    print("rows missing bp/a1/a2:", dict(miss) or "none")
    print("empty (arm, locus):", empty or "none")
    print("CS with no common_purity / total:",
          {m: f"{cp_miss[m]}/{cp_total[m]}" for m in cp_total})


if __name__ == "__main__":
    main()
