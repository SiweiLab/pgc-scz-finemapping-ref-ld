#!/usr/bin/env python3
"""Readers for each method's raw output, shared by export_results.py (which sets OUT and the
label globals before calling them): PIPs, credible sets and per-locus status, in the schema the
figure code reads.

  primary_<arm>_pip.csv.gz        method, locus, rsid, bp, a1, a2, pip
  primary_<arm>_cs_purity.csv.gz  ... pip, cs_ids_all, n_cs_member, native_purity,
                                  native_purity_source, common_purity
  primary_status.csv.gz           method, locus, n_cs, converged, lambda, n_variants,
                                  n_cs_ungated (SuSiEx only)

SuSiEx's L sweep ran mult-step off, a different setting from the primary SuSiEx-ours arm, so it
is labelled SuSiEx-nomultistep and never reads the primary output.
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
import run_finemap as D  # noqa: E402
from run_finemap import INPUTS, RSCRIPT, UNION, WORK, out_dir  # noqa: E402
from export_cs_purity import max_of_mins  # noqa: E402

MIN_PURITY_F = float(MIN_PURITY)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = None  # set by export_results.py
NOAMBIG = False
SUFFIX = "-noambig"
MACHINE_ALT = False
TAG_RE = re.compile(r"^chr[\dXYM]+:(\d+):([ACGT]+):([ACGT]+)$", re.IGNORECASE)

EUR_LD = f"{INPUTS}/ld_per_anc_eur"
EAS_LD = f"{INPUTS}/ld_per_anc_eas"

SUSIEX = "SuSiEx-ours"   # method label used for our SuSiEx runs


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
        # loci every requested arm has an output for (two loci were rerun for some arms)
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


DRIVER_ARM = {SUSIEX: "SuSiEx"}          # export label -> run_finemap arm name


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
