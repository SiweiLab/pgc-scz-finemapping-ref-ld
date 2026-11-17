#!/usr/bin/env python3
"""Per-variant credible-set membership + purity for MultiSuSiE,
MESuSiE, SuSiE-merged and SuSiEx-ours, on both each tool's own ("native") purity
composition and one shared ("common") composition, so the two are separable.

Background:
MultiSuSiE and SuSiEx already have purity on disk; MESuSiE and SuSiE-merged compute it
internally and discard it before writing, but it is recoverable from on-disk CS
membership + the frozen LD panel each arm actually used -- no refit. The tools do not
compose cross-ancestry purity the same way: SuSiEx is max-of-mins
(max_a(min_p |r|_a)), MultiSuSiE/MESuSiE are min-of-maxes (min_p(max_a |r|_a)); by the
minimax inequality SuSiEx's number is always <= the others' on the same panels, so a
shared min_purity=0.5 is a harsher gate in SuSiEx. `common_purity` fixes a single
composition (max-of-mins, matching SuSiEx's own) over one shared pair of LD panels so
all four arms are comparable on identical terms; `native_purity` keeps each arm's own
number so its own gate is still visible.

Native LD source per arm (whatever that arm's own fit actually saw):
  MultiSuSiE -- min_abs_corr already written to _cs.tsv (on_disk, no computation).
  SuSiEx-ours -- CS_PURITY already written to _susiex.summary (on_disk).
  MESuSiE    -- harmonized_v1/{locus}.{EUR,EAS}.ma+.ld.gz (recomputed, min-of-maxes).
  SuSiE-merged -- finemap_in/{locus}.z+.ld, the single merged panel (recomputed,
                  single-panel min -- no ancestry axis to compose).

Common LD source, identical for every arm: per_anc_ld_EUR/ + (per_anc_ld_qcexcl_
derived_239loci|known11)/EAS/ -- the b6b7 backup's own per-ancestry reference panels,
independent of which harmonized/derived tree any given arm's fit actually used.
max-of-mins over these two, always. For SuSiE-merged this is a genuine counterfactual
(that arm never saw per-ancestry LD at all) -- reported as such, not as its own gate.

A CS's common_purity requires EVERY member to resolve in an ancestry's panel before
that ancestry contributes to the max -- a partially-resolved ancestry is dropped
entirely rather than computed on a subset (per the brief: no partial-set substitution).
If neither ancestry can fully resolve a CS, common_purity is left blank for that CS's
members and counted as a miss, never silently dropped or set to 0.

Per-variant PIP is copied byte-for-byte from the same source rows/logic
`export_full_pip_vectors.py` uses, so it stays byte-identical to
`results/benchmark/pip_vectors_full250_qcexcl.csv.gz` -- this export must not move a PIP.
"""
import csv
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_params import MIN_PURITY  # noqa: E402  -- the only place this literal may live

MIN_PURITY_F = float(MIN_PURITY)

BACKUP = os.path.expanduser("~/work/full250_qcexcl_backup")
B6B7 = os.path.expanduser("~/workspace/dataproc-staging-post-gwas/b6b7_backup_2026-09-09")
HARM_V1 = os.path.expanduser("~/work/full250_qcexcl_derived/harmonized_v1")
FINEMAP_IN = f"{B6B7}/finemap_in"
EUR_LD_DIR = f"{B6B7}/per_anc_ld_EUR"
EAS_LD_DIRS = [f"{B6B7}/per_anc_ld_qcexcl_derived_239loci/EAS",
               f"{B6B7}/per_anc_ld_qcexcl_known11/EAS"]
TASK_A = os.path.expanduser("~/work/benchmark")
OUT_DIR = os.path.expanduser("~/work/cs_purity")
os.makedirs(OUT_DIR, exist_ok=True)

TAG_RE = re.compile(r"^(chr[\dXYM]+):(\d+):([ACGT]+):([ACGT]+)$", re.IGNORECASE)


def read_tsv(path, delim="\t"):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter=delim))


def parse_tag(tag):
    m = TAG_RE.match(tag)
    if m:
        return int(m.group(2)), m.group(3), m.group(4)
    return None, None, None


def loci_list():
    names = os.listdir(f"{BACKUP}/multisusie_out")
    return sorted({n.rsplit("_multisusie_pip.tsv", 1)[0] for n in names
                   if n.endswith("_multisusie_pip.tsv")})


# ---------------------------------------------------------------------------
# rsid/bp/a1/a2 crosswalk -- identical logic to export_full_pip_vectors.py's
# build_ref_table, extended with a bp->rsid reverse index for the SuSiE-merged
# common_purity join (that arm's own identifier is a chr:pos:a1:a2 tag, not a rsid).
# ---------------------------------------------------------------------------
def build_ref_table(locus):
    ref = {}
    for row in read_tsv(f"{BACKUP}/multisusie_out/{locus}_multisusie_pip.tsv"):
        ref[row["rsid"]] = (row["pos"], row["effect_allele"], row["other_allele"])
    fm_path = f"{BACKUP}/finemap_in_snp_only/{locus}.snp"
    if os.path.exists(fm_path):
        with open(fm_path) as f:
            for row in csv.DictReader(f, delimiter=" ", skipinitialspace=True):
                if row["rsid"] not in ref:
                    ref[row["rsid"]] = (row["position"], row["allele1"], row["allele2"])
    bp_index = {}
    for rsid, (bp, a1, a2) in ref.items():
        bp_index.setdefault(str(bp), []).append((rsid, a1.upper(), a2.upper()))
    return ref, bp_index


def resolve_tag_to_rsid(bp, a1, a2, bp_index):
    """SuSiE-merged's chr:pos:a1:a2 tag -> rsid, via the bp crosswalk, disambiguated
    by allele SET (order-agnostic: a strand/effect flip can swap a1/a2). The
    crosswalk's finemap_in_snp_only fallback carries tag-shaped "rsids" for variants
    with no true rsID (see export_full_pip_vectors.py's own docstring) -- those are
    useless as a target (the per-ancestry panels are rsID-keyed), so they are excluded
    from candidates rather than left to cause a spurious ambiguity against a real rsID
    at the same position/alleles."""
    cands = bp_index.get(str(bp), [])
    a_set = {a1.upper(), a2.upper()}
    matches = [rsid for rsid, ca1, ca2 in cands
               if {ca1, ca2} == a_set and not TAG_RE.match(rsid)]
    if len(matches) == 1:
        return matches[0]
    return None  # 0 or ambiguous >1 -- a miss, not a guess


# ---------------------------------------------------------------------------
# LD panel loaders
# ---------------------------------------------------------------------------
def load_npz_panel(ld_dir, locus, id_col):
    r_path = f"{ld_dir}/{locus}_R.npz"
    s_path = f"{ld_dir}/{locus}_snps.txt.gz"
    if not (os.path.exists(r_path) and os.path.exists(s_path)):
        return None, None
    R = np.load(r_path)["R"].astype(np.float64)
    snps = pd.read_csv(s_path, sep="\t")
    idx = {v: i for i, v in enumerate(snps[id_col].values) if isinstance(v, str) and v}
    return R, idx


def load_common_panels(locus):
    R_eur, eur_idx = load_npz_panel(EUR_LD_DIR, locus, "MarkerName")
    R_eas, eas_idx = None, None
    for d in EAS_LD_DIRS:
        R_eas, eas_idx = load_npz_panel(d, locus, "SNP_ID")
        if R_eas is not None:
            break
    return R_eur, eur_idx, R_eas, eas_idx


def load_harmonized_v1_panel(locus):
    """MESuSiE's own native LD -- EUR+EAS share one SNP order (harmonize_loci.py's
    guarantee, verified in task1_inventory_report.md)."""
    eur_ma = f"{HARM_V1}/{locus}.EUR.ma"
    eas_ld = f"{HARM_V1}/{locus}.EAS.ld.gz"
    eur_ld = f"{HARM_V1}/{locus}.EUR.ld.gz"
    if not (os.path.exists(eur_ma) and os.path.exists(eur_ld) and os.path.exists(eas_ld)):
        return None, None, None
    snp_order = pd.read_csv(eur_ma, sep="\t")["SNP"].tolist()
    idx = {s: i for i, s in enumerate(snp_order)}
    R_eur = pd.read_csv(eur_ld, sep="\t", header=None).values.astype(np.float64)
    R_eas = pd.read_csv(eas_ld, sep="\t", header=None).values.astype(np.float64)
    return idx, R_eur, R_eas


def load_finemap_merged_panel(locus):
    z_path = f"{FINEMAP_IN}/{locus}.z"
    ld_path = f"{FINEMAP_IN}/{locus}.ld"
    if not (os.path.exists(z_path) and os.path.exists(ld_path)):
        return None, None
    tags = pd.read_csv(z_path, sep=r"\s+")["rsid"].tolist()
    idx = {t: i for i, t in enumerate(tags)}
    R = pd.read_csv(ld_path, sep=r"\s+", header=None).values.astype(np.float64)
    return idx, R


# ---------------------------------------------------------------------------
# Purity formulas
# ---------------------------------------------------------------------------
def _upper_tri_min(R, idxs):
    if len(idxs) <= 1:
        return 1.0
    sub = np.abs(R[np.ix_(idxs, idxs)])
    iu = np.triu_indices(len(idxs), k=1)
    vals = sub[iu]
    return float(np.min(vals)) if vals.size else 1.0


def resolve_indices(members, idx_map):
    idxs = []
    for m in members:
        if m not in idx_map:
            return None
        idxs.append(idx_map[m])
    return idxs


def max_of_mins(members, eur_idx, eas_idx, R_eur, R_eas):
    """The common_purity composition, and SuSiEx's own native composition."""
    if len(members) <= 1:
        return 1.0, True, True
    eur_idxs = resolve_indices(members, eur_idx) if eur_idx is not None else None
    eas_idxs = resolve_indices(members, eas_idx) if eas_idx is not None else None
    vals = []
    if eur_idxs is not None:
        vals.append(_upper_tri_min(R_eur, eur_idxs))
    if eas_idxs is not None:
        vals.append(_upper_tri_min(R_eas, eas_idxs))
    if not vals:
        return None, eur_idxs is not None, eas_idxs is not None
    return max(vals), eur_idxs is not None, eas_idxs is not None


def min_of_maxes(members, shared_idx, R_eur, R_eas):
    """MultiSuSiE/MESuSiE's own native composition -- both ancestries share one
    SNP-order index by construction (harmonize_loci.py), so this is all-or-nothing
    per CS, not per pair."""
    if len(members) <= 1:
        return 1.0
    idxs = resolve_indices(members, shared_idx)
    if idxs is None:
        return None
    sub_eur = np.abs(R_eur[np.ix_(idxs, idxs)])
    sub_eas = np.abs(R_eas[np.ix_(idxs, idxs)])
    pmax = np.maximum(sub_eur, sub_eas)
    iu = np.triu_indices(len(idxs), k=1)
    vals = pmax[iu]
    return float(np.min(vals)) if vals.size else 1.0


def merged_only_min(members, idx_map, R):
    idxs = resolve_indices(members, idx_map)
    if idxs is None:
        return None
    return _upper_tri_min(R, idxs)


# ---------------------------------------------------------------------------
# Per-arm row builders. Each returns:
#   pip_rows: {rsid: (bp, a1, a2, pip_string)}   (pip_string copied verbatim / same
#             computation as export_full_pip_vectors.py -- never re-derived from a
#             different formula, to guarantee byte-identical output)
#   cs_members: {cs_id: [rsid, ...]}
#   native_purity_by_cs: {cs_id: float or None}
#   native_source: "on_disk" | "recomputed"
# ---------------------------------------------------------------------------
def multisusie_arm(locus, ref, bp_index):
    pip_rows = {}
    for row in read_tsv(f"{BACKUP}/multisusie_out/{locus}_multisusie_pip.tsv"):
        pip_rows[row["rsid"]] = (row["pos"], row["effect_allele"], row["other_allele"], row["pip"])
    cs_members, native = {}, {}
    for row in read_tsv(f"{BACKUP}/multisusie_out/{locus}_multisusie_cs.tsv"):
        cs_members.setdefault(row["cs_id"], []).append(row["rsid"])
        native[row["cs_id"]] = float(row["min_abs_corr"])
    return pip_rows, cs_members, native, "on_disk"


def mesusie_arm(locus, ref, bp_index):
    pip_rows = {}
    n_backfill_miss = 0
    for row in read_tsv(f"{BACKUP}/mesusie_out/{locus}_mesusie_pip.tsv"):
        bp, a1, a2 = ref.get(row["rsid"], ("", "", ""))
        if bp == "":
            n_backfill_miss += 1
        pip_rows[row["rsid"]] = (bp, a1, a2, row["pip"])
    cs_members = {}
    for row in read_tsv(f"{BACKUP}/mesusie_out/{locus}_mesusie_config.tsv"):
        cs_members.setdefault(row["cs_id"], []).append(row["rsid"])
    shared_idx, R_eur, R_eas = load_harmonized_v1_panel(locus)
    native = {}
    for cs_id, members in cs_members.items():
        if shared_idx is None:
            native[cs_id] = None
            continue
        native[cs_id] = min_of_maxes(members, shared_idx, R_eur, R_eas)
    return pip_rows, cs_members, native, "recomputed", n_backfill_miss


def susie_merged_arm(locus):
    pip_rows = {}
    cs_members = {}
    for row in read_tsv(f"{BACKUP}/susie_merged_out/{locus}_susie_pip.tsv"):
        bp, a1, a2 = parse_tag(row["rsid"])
        pip_rows[row["rsid"]] = (bp if bp is not None else row.get("position", ""),
                                  a1 or "", a2 or "", row["PIP"])
        if row["in_cs"] != ".":
            cs_members.setdefault(row["in_cs"], []).append(row["rsid"])
    idx_map, R = load_finemap_merged_panel(locus)
    native = {}
    for cs_id, members in cs_members.items():
        native[cs_id] = merged_only_min(members, idx_map, R) if idx_map is not None else None
    return pip_rows, cs_members, native, "recomputed"


def susiex_arm(locus, ref, bp_index):
    pip_rows = {}
    n_backfill_miss = 0
    snp_f = read_tsv(f"{BACKUP}/susiex_out/{locus}_susiex.snp")
    for row in snp_f:
        cs_cols = [k for k in row if k.startswith("PIP(CS")]
        pip = max((float(row[c]) for c in cs_cols), default=0.0)
        bp, a1, a2 = ref.get(row["SNP"], (row.get("BP", ""), "", ""))
        if a1 == "":
            n_backfill_miss += 1
        pip_rows[row["SNP"]] = (bp, a1, a2, pip)
    cs_members = {}
    for row in read_tsv(f"{BACKUP}/susiex_out/{locus}_susiex.cs"):
        cs_members.setdefault(row["CS_ID"], []).append(row["SNP"])
    native = {}
    summ_path = f"{BACKUP}/susiex_out/{locus}_susiex.summary"
    if os.path.exists(summ_path):
        with open(summ_path) as f:
            lines = f.readlines()
        if len(lines) > 1 and lines[1].strip() != "NULL":
            for row in csv.DictReader(lines[1:], delimiter="\t"):
                native[row["CS_ID"]] = float(row["CS_PURITY"])
    return pip_rows, cs_members, native, "on_disk", n_backfill_miss


# ---------------------------------------------------------------------------
# Assembly: turn (pip_rows, cs_members, native) into per-variant output rows
# ---------------------------------------------------------------------------
def assemble(method, locus, pip_rows, cs_members, native_by_cs, native_source,
             R_eur_c, eur_idx_c, R_eas_c, eas_idx_c, merged_tag_lookup=None,
             bp_index=None):
    """merged_tag_lookup: only for SuSiE-merged -- maps its tag rsid -> common-panel
    rsid via the bp/allele crosswalk, for the common_purity computation only."""
    variant_cs = {}
    for cs_id, members in cs_members.items():
        for m in members:
            variant_cs.setdefault(m, []).append(cs_id)

    common_by_cs = {}
    ancestry_hits = {"eur_only": 0, "eas_only": 0, "both": 0, "neither": 0}
    miss_members = 0
    total_members = 0
    for cs_id, members in cs_members.items():
        if merged_tag_lookup is not None:
            resolved = [merged_tag_lookup.get(m) for m in members]
            total_members += len(members)
            miss_members += sum(1 for r in resolved if r is None)
            if any(r is None for r in resolved):
                common_by_cs[cs_id] = (None, False, False)
                continue
            members_for_common = resolved
        else:
            members_for_common = members
            total_members += len(members)
        val, has_eur, has_eas = max_of_mins(members_for_common, eur_idx_c, eas_idx_c,
                                             R_eur_c, R_eas_c)
        common_by_cs[cs_id] = (val, has_eur, has_eas)
        if has_eur and has_eas:
            ancestry_hits["both"] += 1
        elif has_eur:
            ancestry_hits["eur_only"] += 1
        elif has_eas:
            ancestry_hits["eas_only"] += 1
        else:
            ancestry_hits["neither"] += 1
            miss_members += len(members) if merged_tag_lookup is None else 0

    rows = []
    for rsid, (bp, a1, a2, pip) in pip_rows.items():
        cs_ids = variant_cs.get(rsid, [])
        native_vals = [native_by_cs[c] for c in cs_ids if native_by_cs.get(c) is not None]
        common_vals = [common_by_cs[c][0] for c in cs_ids if common_by_cs.get(c, (None,))[0] is not None]
        rows.append({
            "method": method, "locus": locus, "rsid": rsid, "bp": bp, "a1": a1, "a2": a2,
            "pip": pip,
            "cs_ids_all": ";".join(cs_ids),
            "n_cs_member": len(cs_ids),
            "native_purity": min(native_vals) if native_vals else "",
            "native_purity_source": native_source if cs_ids else "",
            "common_purity": min(common_vals) if common_vals else "",
        })
    stats = {"n_cs": len(cs_members), "ancestry_hits": ancestry_hits,
              "miss_members": miss_members, "total_members": total_members}
    return rows, stats


def write_gz(path, rows, header):
    with gzip.open(path, "wt", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)


HEADER = ["method", "locus", "rsid", "bp", "a1", "a2", "pip", "cs_ids_all",
          "n_cs_member", "native_purity", "native_purity_source", "common_purity"]

def main():
    loci = loci_list()
    print(f"{len(loci)} loci", flush=True)

    out = {"MultiSuSiE": [], "MESuSiE": [], "SuSiE-merged": [], "SuSiEx-ours": []}
    backfill_miss = {"MESuSiE": 0, "SuSiEx-ours": 0}
    cs_stats = {k: {"n_cs": 0, "ancestry_hits": {"eur_only": 0, "eas_only": 0, "both": 0, "neither": 0},
                     "miss_members": 0, "total_members": 0} for k in out}

    for i, locus in enumerate(loci, 1):
        ref, bp_index = build_ref_table(locus)
        R_eur_c, eur_idx_c, R_eas_c, eas_idx_c = load_common_panels(locus)

        # MultiSuSiE
        pip_rows, cs_members, native, src = multisusie_arm(locus, ref, bp_index)
        rows, st = assemble("MultiSuSiE", locus, pip_rows, cs_members, native, src,
                             R_eur_c, eur_idx_c, R_eas_c, eas_idx_c)
        out["MultiSuSiE"].extend(rows)
        for k in st["ancestry_hits"]:
            cs_stats["MultiSuSiE"]["ancestry_hits"][k] += st["ancestry_hits"][k]
        cs_stats["MultiSuSiE"]["n_cs"] += st["n_cs"]
        cs_stats["MultiSuSiE"]["miss_members"] += st["miss_members"]
        cs_stats["MultiSuSiE"]["total_members"] += st["total_members"]

        # MESuSiE
        pip_rows, cs_members, native, src, n_miss = mesusie_arm(locus, ref, bp_index)
        backfill_miss["MESuSiE"] += n_miss
        rows, st = assemble("MESuSiE", locus, pip_rows, cs_members, native, src,
                             R_eur_c, eur_idx_c, R_eas_c, eas_idx_c)
        out["MESuSiE"].extend(rows)
        for k in st["ancestry_hits"]:
            cs_stats["MESuSiE"]["ancestry_hits"][k] += st["ancestry_hits"][k]
        cs_stats["MESuSiE"]["n_cs"] += st["n_cs"]
        cs_stats["MESuSiE"]["miss_members"] += st["miss_members"]
        cs_stats["MESuSiE"]["total_members"] += st["total_members"]

        # SuSiE-merged
        pip_rows, cs_members, native, src = susie_merged_arm(locus)
        tag_lookup = {}
        for cs_id, members in cs_members.items():
            for m in members:
                if m in tag_lookup:
                    continue
                bp, a1, a2 = parse_tag(m)
                tag_lookup[m] = resolve_tag_to_rsid(bp, a1, a2, bp_index) if bp is not None else None
        rows, st = assemble("SuSiE-merged", locus, pip_rows, cs_members, native, src,
                             R_eur_c, eur_idx_c, R_eas_c, eas_idx_c,
                             merged_tag_lookup=tag_lookup, bp_index=bp_index)
        out["SuSiE-merged"].extend(rows)
        for k in st["ancestry_hits"]:
            cs_stats["SuSiE-merged"]["ancestry_hits"][k] += st["ancestry_hits"][k]
        cs_stats["SuSiE-merged"]["n_cs"] += st["n_cs"]
        cs_stats["SuSiE-merged"]["miss_members"] += st["miss_members"]
        cs_stats["SuSiE-merged"]["total_members"] += st["total_members"]

        # SuSiEx-ours
        pip_rows, cs_members, native, src, n_miss = susiex_arm(locus, ref, bp_index)
        backfill_miss["SuSiEx-ours"] += n_miss
        rows, st = assemble("SuSiEx-ours", locus, pip_rows, cs_members, native, src,
                             R_eur_c, eur_idx_c, R_eas_c, eas_idx_c)
        out["SuSiEx-ours"].extend(rows)
        for k in st["ancestry_hits"]:
            cs_stats["SuSiEx-ours"]["ancestry_hits"][k] += st["ancestry_hits"][k]
        cs_stats["SuSiEx-ours"]["n_cs"] += st["n_cs"]
        cs_stats["SuSiEx-ours"]["miss_members"] += st["miss_members"]
        cs_stats["SuSiEx-ours"]["total_members"] += st["total_members"]

        if i % 25 == 0:
            print(f"[progress] {i}/{len(loci)}", flush=True)

    for method, fname in [("MultiSuSiE", "multisusie_cs_purity.csv.gz"),
                           ("MESuSiE", "mesusie_cs_purity.csv.gz"),
                           ("SuSiE-merged", "susie_merged_cs_purity.csv.gz"),
                           ("SuSiEx-ours", "susiex_cs_purity.csv.gz")]:
        write_gz(f"{OUT_DIR}/{fname}", out[method], HEADER)
        print(f"wrote {OUT_DIR}/{fname}: {len(out[method])} rows", flush=True)

    print("\n=== backfill miss counts ===")
    for k, v in backfill_miss.items():
        print(f"{k}: {v}")

    print("\n=== CS-level common_purity resolvability ===")
    for k, v in cs_stats.items():
        print(f"{k}: n_cs={v['n_cs']} ancestry_hits={v['ancestry_hits']} "
              f"miss_members={v['miss_members']}/{v['total_members']}")


if __name__ == "__main__":
    main()
