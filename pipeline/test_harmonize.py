#!/usr/bin/env python3
"""End-to-end test of harmonize.py on one synthetic locus (no bucket data).

Five master variants:
  v1 1:200:A:G  EUR ss+LD, EAS ss+LD, core   EAS sumstats and EAS LD both stored swapped
  v2 1:300:C:T  EUR ss+LD, EAS LD only       (EAS z missing -> NaN / censored)
  v3 1:400:A:T  palindromic, everywhere
  v4 1:500:G:T  EAS ss+LD only               (no EUR LD -> not in MultiSuSiE, merged or EUR)
  v5 1:600:A:C  EUR ss+LD, EAS ss+LD, not core
Run: python test_harmonize.py
"""
import gzip
import json
import os
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import harmonize as H  # noqa: E402
import ld_io  # noqa: E402

L = "chr1_100_1000"
MASTER = pd.DataFrame(dict(
    locus=L, vkey=["1:200:A:G", "1:300:C:T", "1:400:A:T", "1:500:G:T", "1:600:A:C"],
    rsid=["rs1", "rs2", "rs3", "rs4", "rs5"], chr="1", bp=["200", "300", "400", "500", "600"],
    ref_a1=["A", "C", "A", "G", "A"], ref_a2=["G", "T", "T", "T", "C"],
    palindromic=[False, False, True, False, False],
    in_eur_ss=[True, True, True, False, True], in_eas_ss=[True, False, True, True, True],
    in_core_ss=[True, True, True, False, False]))

# Signed ("true", ref_a1-coded) LD per ancestry, over its own stored variants.
EUR_V = ["1:200:A:G", "1:300:C:T", "1:400:A:T", "1:600:A:C"]
EUR_R = np.array([[1, .6, .2, .1], [.6, 1, .3, .0], [.2, .3, 1, .4], [.1, 0, .4, 1]], dtype=np.float32)
EAS_V = ["1:200:A:G", "1:300:C:T", "1:400:A:T", "1:500:G:T", "1:600:A:C"]
EAS_R = np.array([[1, .5, .1, .2, .3], [.5, 1, .2, .0, .1], [.1, .2, 1, .1, .2],
                  [.2, 0, .1, 1, .0], [.3, .1, .2, 0, 1]], dtype=np.float32)
EAS_SIGN = np.array([-1, 1, 1, 1, 1])  # rs1 stored with alleles swapped in the EAS panel


def _write(tmp):
    mp = os.path.join(tmp, "master.tsv.gz")
    MASTER.to_csv(mp, sep="\t", index=False)
    os.environ["MASTER_TABLE"] = mp
    ld_io._master_cache.clear()
    m = MASTER.set_index("vkey")
    for src, V, R, sg in (("eur_gctb", EUR_V, EUR_R, np.ones(len(EUR_V), int)), ("eas_aou", EAS_V, EAS_R, EAS_SIGN)):
        d = os.path.join(tmp, "ld", f"ld_per_anc_{src}", L)
        os.makedirs(d)
        np.savez(os.path.join(d, "R.npz"), R=(R * np.outer(sg, sg)).astype(np.float32))  # stored in source coding
        a1 = [m.at[v, "ref_a1"] if s == 1 else m.at[v, "ref_a2"] for v, s in zip(V, sg)]
        a2 = [m.at[v, "ref_a2"] if s == 1 else m.at[v, "ref_a1"] for v, s in zip(V, sg)]
        pd.DataFrame(dict(vkey=V, src_a1=a1, src_a2=a2, sign=sg)).to_csv(os.path.join(d, "snps.tsv.gz"), sep="\t", index=False)
    rows = {"eur": [("200", "A", "G", 2.0, .10, 50000), ("300", "C", "T", 1.0, .20, 50000),
                    ("400", "A", "T", -1.0, .30, 50000), ("600", "A", "C", 0.5, .40, 60000)],
            "eas": [("200", "G", "A", 3.0, .70, 14000), ("400", "A", "T", 1.5, .30, 15000),
                    ("500", "G", "T", 4.0, .50, 14000), ("600", "A", "C", 0.1, .40, 16000)],
            "core": [("200", "A", "G", 2.5, .20, 70000), ("300", "T", "C", 1.0, .80, 70000),
                     ("400", "A", "T", -0.5, .30, 70000)]}
    for anc, rr in rows.items():
        d = os.path.join(tmp, "ss", f"{anc}_sumstats")
        os.makedirs(d)
        pd.DataFrame(rr, columns=["POS", "Allele1", "Allele2", "Zscore", "Freq1", "Weight"]).to_csv(
            os.path.join(d, f"{L}_sumstats.tsv.gz"), sep="\t", index=False)
    return dict(ld_root=os.path.join(tmp, "ld"), ss_root=os.path.join(tmp, "ss"))


def run(drop):
    tmp = tempfile.mkdtemp()
    roots = _write(tmp)
    out = os.path.join(tmp, "out")
    log = H.harmonize(L, out, H.ARMS, drop, "eur_gctb", roots["ld_root"], roots["ss_root"])
    return out, log


def test_domains_and_signs():
    out, log = run(False)
    assert log["domains"]["susiex_EUR"] == 4 and log["domains"]["susiex_EAS"] == 4
    assert log["domains"]["multisusie"] == 4          # rs1 rs2 rs3 rs5; rs4 lacks EUR LD
    assert log["domains"]["mesusie"] == 3             # rs1 rs3 rs5
    assert log["domains"]["merged"] == 3              # LD both & core: rs1 rs2 rs3
    assert log["domains"]["eur"] == 4
    # EAS sumstats for rs1 are G/A (swapped vs ref A/G): z 3.0 -> -3.0, freq .7 -> .3
    sst = pd.read_csv(f"{out}/union/{L}_susiex/EAS.sst.txt", sep="\t")
    r1 = sst.set_index("snp").loc["rs1"]
    assert (r1["A1"], r1["A2"], r1["beta"]) == ("A", "G", -3.0)
    frq = pd.read_csv(f"{out}/union/{L}_susiex/EAS_frq.frq", sep=" ").set_index("SNP")
    assert np.isclose(frq.at["rs1", "MAF"], 0.3) and frq.at["rs1", "NCHROBS"] == 2 * 14500
    # EAS LD was stored with rs1 swapped; after ld_io, R(rs1, rs2) must be the true +0.5
    R = np.fromfile(f"{out}/union/{L}_susiex/EAS.ld.bin", dtype=np.float32).reshape(4, 4)
    assert list(sst["snp"]) == ["rs1", "rs3", "rs4", "rs5"] and np.isclose(R[0, 1], 0.1) and np.isclose(R[0, 3], 0.3)


def test_npz_missing_and_censored():
    out, log = run(False)
    d = np.load(f"{out}/union/{L}_harmonized.npz", allow_pickle=True)
    c = np.load(f"{out}/union_censored_bothld/{L}_harmonized.npz", allow_pickle=True)
    assert list(d["rsids"]) == ["rs1", "rs2", "rs3", "rs5"] and list(d["ancestries"]) == ["EUR", "EAS"]
    assert np.isnan(d["z"][1, 1]) and d["z"][1, 0] == -3.0 and d["z"][0, 1] == 1.0
    assert np.isclose(d["R"][1][0, 1], 0.5)            # uncensored: real EAS LD kept for rs2
    assert list(c["rsids"]) == list(d["rsids"])        # both-LD sensitivity: same domain
    assert c["z"][1, 1] == 0.0 and np.all(c["R"][1][1, :] == 0) and np.all(c["R"][1][:, 1] == 0)
    assert np.isclose(c["R"][1][0, 2], 0.1) and np.array_equal(c["R"][0], d["R"][0])  # EUR untouched
    assert np.isclose(d["n"][0], 50000) and np.isclose(d["n"][1], 14500)


def test_censored_full_union():
    out, log = run(False)
    u = np.load(f"{out}/union_censored/{L}_harmonized.npz", allow_pickle=True)
    assert list(u["rsids"]) == ["rs1", "rs2", "rs3", "rs4", "rs5"]   # rs4: EAS-only, no EUR LD
    assert u["z"][0, 3] == 0.0 and np.all(u["R"][0][3, :] == 0) and np.all(u["R"][0][:, 3] == 0)
    assert u["z"][1, 3] == 4.0 and np.isclose(u["R"][1][0, 3], 0.2) and u["R"][1][3, 3] == 1.0
    assert u["z"][1, 1] == 0.0 and np.all(u["R"][1][1, :] == 0)       # rs2: no EAS z
    assert np.isclose(u["R"][0][0, 1], 0.6) and np.isclose(u["R"][1][0, 4], 0.3)
    assert log["domains"]["union_censored"] == 5 and log["domains"]["union_censored_bothld"] == 4
    assert log["domains"]["union_censored_missing"] == {"EUR": 1, "EAS": 1}


def test_merged_weights_and_formula():
    out, log = run(False)
    w = log["merge_weights"]                            # over LD-both = rs1 rs2 rs3 rs5
    assert w["EUR"] == 50000.0 and w["EAS"] == 15000.0  # median of (50k,50k,50k,60k); (14k,15k,16k)
    z = pd.read_csv(f"{out}/merged/{L}.z", sep=" ")
    assert list(z["rsid"]) == ["rs1", "rs2", "rs3"]
    assert z.set_index("rsid").at["rs2", "beta"] == -1.0  # core T/C vs ref C/T
    R = np.loadtxt(f"{out}/merged/{L}.ld")
    expect = (50000 * 0.6 + 15000 * 0.5) / 65000
    assert abs(R[0, 1] - expect) < 1e-6 and np.allclose(np.diag(R), 1)


def test_drop_ambiguous():
    out, log = run(True)
    assert log["domains"]["mesusie"] == 2 and log["domains"]["merged"] == 2 and log["domains"]["eur"] == 3
    assert log["merge_weights"] == {"EUR": 50000.0, "EAS": 15000.0}  # weights taken before the filter
    ma = pd.read_csv(f"{out}/intersect/{L}.EUR.ma", sep="\t")
    assert "rs3" not in set(ma["SNP"])
    with gzip.open(f"{out}/intersect/{L}.EAS.ld.gz", "rt") as f:
        assert np.loadtxt(f).shape == (2, 2)


def test_eur_arm_files():
    out, _ = run(False)
    z = pd.read_csv(f"{out}/eur/{L}.z", sep=" ")
    R = np.loadtxt(f"{out}/eur/{L}.ld")
    snps = pd.read_csv(f"{out}/eur/{L}_snps.txt.gz", sep="\t")
    assert len(z) == R.shape[0] == len(snps) == 4 and list(snps["SNP_ID"]) == list(z["rsid"])
    master = open(f"{out}/eur/{L}.master").read().splitlines()[1].split(";")
    assert master[0].endswith(f"{L}.z") and int(master[-1]) == H.P.N_GWAS_PUBLISHED_EUR
    with open(f"{out}/_logs/{L}.json") as f:
        assert json.load(f)["n_valid"] == {"EUR": 4, "EAS": 4}


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in tests:
        fn()
        print("PASS", fn.__name__)
    print(f"{len(tests)} passed")
