#!/usr/bin/env python3
"""Gate 5 unit tests for v2_io.py, on synthetic files (no bucket data needed).

Covers gate 5's first two items: a swapped variant's z and R row/col both flip, and a variant
absent from a source is absent, not zero. The third (round trip against v1's signed EUR files)
needs v1 on the cluster and runs there. Run: python -m pytest test_v2_io.py  (or python test_v2_io.py)
"""
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2_io  # noqa: E402

L = "chr1_100_1000"
# vkey, rsid, ref_a1, ref_a2, palindromic
MASTER = [("1:200:A:G", "rs1", "A", "G", False),
          ("1:300:C:T", "rs2", "T", "C", False),
          ("1:400:A:T", "rs3", "A", "T", True),
          ("1:500:C:G", "rs4", "G", "C", True)]
R_SRC = np.array([[1.0, 0.5, 0.2],
                  [0.5, 1.0, -0.3],
                  [0.2, -0.3, 1.0]], dtype=np.float32)


def _setup(tmp, snps_rows, ss_rows, R=R_SRC, source="eur_gctb", anc="eur"):
    m = pd.DataFrame(MASTER, columns=["vkey", "rsid", "ref_a1", "ref_a2", "palindromic"])
    m["locus"], m["chr"] = L, "1"
    m["bp"] = m["vkey"].str.split(":").str[1]
    for c in ("in_eur_ss", "in_eas_ss", "in_core_ss"):
        m[c] = True
    mp = os.path.join(tmp, "master.tsv.gz")
    m.to_csv(mp, sep="\t", index=False)
    d = os.path.join(tmp, "ld", f"ld_per_anc_{source}", L)
    os.makedirs(d, exist_ok=True)
    np.savez(os.path.join(d, "R.npz"), R=R)
    pd.DataFrame(snps_rows).to_csv(os.path.join(d, "snps.tsv.gz"), sep="\t", index=False)
    sd = os.path.join(tmp, "ss", f"{anc}_sumstats")
    os.makedirs(sd, exist_ok=True)
    pd.DataFrame(ss_rows).to_csv(os.path.join(sd, f"{L}_sumstats.tsv.gz"), sep="\t", index=False)
    v2_io._master_cache.clear()
    return dict(ld_root=os.path.join(tmp, "ld"), master=mp), dict(ss_root=os.path.join(tmp, "ss"), master=mp)


# rs1 same orientation, rs2 swapped (src C/T vs ref T/C), rs3 same; rs4 is absent from the LD.
SNPS = [dict(vkey="1:200:A:G", src_a1="A", src_a2="G", sign=1),
        dict(vkey="1:300:C:T", src_a1="C", src_a2="T", sign=-1),
        dict(vkey="1:400:A:T", src_a1="A", src_a2="T", sign=1)]
SS = [dict(POS=200, Allele1="A", Allele2="G", Zscore=2.0, Effect=0.2, Freq1=0.1, Weight=1000),
      dict(POS=300, Allele1="C", Allele2="T", Zscore=3.0, Effect=0.3, Freq1=0.3, Weight=1000),
      dict(POS=400, Allele1="A", Allele2="T", Zscore=-1.0, Effect=-0.1, Freq1=0.4, Weight=1000),
      dict(POS=600, Allele1="A", Allele2="C", Zscore=9.0, Effect=0.9, Freq1=0.5, Weight=1000)]  # not in master


def test_swapped_variant_flips_R_row_and_col():
    with tempfile.TemporaryDirectory() as t:
        kl, _ = _setup(t, SNPS, SS)
        R, s = v2_io.load_ld("eur_gctb", L, **kl)
        expect = R_SRC * np.outer([1, -1, 1], [1, -1, 1])
        assert np.allclose(R, expect)
        assert R[0, 1] == -0.5 and R[1, 2] == 0.3 and R[0, 2] == np.float32(0.2)
        assert np.allclose(np.diag(R), 1.0)
        assert list(s.index) == ["1:200:A:G", "1:300:C:T", "1:400:A:T"]
        assert list(s["rsid"]) == ["rs1", "rs2", "rs3"]


def test_swapped_variant_flips_z_beta_freq():
    with tempfile.TemporaryDirectory() as t:
        _, ks = _setup(t, SNPS, SS)
        z = v2_io.load_ss("eur", L, **ks)
        assert z.loc["1:200:A:G", "z"] == 2.0 and z.loc["1:300:C:T", "z"] == -3.0
        assert z.loc["1:300:C:T", "beta"] == -0.3 and np.isclose(z.loc["1:300:C:T", "freq"], 0.7)
        assert z.loc["1:200:A:G", "freq"] == 0.1 and z.loc["1:300:C:T", "n"] == 1000
        assert z.attrs["dropped"] == {"outside_master": 1, "allele_mismatch": 0}


def test_absent_is_absent_not_zero():
    with tempfile.TemporaryDirectory() as t:
        kl, ks = _setup(t, SNPS, SS)
        R, s = v2_io.load_ld("eur_gctb", L, **kl)
        z = v2_io.load_ss("eur", L, **ks)
        assert "1:500:C:G" not in s.index and R.shape == (3, 3)
        assert "1:500:C:G" not in z.index and "1:600:A:C" not in z.index
        try:
            v2_io.subset(R, s, ["1:200:A:G", "1:500:C:G"])
        except v2_io.V2IOError:
            pass
        else:
            raise AssertionError("subset zero-filled an absent variant")
        Rs, ss = v2_io.subset(R, s, ["1:400:A:T", "1:200:A:G"])
        assert list(ss.index) == ["1:400:A:T", "1:200:A:G"] and Rs[0, 1] == np.float32(0.2)


def test_wrong_stored_sign_is_refused():
    bad = [dict(r) for r in SNPS]
    bad[1]["sign"] = 1
    with tempfile.TemporaryDirectory() as t:
        kl, _ = _setup(t, bad, SS)
        try:
            v2_io.load_ld("eur_gctb", L, **kl)
        except v2_io.V2IOError as e:
            assert "stored sign disagrees" in str(e)
        else:
            raise AssertionError("a wrong stored sign was accepted")


def test_minus_strand_rows_sign_through_the_complement():
    # AoU stores + strand alleles. rs1 on the minus strand: src T/C, complement A/G = ref -> +1;
    # rs2 minus: src G/A, complement C/T = ref swapped -> -1.
    rows = [dict(vkey="1:200:A:G", src_a1="T", src_a2="C", sign=1, strand="-"),
            dict(vkey="1:300:C:T", src_a1="G", src_a2="A", sign=-1, strand="-"),
            dict(vkey="1:400:A:T", src_a1="A", src_a2="T", sign=1, strand="+")]
    with tempfile.TemporaryDirectory() as t:
        kl, _ = _setup(t, rows, SS, source="eas_aou")
        R, s = v2_io.load_ld("eas_aou", L, **kl)
        assert list(s["sign"]) == [1, -1, 1] and R[0, 1] == -0.5


def test_allele_mismatch_is_dropped_not_signed():
    ss = [dict(r) for r in SS]
    ss[0]["Allele1"], ss[0]["Allele2"] = "A", "C"  # 1:200:A:C is not a master vkey -> outside
    ss[1]["Allele1"] = "c"  # lower case still matches
    with tempfile.TemporaryDirectory() as t:
        _, ks = _setup(t, SNPS, ss)
        z = v2_io.load_ss("eur", L, **ks)
        assert "1:200:A:G" not in z.index and z.loc["1:300:C:T", "z"] == -3.0


def test_duplicate_vkey_is_refused():
    dup = SNPS + [dict(SNPS[0])]
    R4 = np.eye(4, dtype=np.float32)
    with tempfile.TemporaryDirectory() as t:
        kl, _ = _setup(t, dup, SS, R=R4)
        try:
            v2_io.load_ld("eur_gctb", L, **kl)
        except v2_io.V2IOError as e:
            assert "duplicate vkey" in str(e)
        else:
            raise AssertionError("duplicate vkey accepted")


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in tests:
        fn()
        print("PASS", fn.__name__)
    print(f"{len(tests)} passed")
