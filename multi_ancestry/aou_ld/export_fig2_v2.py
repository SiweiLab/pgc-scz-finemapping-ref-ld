#!/usr/bin/env python3
"""B64: Fig. 2a-c data on pipeline v2, in B32's schema so metrics/plot_fig2_locus.py reads it unchanged.

Per locus, into OUT (default task_b64_out/):
  fig2_<L>_assoc.csv.gz      rsid, bp, a1, a2, z_eur, p_eur, z_eas, p_eas, z_merged, p_merged
                             z signed to the master's ref_a1 (= a1) by v2_io; z_merged is the core (trans-ancestry) meta z.
  fig2_<L>_r2_to_lead.csv.gz rsid, bp, lead_rsid, r2_eur, r2_eas, r2_merged [, focus_rsid, r2f_eur, r2f_eas, r2f_merged]
                             r2_eur = GCTB-10k, r2_eas = AoU EAS, r2_merged = (alpha R_EUR + (1 - alpha) R_EAS)^2 with the
                             harmonize_v2 merge weight. NaN where the variant is absent from that panel (never 0).
The variant set is the master's NON-palindromic variants at the locus (the primary analysis, D1). lead = Yuan's top PIP.
Usage (cluster; V2_LD_ROOT/V2_SS_ROOT as B50):
  python export_fig2_v2.py --tree ~/work/v2_harmonized/noambig \
      chr12_110129434_111489434 chr18_52716708_53630258 chr6_96319894_96477294:rs59107957
  (an optional :rsid after a locus adds the r2-to-focus columns)
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v2_io  # noqa: E402

BENCH = os.path.join(os.path.dirname(HERE), "metrics", "pip_vectors_susiex_published_full.csv.gz")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("loci", nargs="+")
    ap.add_argument("--tree", required=True, help="harmonize_v2 tree, for the merge weight (_logs)")
    ap.add_argument("--out", default=os.path.join(HERE, "task_b64_out"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    bench = pd.read_csv(BENCH, dtype={"rsid": str})
    for item in a.loci:
        L, focus = (item.split(":", 1) + [None])[:2]
        m = v2_io.load_master(L)
        m = m[~m["palindromic"]]
        w = json.load(open(os.path.join(os.path.expanduser(a.tree), "_logs", f"{L}.json")))["merge_weights"]
        alpha = w["EUR"] / (w["EUR"] + w["EAS"])
        ss = {k: v2_io.load_ss(k, L)["z"] for k in ("eur", "eas", "core")}
        ld = {k: v2_io.load_ld(k, L) for k in ("eur_gctb", "eas_aou")}
        b = bench[bench.locus == L]
        lead = b.loc[b.pip.idxmax(), "rsid"]
        rs2vk = dict(zip(m["rsid"], m.index))

        def row_r(src, rsid):
            R, s = ld[src]
            v = rs2vk.get(rsid)
            if v is None or v not in s.index:
                return pd.Series(np.nan, index=m.index)
            return pd.Series(R[s.index.get_loc(v)], index=s.index).reindex(m.index)

        out = pd.DataFrame({"rsid": m["rsid"].to_numpy(), "bp": m["bp"].to_numpy(),
                            "a1": m["ref_a1"].to_numpy(), "a2": m["ref_a2"].to_numpy()})
        for k, col in (("eur", "eur"), ("eas", "eas"), ("core", "merged")):
            z = ss[k].reindex(m.index).to_numpy()
            out[f"z_{col}"], out[f"p_{col}"] = z, 2 * norm.sf(np.abs(z))
        out.to_csv(os.path.join(a.out, f"fig2_{L}_assoc.csv.gz"), index=False)

        r2 = pd.DataFrame({"rsid": m["rsid"].to_numpy(), "bp": m["bp"].to_numpy(), "lead_rsid": lead})
        for tag, rsid in (("", lead), ("f", focus)):
            if rsid is None:
                continue
            e, s_ = row_r("eur_gctb", rsid), row_r("eas_aou", rsid)
            r2[f"r2{tag}_eur"], r2[f"r2{tag}_eas"] = (e ** 2).to_numpy(), (s_ ** 2).to_numpy()
            r2[f"r2{tag}_merged"] = ((alpha * e + (1 - alpha) * s_) ** 2).to_numpy()
        if focus:
            r2["focus_rsid"] = focus
        r2.to_csv(os.path.join(a.out, f"fig2_{L}_r2_to_lead.csv.gz"), index=False)
        print(L, "lead", lead, "focus", focus, "variants", len(out), "alpha", round(alpha, 5), flush=True)


if __name__ == "__main__":
    main()
