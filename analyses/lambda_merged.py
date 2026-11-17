"""SuSiE-RSS lambda per locus on the merged LD: susieR::estimate_s_rss on the core z and the merged LD, as
recorded by the SuSiE-merged fit of each run. Writes results/lambda_merged/lambda_merged_<variant>.csv
(locus, lambda); the paper uses the noambig file (Supplementary Tables).

  python lambda_merged.py            # reads ~/work/primary_v2_{noambig,primary}/SuSiE-merged
"""
import os, re
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results", "lambda_merged")


def lambdas(work):
    rows = []
    d = os.path.expanduser(f"{work}/SuSiE-merged")
    for f in os.listdir(d):
        if f.endswith("_summary.tsv"):
            s = pd.read_csv(os.path.join(d, f), sep="\t")
            rows.append((s["locus"].iloc[0], float(s["lambda"].iloc[0])))
    if not rows:  # fall back to logs
        for f in os.listdir(f"{d}/_logs"):
            m = re.search(r"lambda: ([\d.eE+-]+|NA)", open(f"{d}/_logs/{f}").read())
            rows.append((f[:-4], float(m.group(1)) if m and m.group(1) != "NA" else np.nan))
    return pd.Series(dict(rows), name="lambda")


os.makedirs(OUT, exist_ok=True)
for tag in ("noambig", "primary"):
    s = lambdas(f"~/work/primary_v2_{tag}").dropna().sort_index()
    s.to_frame().to_csv(os.path.join(OUT, f"lambda_merged_{tag}.csv"))
    print(f"== {tag}: {len(s)} loci; median {s.median():.4g}; > 0.01 at {(s > 0.01).sum()}")
