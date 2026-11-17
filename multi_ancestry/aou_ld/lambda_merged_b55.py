"""SuSiE-RSS lambda (susieR::estimate_s_rss on core z + merged LD, as recorded by SuSiE-merged) per locus:
v2 noambig / v2 primary (this run) vs v1 (b59d67a: B30 primary, B34 noambig). Hypothesis (2026-09-30):
CARMA not finishing at chr3 may be z/LD mismatch in the v2 merged LD."""
import gzip, io, os, re, subprocess, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
def v1(path):
    """v1 comparison column only (not used in the paper): read from the original repository's history at b59d67a.
    That history is not part of this release, so the column comes back empty and only lambda_v2 is filled."""
    root = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True, cwd=HERE).stdout.strip()
    b = subprocess.run(["git", "show", f"b59d67a:Post-hoc_analyses/fine-mapping/multi_ancestry/aou_ld/{path}"], capture_output=True, cwd=root or HERE).stdout
    if not b:
        return pd.DataFrame({"method": pd.Series(dtype=str), "locus": pd.Series(dtype=str), "lambda": pd.Series(dtype=float)})
    return pd.read_csv(io.BytesIO(gzip.decompress(b)))
def v2(work):
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
    return pd.Series(dict(rows))
out = {}
for tag, work, p1 in (("noambig", "~/work/primary_v2_noambig", "task_b34_out/primary_status.csv.gz"),
                      ("primary", "~/work/primary_v2_primary", "task_b30_out/primary_status.csv.gz")):
    s1 = v1(p1); s1 = s1[s1.method.str.startswith("SuSiE-merged")].set_index("locus")["lambda"].astype(float)
    s2 = v2(work)
    d = pd.DataFrame({"lambda_v1": s1, "lambda_v2": s2}).dropna(how="all")
    d["ratio"] = d.lambda_v2 / d.lambda_v1
    d.to_csv(os.path.join(HERE, "task_b55_out", "lambda_merged", f"lambda_merged_{tag}.csv"))
    print(f"== {tag}: loci v1 {d.lambda_v1.notna().sum()}, v2 {d.lambda_v2.notna().sum()}")
    print(d[["lambda_v1", "lambda_v2"]].describe(percentiles=[.5, .9, .99]).T[["50%", "90%", "99%", "max"]].to_string())
    print(f"  v2 > v1 at {(d.lambda_v2 > d.lambda_v1).sum()} loci; ratio median {d.ratio.median():.2f}; v2/v1 > 2 at {(d.ratio > 2).sum()}; "
          f"v2 > 0.01 at {(d.lambda_v2 > 0.01).sum()} (v1 {(d.lambda_v1 > 0.01).sum()})")
    print("  top 8 by v2 lambda:\n" + d.sort_values("lambda_v2", ascending=False).head(8).to_string(float_format=lambda x: f"{x:.4g}"))
    print("  chr3_48182866:", d.loc["chr3_48182866_53553354"].round(5).to_dict() if "chr3_48182866_53553354" in d.index else "n/a")
