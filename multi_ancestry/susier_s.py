#!/usr/bin/env python3
"""Exact Python port of susieR::estimate_s_rss(z, R, n, r_tol=1e-8, method="null-mle"), susieR 0.14.2.

Why this file exists: ld_diagnostic.estimate_s_and_kriging, labelled "susieR estimate_s_rss" since 2026-06-12,
is NOT that function. It fits z ~ N(0, sigma^2 [(1-s)R + sI]) with a free scale sigma^2 profiled out, and it
never applies the n adjustment. susieR fits z ~ N(0, (1-s)R + sI) with no scale, after
z <- z * sqrt((n-1)/(z^2+n-2)) when n is given. This port matches susieR line for line;
check_against_r() verifies it against the installed R package.
"""
import numpy as np
from scipy.optimize import minimize_scalar


def estimate_s_rss(z, R, n=None, r_tol=1e-8):
    z = np.where(np.isfinite(z), np.asarray(z, dtype=np.float64), 0.0)
    d, V = np.linalg.eigh(np.asarray(R, dtype=np.float64))
    d = np.where(d < r_tol, 0.0, d)
    if n is not None:
        if n <= 1:
            raise ValueError("n must be greater than 1")
        z = np.sqrt((n - 1) / (z ** 2 + n - 2)) * z
    ztv = V.T @ z

    def negll(s):
        v = (1 - s) * d + s
        return 0.5 * np.sum(np.log(v)) + 0.5 * np.sum(ztv ** 2 / v)

    # R: optim(0.5, method="Brent", lower=0, upper=1); bounded Brent never evaluates the endpoints either
    return float(minimize_scalar(negll, bounds=(0.0, 1.0), method="bounded",
                                 options={"xatol": 1e-10, "maxiter": 1000}).x)


def check_against_r(seed=1):
    """Compare with susieR on synthetic loci (needs Rscript + susieR). Prints max abs difference."""
    import os
    import subprocess
    import tempfile
    rng = np.random.default_rng(seed)
    cases = []
    for p, m, mism in [(200, 300, 0.0), (300, 150, 0.05), (150, 1000, 0.2), (400, 100, 0.0)]:
        X = rng.standard_normal((m, p)).cumsum(axis=1)          # correlated, rank-deficient when m < p
        R = np.corrcoef(X, rowvar=False)
        b = np.zeros(p); b[rng.integers(p)] = 8.0
        z = R @ b + rng.multivariate_normal(np.zeros(p), R + 1e-8 * np.eye(p), method="eigh")
        z[rng.random(p) < mism] *= -1                            # allele flips
        cases.append((z, R, 50000.0))
    diffs = []
    with tempfile.TemporaryDirectory() as td:
        for i, (z, R, n) in enumerate(cases):
            np.savetxt(os.path.join(td, f"z{i}.txt"), z)
            np.savetxt(os.path.join(td, f"R{i}.txt"), R)
            for use_n in (True, False):
                nn = f", n={n}" if use_n else ""
                code = (f"suppressMessages(library(susieR)); z<-scan('{td}/z{i}.txt',quiet=TRUE); "
                        f"R<-as.matrix(read.table('{td}/R{i}.txt')); cat(sprintf('%.12g', "
                        f"suppressWarnings(estimate_s_rss(z,R{nn}))))")
                r_s = float(subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=True).stdout)
                py_s = estimate_s_rss(z, R, n if use_n else None)
                diffs.append(abs(r_s - py_s))
                print(f"case {i} n={'yes' if use_n else 'no '}  R={r_s:.6g}  py={py_s:.6g}")
    print("max |R - py| =", max(diffs))


if __name__ == "__main__":
    check_against_r()
