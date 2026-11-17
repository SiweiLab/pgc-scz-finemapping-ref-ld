# Tool versions — as run, and latest (2026-09-23)

Every version and default here was read from the tool itself: the binary's output,
`formals()` on the installed package, or source at the pinned commit. None comes from docs.
Constants live in `pipeline/run_params.py` (the "Version pins" section).

## Versions

| tool | as run (original + all B-tasks) | latest upstream, 2026-09-23 | same? | rerun uses |
|---|---|---|---|---|
| FINEMAP | **v1.4.2** (md5 `b5741e43…`, official build) | v1.4.2 (christianbenner.com) | ✓ | v1.4.2 |
| susieR | **0.14.2** (CRAN) | CRAN 0.14.2; GitHub master 0.16.6 (`8e56a8e`, unreleased, API rewrite) | ✓ (release) | 0.14.2 — see below |
| MultiSuSiE | `10351eb` (pip from GitHub, no tags) | `10351eb` (last commit 2026-02-15) | ✓ | `10351eb` |
| MESuSiE | `5a3dd32` (install_github, no tags) | `5a3dd32` (last commit 2025-08-02) | ✓ | `5a3dd32` |
| SuSiEx | `v1.1.2-3-g1db2f58` | `1db2f58` | ✓ | `1db2f58` |
| RSparsePro | `5990896` (git clone 2026-09-08) | `5990896` (last commit 2026-06-23) | ✓ | `5990896` |
| CARMA | `856f0fb`, R 4.5.3 conda env (pilot only) | `856f0fb` | ✓ | `856f0fb` |
| MACHINE | `f5abc5b` (pilots only; not in the paper) | `f5abc5b` | ✓ | — |

For the original runs, "as run" is inferred for the untagged GitHub tools: upstream has no
commit between that tool's last commit and our install dates (checked 2026-09-08
against the archived July installs), so HEAD then was HEAD now. **No tool's installed version differs
from the latest release**, so "use the latest" changes no binary or package for the rerun.

## Two corrections to earlier records

- **FINEMAP was v1.4.2 throughout, not v1.4.1 (an earlier note said v1.4.1).** v1.4.2 prints `v1.4.2` on
  stdout but writes `Welcome to FINEMAP v1.4.1` into its `.log_sss` file; the binary contains
  both strings. That note read the file banners. The surviving stdout captures (qcexcl
  `.finemap.log`) and every install record say v1.4.2.
- **FINEMAP `--n-conv-sss` defaults to 100 in both versions, not 1000.** Both binaries' `--help`
  says "Default is 1000", but both run with 100 when the flag is omitted (SETTINGS block, their
  own example data). Every run used 100. So "10× stricter" was the help text, not the code.
  Pinned now as `N_CONV_SSS_FINEMAP = 100`, to be passed explicitly.
  **The v2 runs reported in the paper pass `--n-conv-sss 1000` explicitly** (`N_CONV_SSS_FINEMAP_PRIMARY` in `pipeline/run_params.py`); the 100 above describes the earlier runs.

## Defaults the drivers never pass (as run)

| tool | driver passes | unpassed defaults that applied |
|---|---|---|
| FINEMAP | `--sss --n-causal-snps 5 [--prob-cred-set 0.99]` | `n-conv-sss 100`, `prob-conv-sss-tol 0.001`, `n-iter 100000`, `corr-config 0.95`, `prior-std 0.05`, `n-configs-top 50000`, `prior-k0 0`; **no seed exists** |
| susieR `susie_rss` | `z, R, n, L=10, coverage, min_abs_corr` | via `susie_suff_stat`: `scaled_prior_variance 0.2`, `estimate_prior_method "optim"`, `max_iter 100`, `tol 1e-3`, `refine FALSE`, `n_purity 100`, `estimate_residual_variance FALSE`, `check_prior TRUE`, `z_ld_weight 0` |
| MultiSuSiE `multisusie_rss` | `L=10, coverage=0.99, min_abs_corr=0.5, z_list` | `rho 0.75`, `estimate_prior_method "early_EM"`, `pop_spec_effect_priors True`, `iter_before_zeroing_effects 5`, `max_iter 100`, `tol 1e-3`, **`float_type float32`**, `single_population_mac_thresh 20`. Library defaults are `coverage 0.95`, **`min_abs_corr 0`** — 0.5 comes only from our wrapper |
| MESuSiE `meSuSie_core` | `L=10` | `optim_method "optim"`, `estimate_residual_variance FALSE`, `max_iter 100`, `cor_method "min.abs.corr"`, `cor_threshold 0.5`, `ancestry_weight NULL`. `meSuSie_get_cs` still hardcodes coverage 0.95 upstream → `get_cs_fixed` still needed |
| SuSiEx | `n_sig 5, level 0.99, pval_thresh 1, mult-step True, keep-ambig False` | `maf 0.005`, `min_purity 0.5`, `tol 1e-4`, `max_iter 100`, `threads 1` |
| RSparsePro | `K 10, cthres 0.95, eincre 1.5` | `maxite 100`, `eps 1e-5`, `ubound 100000`, `minldthres 0.7`, `maxldthres 0.2`, `varemax 100`, `varemin 1e-3` |
| CARMA | `z.list, ld.list, lambda.list = 1/sqrt(p), output.labels` | all of `formals(CARMA)`; `lambda.list` is **required** despite a `NULL` default |
| single-anc SuSiE-RSS (`susie/run_susie_rss.R`) | hardcoded `L=10, coverage=0.95, min_abs_corr=0.1` | as susieR above. Note the 0.95/0.1 **differs from the 0.99/0.5** of every other arm |

## Defaults that differ between versions

- **FINEMAP v1.4.1 → v1.4.2**: none. Identical option lists, identical effective SETTINGS.
- **susieR 0.14.2 → 0.16.6 (unreleased master)**: many. `susie_rss` loses `prior_variance`
  (50), `z_ld_weight` and `check_prior`, and gains `scaled_prior_variance` 0.2 and
  `R_mismatch`, `unmappable_effects`, `L_greedy` and `convergence_method`. `susie` changes
  `tol` 1e-3→NULL, `n_purity` 100→"auto" and `residual_variance_lowerbound`. Not a drop-in
  upgrade. Our drivers' `check_prior = FALSE` fallback would error.
- **RSparsePro**: its 2026 source change (`34bef1b`) left every argparse default untouched.
- MultiSuSiE, MESuSiE, SuSiEx, CARMA, MACHINE: a single version in play, so nothing to diff.
