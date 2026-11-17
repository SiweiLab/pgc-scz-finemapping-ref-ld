#!/usr/bin/env bash
# Rebuild every figure panel, composite and supplementary table from the committed results (no cluster, no refits).
# Verified: every output is pixel-identical (PNG) or exactly equal (CSV, XLSX, caption .txt) to the paper's.
#
#   bash make_figures.sh                    # writes figures/output/ (the default)
#   OUT=/tmp/fig_check bash make_figures.sh # writes elsewhere
#
# External inputs, not distributed (see ../inputs/external/README.md):
#   ../inputs/external/202.EUR_ASN.SuSiEx_spec/   Yuan et al. SuSiEx output (bench_sensfdr.BENCH_RAW; the 193 loci)
#   ../results/benchmark/pip_vectors_susiex_published_full.csv.gz   rebuilt from it (bench_sensfdr.BENCH)
#   ../inputs/external/*.xlsx                     Trubetskoy 2022 ST11, Mawson 2026 media-2
# Outside the repo: Arial Bold (macOS) or matplotlib's DejaVu Sans Bold (paths.bold_font), for panel letters.
# Repo inputs: ../results/, ../inputs/loci_*.txt, and hand_made/ (fig1A.png, "supp fig1A.png", two captions), which
# no script makes and which are copied into OUT first.
#
# Every script that reads or writes the output folder takes it from paths.V2 (env FIG_OUT), so OUT redirects all of them.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

OUT="${OUT:-output}"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
X=../results/finemap

export FIG_OUT="$OUT"   # paths.py: every script that reads/writes the output folder uses this

cp hand_made/* "$OUT/"

# ---------------------------------------------------------------- panels
# Fig 1b (and ST1 AUPRC): in-CS precision-recall, 193 loci, no-palindromic
python3 plot_sensfdr_curve.py --pr --suffix=-noambig --export $X/main --out "$OUT/fig1b_pr_curve_inCS.png"
# Supp Fig 2a source (ST1 AUPRC): all variants scored, no-palindromic
python3 plot_sensfdr_curve.py --pr --all-variants --suffix=-noambig --export $X/main --out "$OUT/fig1b_pr_curve_allvar.png"
# Supp Fig 2b source (ST1 AUPRC): in-CS, A/T and C/G kept
python3 plot_sensfdr_curve.py --pr --export $X/with_ambiguous --suffix= --out "$OUT/fig1b_pr_curve_inCS_withpal.png"
# Supp Fig 2c source (ST1 AUPRC): all variants, A/T and C/G kept
python3 plot_sensfdr_curve.py --pr --all-variants --export $X/with_ambiguous --suffix= --out "$OUT/fig1b_pr_curve_allvar_withpal.png"

# Fig 1c / Supp Fig 4 / ST1: exact lead recovery (no r2 proxy); C = native CSs, E = common purity >= 0.5
for st in 122 175 193; do
    # Fig 1c / Supp Fig 4 (scenario C, stratum $st)
    python3 recovery_rate.py --export $X/main --suffix=-noambig --stratum $st --no-proxy-rescue \
        --out "$OUT/fig1c_${st}_C.csv"
    # Fig 1c / Supp Fig 4 (scenario E, stratum $st)
    python3 recovery_rate.py --export $X/main --suffix=-noambig --stratum $st --no-proxy-rescue --purity 0.5 \
        --out "$OUT/fig1c_${st}_E.csv"
done

# Fig 1c: recovery bars, 122 loci, scenario E
python3 plot_fig1c_recovery.py --fig "$OUT" --stratum 122 \
    --note "Exports: pipeline v2, no-palindromic run (the primary analysis; A/T and C/G variants excluded from every arm, as the benchmark excludes them)." \
    --out "$OUT/fig1c_recovery.png"
# Supp Fig 4: recovery x 3 strata x 2 purity conventions (+ _group_means.csv for ST1)
python3 plot_supp_fig2_recovery.py --fig "$OUT" --out "$OUT/supp_fig2_recovery_strata.png"

# Fig 2a: chr12 ATP2A2 locus
python3 plot_fig2_locus.py --locus chr12_110129434_111489434 --locus-data ../results/fig2_loci --export $X/main \
    --suffix=-noambig --pip-dots 0.1 --out "$OUT/fig2a_chr12_110129434_111489434.png"
# Fig 2b: chr12:29.90-29.95 Mb (SuSiEx unsupported calls)
python3 plot_fig2_locus.py --locus chr12_29898169_29948269 --locus-data ../results/fig2_loci --export $X/main \
    --suffix=-noambig --pip-dots 0.1 --focus 12:29935396_C_A,rs654825 --out "$OUT/fig2b_chr12_29898169_29948269.png"
# Fig 2c: chr3:36.8 Mb (near TRANK1)
python3 plot_fig2_locus.py --locus chr3_36834083_36998583 --locus-data ../results/fig2_loci --export $X/main \
    --suffix=-noambig --pip-dots 0.1 --focus rs75968099 --out "$OUT/fig2c_chr3_36834083_36998583.png"
# Fig 2d: the 13 benchmark leads with PIP > 0.9, transposed
python3 plot_fig2d.py --export $X/main --suffix=-noambig --per-locus "$OUT/fig1c_193_C_per_locus.csv" --transpose \
    --out "$OUT/fig2d_top13_wide.png"

# ST4 / Supp Fig 11: in-CS FP/TP at PIP > 0.5 / 0.9
python3 fp_in_cs.py --export $X/main --suffix=-noambig --out "$OUT/fp_in_cs.csv"
# Supp Fig 6: published in-sample merged-LD fine-mapping vs reference LD (caption .txt is hand-written)
python3 plot_supp_published_vs_ref.py --out "$OUT/supp_published_vs_ref.png"
# Supp Fig 5: published in-sample analyses as the benchmark (one variable at a time); coverage-0.95 runs for the 0.95 benchmarks
python3 plot_supp_published_benchmark.py --out "$OUT/supp_published_benchmark.png"
# Supp Fig 7 / ST1: resolution at the 41 benchmark loci with lead PIP > 0.5
python3 plot_supp_resolution.py --export $X/main --suffix=-noambig --out "$OUT/supp_resolution_pip05.png"
# Supp Fig 3a: L sweep, in-CS
python3 plot_Lsweep_curve.py --pr --suffix=-noambig --export $X/L_sweep --primary-export $X/main \
    --out "$OUT/supp_Lsweep_pr_inCS.png"
# Supp Fig 3b: L sweep, all variants
python3 plot_Lsweep_curve.py --pr --all-variants --suffix=-noambig --export $X/L_sweep --primary-export $X/main \
    --out "$OUT/supp_Lsweep_pr_allvar.png"
# Supp Fig 9 numbers: SuSiEx UKB vs AoU EUR LD (writes supp_susiex_eur_panel.csv; needs fig1c_193_C_per_locus.csv)
python3 susiex_eur_panel_stats.py
# Supp Fig 9: SuSiEx UKB vs AoU EUR LD
python3 plot_supp_susiex_eur_panel.py --out "$OUT/supp_susiex_eur_panel.png"
# Supp Fig 11: FINEMAP/SuSiE-RSS merged vs EUR-only (needs fig1c_*_E.csv, fp_in_cs.csv)
python3 plot_supp_eur_only.py --out "$OUT/supp_eur_only.png"
# Supp Fig 10: EUR-only SuSiE-RSS under three EUR panels
python3 plot_eur_panels.py --out "$OUT/eur_panels.png"
# Supp Fig 8a/b: kriging at the Fig 2b locus (+ .csv)
python3 plot_supp_kriging.py --out "$OUT/supp_kriging_chr12_29898169.png"
# Supp Fig 8c: in-CS precision-recall with loci at per-ancestry s > 0.01 (susieR estimate_s_rss) removed; then compose 8 = a/b kriging + c
python3 plot_sensfdr_curve.py --pr --suffix=-noambig --export $X/main --s-flag 0.01 --out "$OUT/supp_pr_inCS_sflag.png"
python3 compose_supp8.py
# Supp Fig 1b: variants each method fitted
python3 plot_variant_domains.py --v2 --arms-only --out "$OUT/supp_fig1b_arm_domains.png"

# ---------------------------------------------------------------- composites
# Fig 1: a (fig1A.png, hand-made) + 1b + 1c
python3 compose_fig1.py
# Fig 2: a-c in one row, d below
python3 compose_fig2.py
# Supp Fig 2: three Fig 1b variants in one row
python3 compose_supp4.py
# Supp Fig 3: L sweep a/b
python3 compose_supp5.py
# Supp Fig 1: a ("supp fig1A.png", hand-made) + b
python3 compose_supp1.py --dir "$OUT"

# ---------------------------------------------------------------- tables
# ST3: CARMA flagged variants (.xlsx and per-locus .csv)
python3 supp_table_carma_outliers.py
# ST1-ST4: Supplementary_Tables.xlsx (reads fig1b curves, fig1c, resolution, group means, fp_in_cs, CARMA xlsx)
python3 build_supp_tables.py

echo "[done] all outputs in $OUT"
