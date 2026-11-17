#!/usr/bin/env bash
# Regenerate every v2 manuscript figure/table in figures/v2/ from the laptop-side exports (no cluster, no refits).
# Reconstructed 2026-10-06 and verified: with OUT=/tmp/regen_check every output below is
# pixel-identical (PNG) / exactly equal (CSV, XLSX, caption .txt) to the committed figures/v2 files.
#
#   bash regenerate_v2.sh                                   # writes figures/v2 (the default)
#   OUT=/tmp/regen_check bash regenerate_v2.sh   # writes elsewhere, repo untouched
#
# External inputs: ../inputs/ (Trubetskoy 2022 ST11, Mawson 2026 media-2; see ../inputs/README.md), including
#   ../inputs/202.EUR_ASN.SuSiEx_spec/                           Yuan et al. SuSiEx output (bench_sensfdr.KAI; the 193 loci).
#                                                                 Git-ignored, not distributed: available from the authors of Yuan et al. on request.
# Outside the repo:
#   Arial Bold (macOS) or matplotlib DejaVu Sans Bold (paths.bold_font)   panel letters in the compose_*.py composites
# Repo inputs: ../aou_ld/task_v2_out/{noambig,primary,lsweep/noambig,extra_susiex_euraou}, ../aou_ld/task_{a,b53,b55,
#   b57,b64,b67,b70}_out, and metrics/*.csv(.gz), loci_*.txt.
# Hand-made inputs that no script makes (copied into OUT when OUT is not figures/v2):
#   figures/v2/fig1A.png (Fig 1a), figures/v2/supp fig1A.png (Supp Fig 1a), figures/v2/supp_published_{vs_ref,benchmark}_caption.txt
#
# Scripts that read or write the v2 folder take it from paths.V2 (env V2_OUT), so OUT redirects all of them.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

OUT="${OUT:-figures/v2}"
mkdir -p "$OUT"
V2="$(cd figures/v2 && pwd)"
OUT="$(cd "$OUT" && pwd)"
X=../aou_ld/task_v2_out

export V2_OUT="$OUT"   # paths.py: every script that reads/writes the v2 folder uses this

if [ "$OUT" != "$V2" ]; then
    cp "$V2/fig1A.png" "$V2/supp fig1A.png" "$V2/supp_published_vs_ref_caption.txt" "$V2/supp_published_benchmark_caption.txt" "$OUT/"
fi

# ---------------------------------------------------------------- panels
# Fig 1b (and ST1 AUPRC): in-CS precision-recall, 193 loci, no-palindromic
python3 plot_sensfdr_curve.py --pr --suffix=-noambig --export $X/noambig --out "$OUT/fig1b_pr_curve_inCS.png"
# Supp Fig 5a source (ST1 AUPRC): all variants scored, no-palindromic
python3 plot_sensfdr_curve.py --pr --all-variants --suffix=-noambig --export $X/noambig --out "$OUT/fig1b_pr_curve_allvar.png"
# Supp Fig 5b source (ST1 AUPRC): in-CS, A/T and C/G kept
python3 plot_sensfdr_curve.py --pr --export $X/primary --suffix= --out "$OUT/fig1b_pr_curve_inCS_withpal.png"
# Supp Fig 5c source (ST1 AUPRC): all variants, A/T and C/G kept
python3 plot_sensfdr_curve.py --pr --all-variants --export $X/primary --suffix= --out "$OUT/fig1b_pr_curve_allvar_withpal.png"

# Fig 1c / Supp Fig 2 / ST1: exact lead recovery (no r2 proxy); C = native CSs, E = common purity >= 0.5
for st in 122 175 193; do
    # Fig 1c / Supp Fig 2 (scenario C, stratum $st)
    python3 recovery_rate.py --export $X/noambig --suffix=-noambig --stratum $st --no-proxy-rescue \
        --out "$OUT/fig1c_${st}_C.csv"
    # Fig 1c / Supp Fig 2 (scenario E, stratum $st)
    python3 recovery_rate.py --export $X/noambig --suffix=-noambig --stratum $st --no-proxy-rescue --purity 0.5 \
        --out "$OUT/fig1c_${st}_E.csv"
done

# Fig 1c: recovery bars, 122 loci, scenario E
python3 plot_fig1c_recovery.py --fig "$OUT" --stratum 122 \
    --note "Exports: pipeline v2, no-palindromic run (the primary analysis; A/T and C/G variants excluded from every arm, as the benchmark excludes them)." \
    --out "$OUT/fig1c_recovery.png"
# Supp Fig 2: recovery x 3 strata x 2 purity conventions (+ _group_means.csv for ST1)
python3 plot_supp_fig2_recovery.py --fig "$OUT" --out "$OUT/supp_fig2_recovery_strata.png"

# Fig 2a: chr12 ATP2A2 locus
python3 plot_fig2_locus.py --locus chr12_110129434_111489434 --b32 ../aou_ld/task_b64_out --export $X/noambig \
    --suffix=-noambig --pip-dots 0.1 --out "$OUT/fig2a_chr12_110129434_111489434.png"
# Fig 2b: chr12:29.90-29.95 Mb (SuSiEx unsupported calls)
python3 plot_fig2_locus.py --locus chr12_29898169_29948269 --b32 ../aou_ld/task_b64_out --export $X/noambig \
    --suffix=-noambig --pip-dots 0.1 --focus 12:29935396_C_A,rs654825 --out "$OUT/fig2b_chr12_29898169_29948269.png"
# Fig 2c: chr3:36.8 Mb (near TRANK1)
python3 plot_fig2_locus.py --locus chr3_36834083_36998583 --b32 ../aou_ld/task_b64_out --export $X/noambig \
    --suffix=-noambig --pip-dots 0.1 --focus rs75968099 --out "$OUT/fig2c_chr3_36834083_36998583.png"
# Fig 2d: the 13 benchmark leads with PIP > 0.9, transposed
python3 plot_fig2d.py --export $X/noambig --suffix=-noambig --per-locus "$OUT/fig1c_193_C_per_locus.csv" --transpose \
    --out "$OUT/fig2d_top13_wide.png"

# ST2 / Supp Fig 8: in-CS FP/TP at PIP > 0.5 / 0.9
python3 fp_in_cs_v2.py --export $X/noambig --suffix=-noambig --out "$OUT/fp_in_cs_v2.csv"
# Supp Fig 3: published in-sample merged-LD fine-mapping vs reference LD (caption .txt is hand-written)
python3 plot_supp_published_vs_ref.py --out "$OUT/supp_published_vs_ref.png"
# published in-sample analyses as the benchmark (one variable at a time); B69 coverage-0.95 runs for the 0.95 benchmarks
python3 plot_supp_published_benchmark.py --out "$OUT/supp_published_benchmark.png"
# Supp Fig 4 / ST1: resolution at the 41 benchmark loci with lead PIP > 0.5
python3 plot_supp_resolution.py --export $X/noambig --suffix=-noambig --out "$OUT/supp_resolution_pip05.png"
# Supp Fig 6a: L sweep, in-CS
python3 plot_Lsweep_curve.py --pr --suffix=-noambig --export $X/lsweep/noambig --primary-export $X/noambig \
    --out "$OUT/supp_Lsweep_pr_inCS.png"
# Supp Fig 6b: L sweep, all variants
python3 plot_Lsweep_curve.py --pr --all-variants --suffix=-noambig --export $X/lsweep/noambig --primary-export $X/noambig \
    --out "$OUT/supp_Lsweep_pr_allvar.png"
# Supp Fig 7 numbers: SuSiEx UKB vs AoU EUR LD (writes supp_susiex_eur_panel.csv; needs fig1c_193_C_per_locus.csv)
python3 supp_susiex_eur_panel.py
# Supp Fig 7: SuSiEx UKB vs AoU EUR LD
python3 plot_supp_susiex_eur_panel.py --out "$OUT/supp_susiex_eur_panel.png"
# Supp Fig 8: FINEMAP/SuSiE-RSS merged vs EUR-only (needs fig1c_*_E.csv, fp_in_cs_v2.csv)
python3 plot_supp_eur_only.py --out "$OUT/supp_eur_only.png"
# Supp Fig 9: EUR-only SuSiE-RSS under three EUR panels
python3 plot_eur_panels_b53.py --b53-dir ../aou_ld/task_b53_out --s-table ../aou_ld/task_b71_out/s_stats_susier.csv --out "$OUT/eur_panels_b53.png"
# Supp Fig 8a/b: kriging at the Fig 2b locus (+ .csv)
python3 plot_supp_kriging.py --out "$OUT/supp_kriging_chr12_29898169.png"
# Supp Fig 8c: in-CS precision-recall with loci at per-ancestry s > 0.01 (B70, susieR estimate_s_rss) removed; then compose 8 = a/b kriging + c
python3 plot_sensfdr_curve.py --pr --suffix=-noambig --export $X/noambig --s-flag 0.01 --out "$OUT/supp_pr_inCS_sflag.png"
python3 compose_supp8.py
# Supp Fig 1b: variants each method fitted
python3 plot_variant_domains.py --v2 --arms-only --out "$OUT/supp_fig1b_arm_domains.png"

# ---------------------------------------------------------------- composites
# Fig 1: a (fig1A.png, hand-made) + 1b + 1c
python3 compose_fig1.py
# Fig 2: a-c in one row, d below
python3 compose_fig2.py
# Supp Fig 5: three Fig 1b variants in one row
python3 compose_supp4.py
# Supp Fig 6: L sweep a/b
python3 compose_supp5.py
# Supp Fig 1: a ("supp fig1A.png", hand-made) + b
python3 compose_supp1.py --dir "$OUT"

# ---------------------------------------------------------------- tables
# ST3: CARMA flagged variants (.xlsx and per-locus .csv)
python3 supp_table_carma_outliers.py
# ST1-ST4: Supplementary_Tables.xlsx (reads fig1b curves, fig1c, resolution, group means, fp_in_cs, CARMA xlsx)
python3 build_supp_tables.py

# ---------------------------------------------------------------- assembly (writes into figures/v2/supp and
# figures/v2/manuscript, always under figures/v2 regardless of OUT -- run by hand once figures/v2 itself is final)
# python3 assemble_manuscript.py

echo "[done] all outputs in $OUT"
