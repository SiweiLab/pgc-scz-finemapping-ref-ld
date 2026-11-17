# External inputs

Files the analysis reads that were not produced by this pipeline.

| File | What it is | Source | Read by |
|---|---|---|---|
| `Trubetskoy2022_SupplementaryTable11.xlsx` | PGC3-SCZ FINEMAP results (ST11a–e); ST11d = marginal PIPs | Trubetskoy *et al.* 2022, *Nature* 604:502, Supplementary Table 11 | `metrics/plot_supp_published_vs_ref.py` (Supplementary Fig. 3) |
| `Yuan2024_SupplementaryTables_MOESM4.xlsx` | Yuan *et al.* supplementary tables (`41588_2024_1870_MOESM4_ESM.xlsx`) | Yuan *et al.* 2024, *Nat. Genet.*, Supplementary Tables | `metrics/derive_confident_stratum.py` → `metrics/loci_confident122.txt` |
| `Mawson2026_media-2_FINEMAP_full.xlsx` | Mawson *et al.* full FINEMAP credible-set table | Mawson *et al.* 2026 (preprint), supplementary file media-2 | `metrics/build_mawson_benchmark.py` → `metrics/mawson_finemap_dedup.csv` |

Mawson *et al.*'s SuSiE-RSS table is used as the derived `metrics/mawson_susie_dedup.csv` (their full SuSiE
supplementary table; see `build_mawson_benchmark.py`).

## Not included: the benchmark

The benchmark is the full in-sample SuSiEx output of Yuan *et al.* 2024 for the PGC3-SCZ EUR + EAS analysis
(`202.EUR_ASN.SuSiEx_spec`, shared by the authors). It is **available from the authors of Yuan *et al.* on
request**. Scripts expect it at `inputs/202.EUR_ASN.SuSiEx_spec/` (`metrics/bench_sensfdr.py`, `KAI`), which is
git-ignored so it is never committed.
