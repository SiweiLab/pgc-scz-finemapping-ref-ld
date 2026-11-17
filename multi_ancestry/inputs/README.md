# External inputs

Files the analysis reads that were not produced by this pipeline. None of them is redistributed here:
download each from its source and place it in this folder under the name given.

| Expected file | What it is | Source | Read by |
|---|---|---|---|
| `Trubetskoy2022_SupplementaryTable11.xlsx` | PGC3-SCZ FINEMAP results (ST11a–e); ST11d = marginal PIPs | Trubetskoy *et al.* 2022, *Nature* 604:502, Supplementary Table 11 (https://doi.org/10.1038/s41586-022-04434-5) | `metrics/plot_supp_published_vs_ref.py`, `metrics/plot_supp_published_benchmark.py` |
| `Yuan2024_SupplementaryTables_MOESM4.xlsx` | Yuan *et al.* supplementary tables (`41588_2024_1870_MOESM4_ESM.xlsx`) | Yuan *et al.* 2024, *Nat. Genet.* 56:1841, Supplementary Tables (https://doi.org/10.1038/s41588-024-01870-z) | `metrics/derive_confident_stratum.py` → `metrics/loci_confident122.txt` |
| `Mawson2026_media-2_FINEMAP_full.xlsx` | Mawson *et al.* full FINEMAP credible-set table (supplementary file media-2) | Mawson *et al.* 2026, preprint (https://doi.org/10.64898/2026.05.19.26353567) | `metrics/build_mawson_benchmark.py` → `metrics/mawson_finemap_dedup.csv` |

The tables derived from these (`metrics/loci_confident122.txt`, `metrics/mawson_finemap_dedup.csv`,
`metrics/mawson_susie_dedup.csv`) are included. `mawson_susie_dedup.csv` was extracted by hand from Mawson
*et al.*'s SuSiE-RSS supplementary table (their media-5); see `metrics/build_mawson_benchmark.py`.

## Not included: the benchmark

The benchmark is the full in-sample SuSiEx output of Yuan *et al.* 2024 for the PGC3-SCZ EUR + EAS analysis
(`202.EUR_ASN.SuSiEx_spec`), together with the full per-variant posterior rebuilt from it
(`metrics/pip_vectors_susiex_published_full.csv.gz`, made by `metrics/reconstruct_published_full_posterior.py`).
Both are **available from the authors of Yuan *et al.* on request**. Scripts expect the output at
`inputs/202.EUR_ASN.SuSiEx_spec/` (`metrics/bench_sensfdr.py`, `KAI`) and the rebuilt posterior at
`metrics/pip_vectors_susiex_published_full.csv.gz`; every concordance metric reads them.
