# External inputs

Files the analysis reads that were not produced by this pipeline: the benchmark (Yuan *et al.*'s SuSiEx output) and three
published supplementary tables.

The benchmark output is included, redistributed with its authors' permission. The three published tables are **not
included**: download each from its source below and save it here under the name in the first column (these names are
git-ignored). `figures/make_figures.sh` needs only the Trubetskoy table and stops with a message if it is missing; the
other two are read only by the `benchmark/` scripts that built tables already committed here.

| Save as | What it is | Source | Read by |
|---|---|---|---|
| `Trubetskoy2022_SupplementaryTable11.xlsx` | PGC3-SCZ FINEMAP results (ST11a–e); ST11d = marginal PIPs | Trubetskoy *et al.* 2022, *Nature* 604:502, Supplementary Table 11 (https://doi.org/10.1038/s41586-022-04434-5) | `figures/plot_supp_published_vs_ref.py`, `figures/plot_supp_published_benchmark.py` |
| `Yuan2024_SupplementaryTables_MOESM4.xlsx` | Yuan *et al.* supplementary tables (`41588_2024_1870_MOESM4_ESM.xlsx`) | Yuan *et al.* 2024, *Nat. Genet.* 56:1841, Supplementary Tables (https://doi.org/10.1038/s41588-024-01870-z) | `benchmark/derive_confident_stratum.py` → `inputs/loci_confident122.txt` |
| `Mawson2026_media-2_FINEMAP_full.xlsx` | Mawson *et al.* full FINEMAP credible-set table (supplementary file media-2) | Mawson *et al.* 2026, preprint (https://doi.org/10.64898/2026.05.19.26353567) | `benchmark/build_mawson_benchmark.py` → `results/benchmark/mawson_finemap_dedup.csv` |

The tables derived from these (`inputs/loci_confident122.txt`, `results/benchmark/mawson_finemap_dedup.csv`,
`results/benchmark/mawson_susie_dedup.csv`) are included. `mawson_susie_dedup.csv` was extracted by hand from Mawson
*et al.*'s SuSiE-RSS supplementary table (their media-5); see `benchmark/build_mawson_benchmark.py`.

## The benchmark

The benchmark is the full in-sample SuSiEx output of Yuan *et al.* 2024 for the PGC3-SCZ EUR + EAS analysis: the
`.snp`, `.cs` and `.summary` files for all 250 loci, zipped as `yuan2024_susiex_output.zip` (shared, and redistributed here, with the authors' permission;
the job script and logs are omitted). `figures/make_figures.sh` unpacks it to `202.EUR_ASN.SuSiEx_spec/` here
(`figures/bench_sensfdr.py`, `BENCH_RAW`). The full per-variant posterior rebuilt from it by
`benchmark/reconstruct_published_posterior.py` is `results/benchmark/pip_vectors_susiex_published_full.csv.gz`
(`BENCH`). Every concordance metric reads them.

As run: SuSiEx 1.1, `--level=0.99 --pval_thresh=1`, n = 58,749 (EUR) and 14,425 (EAS), on their no-strand-ambiguous
summary statistics and in-sample LD.
