#!/usr/bin/env python3
"""export the published benchmark PIPs (SuSiEx in-sample, PGC-FINEMAP
in-sample) in Task A's tidy schema, for all 250 target loci -- not the
restricted `matched` subset compare_published_vs_reference_pip.py uses
elsewhere. Also reports each source's observed PIP floor (Task 0.3).

Output: results/benchmark/published_pip_vectors.csv.gz (method, locus, rsid, bp, a1, a2, pip)
"""
import csv
import gzip
import os
import openpyxl
import pandas as pd

REF = os.path.expanduser("~/work/ref_inputs")
OUT_DIR = os.path.expanduser("~/work/benchmark")
os.makedirs(OUT_DIR, exist_ok=True)

loci_tbl = pd.read_csv(f"{REF}/loci_250_pgc3scz.tsv", sep="\t")
index_snp_to_locus = {
    row.index_snp: f"chr{row.chr}_{row.start}_{row.stop}" for row in loci_tbl.itertuples()
}
target_snps = set(index_snp_to_locus)
print(f"{len(target_snps)} target index_snps from loci_250_pgc3scz.tsv")


def sheet_to_df(wb, sheet, header_row):
    rows = list(wb[sheet].iter_rows(values_only=True))
    df = pd.DataFrame(rows[header_row + 1:], columns=rows[header_row]).dropna(how="all")
    return df


rows_out = []

# ---------------- SuSiEx-published (Yuan et al. Suppl Table 27) ----------------
susiex_wb = openpyxl.load_workbook(f"{REF}/susiex_supp_table.xlsx", read_only=True)
t27 = sheet_to_df(susiex_wb, "Suppl Table 27", header_row=1)
t27.columns = ["index_snp", "loci", "cs_id", "rsid", "bp", "a1", "a2", "frq_a_eur",
               "frq_u_eur", "frq_a_eas", "frq_u_eas", "or_eur", "or_eas", "se_eur",
               "se_eas", "p_eur", "p_eas", "pip", "consequence", "gene",
               "gtex"] + list(t27.columns[21:])
t27["pip"] = pd.to_numeric(t27["pip"], errors="coerce")
sub = t27[t27["index_snp"].isin(target_snps)].dropna(subset=["pip"])
n_dropped = (~t27["index_snp"].isin(target_snps)).sum()
for r in sub.itertuples():
    rows_out.append(("SuSiEx-published", index_snp_to_locus[r.index_snp],
                      r.rsid, r.bp, r.a1, r.a2, r.pip))
susiex_floor = sub["pip"].min()
print(f"SuSiEx-published: {len(sub)} rows, {sub['index_snp'].nunique()} loci "
      f"(of {t27['index_snp'].nunique()} total in the sheet; {n_dropped} rows "
      f"outside our 250-locus target list dropped). Observed PIP floor: {susiex_floor}")

# ---------------- PGC-FINEMAP-published (Trubetskoy et al. ST11a) ----------------
pgc_wb = openpyxl.load_workbook(f"{REF}/pgc3_scz_supp_table11.xlsx", read_only=True)
t11a = sheet_to_df(pgc_wb, "ST11a 95% Credible Sets", header_row=0)
t11a["finemap_posterior_probability"] = pd.to_numeric(
    t11a["finemap_posterior_probability"], errors="coerce")
sub2 = t11a[t11a["index_snp"].isin(target_snps)].dropna(subset=["finemap_posterior_probability"])
n_dropped2 = (~t11a["index_snp"].isin(target_snps)).sum()
for r in sub2.itertuples():
    rows_out.append(("PGC-FINEMAP-published", index_snp_to_locus[r.index_snp],
                      r.rsid, r.position, r.GWAS_effect_allele, r.other_allele,
                      r.finemap_posterior_probability))
finemap_floor = sub2["finemap_posterior_probability"].min()
print(f"PGC-FINEMAP-published: {len(sub2)} rows, {sub2['index_snp'].nunique()} loci "
      f"(of {t11a['index_snp'].nunique()} total in the sheet; {n_dropped2} rows "
      f"outside our 250-locus target list dropped). Observed PIP floor: {finemap_floor}")

out_path = f"{OUT_DIR}/published_pip_vectors.csv.gz"
with gzip.open(out_path, "wt", newline="") as f:
    w = csv.writer(f)
    w.writerow(["method", "locus", "rsid", "bp", "a1", "a2", "pip"])
    w.writerows(rows_out)
print("wrote", out_path, f"({len(rows_out)} rows)")

print()
print("=== Task 0.3 answer ===")
print(f"SuSiEx-published (Suppl Table 27) floor: {susiex_floor} "
      f"-- table title states 'Variants with PIP > 0.1', i.e. explicitly truncated at 0.1.")
print(f"PGC-FINEMAP-published (ST11a) floor: {finemap_floor} "
      f"-- ST11a is the 95% credible set membership table, not an independent PIP>x filter; "
      f"its floor is whatever the smallest CS-member PIP happens to be, not a stated cutoff.")
