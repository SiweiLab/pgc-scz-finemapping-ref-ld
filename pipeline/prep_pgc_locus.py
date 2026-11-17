#!/usr/bin/env python3
"""
Prepare per-locus sumstats from a PGC sumstats VCF-TSV for the GCTB-LD SuSiE pipeline.

The PGC3-SCZ wave3 public release uses a VCF-style TSV (## meta lines then a
tab-separated table) with columns:
    CHROM  ID  POS  A1  A2  FCAS  FCON  IMPINFO  BETA  SE  PVAL  NCAS  NCON  NEFF
where A1 is the effect allele and BETA is the log-OR of A1.

This writes, per locus, a {chrN_start_end}_sumstats.tsv.gz with the columns the
downstream tools expect:
    SNP_ID  MarkerName  Allele1  Allele2  Zscore  Weight  POS
    - SNP_ID     = chrN:POS:A1:A2
    - MarkerName = rsID (PGC ID)            <- matched to GCTB snp.info ID
    - Allele1    = A1 (effect allele)       <- Zscore aligned to this allele
    - Allele2    = A2
    - Zscore     = BETA / SE
    - Weight     = NEFF (effective N)       <- used as susie_rss n=

build_ukb_ld_matrix.py flip-corrects the GCTB R matrix to Allele1, and
run_susie_gctb.R uses the original Zscore, so the Allele1 = A1 convention here is
consistent end to end.

Usage:
  python prep_pgc_locus.py loci.bed pgc_sumstats.vcf.tsv.gz output_dir
      [--build-prefix chr]     # prepend "chr" to CHROM in SNP_ID/output tags (default: yes)
"""
import sys
import os
import argparse
import collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gcs_io import makedirs as gcs_makedirs, open_maybe_gz  # noqa: E402


def open_file(f):
    return open_maybe_gz(f, "rt", encoding="utf-8", errors="replace")


def load_loci(bed_file):
    loci_by_chrom = collections.defaultdict(list)
    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            chrom = parts[0]
            chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
            start, end = int(parts[1]), int(parts[2])
            loci_by_chrom[chrom].append((start, end))
    return loci_by_chrom


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("loci_bed",   help="BED file with loci (chr start end)")
    parser.add_argument("sumstats",   help="PGC sumstats VCF-TSV (gzipped or plain)")
    parser.add_argument("output_dir", help="Directory to save per-locus TSVs")
    args = parser.parse_args()

    gcs_makedirs(args.output_dir)

    loci_by_chrom = load_loci(args.loci_bed)
    n_loci = sum(len(v) for v in loci_by_chrom.values())
    print(f"Loaded {n_loci} loci across {len(loci_by_chrom)} chromosomes")

    out_cols = ["SNP_ID", "MarkerName", "Allele1", "Allele2", "Zscore", "Weight", "Freq1", "POS"]
    loci_rows = {(chrom, s, e): [] for chrom, loci in loci_by_chrom.items() for s, e in loci}

    n_read = n_kept = n_skipped_se = 0
    with open_file(args.sumstats) as fh:
        # Skip ## meta lines, then the next line is the column header.
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            header = line.rstrip("\n").split("\t")
            break
        if header is None:
            sys.exit("ERROR: no header line found (only ## meta lines?)")

        idx = {c: header.index(c) for c in
               ("CHROM", "ID", "POS", "A1", "A2", "BETA", "SE", "NEFF",
                "FCAS", "FCON", "NCAS", "NCON")
               if c in header}
        required = {"CHROM", "ID", "POS", "A1", "A2", "BETA", "SE", "NEFF"}
        missing = required - set(idx)
        if missing:
            sys.exit(f"ERROR: sumstats missing columns {missing}. Found: {header}")
        # Freq1 = sample-size-weighted A1 frequency from case/control freqs (for FINEMAP maf).
        have_freq = {"FCAS", "FCON", "NCAS", "NCON"} <= set(idx)

        for line in fh:
            n_read += 1
            cols = line.rstrip("\n").split("\t")
            chrom_val = cols[idx["CHROM"]]
            chrom = chrom_val if chrom_val.startswith("chr") else f"chr{chrom_val}"
            if chrom not in loci_by_chrom:
                continue
            pos = int(cols[idx["POS"]])
            hit = [(chrom, s, e) for (s, e) in loci_by_chrom[chrom] if s <= pos <= e]
            if not hit:
                continue

            a1 = cols[idx["A1"]].upper()
            a2 = cols[idx["A2"]].upper()
            try:
                beta = float(cols[idx["BETA"]])
                se = float(cols[idx["SE"]])
            except ValueError:
                continue
            if se <= 0:
                n_skipped_se += 1
                continue
            z = beta / se
            weight = cols[idx["NEFF"]]
            rsid = cols[idx["ID"]]
            snp_id = f"{chrom}:{pos}:{a1}:{a2}"

            freq1 = "NA"
            if have_freq:
                try:
                    ncas = float(cols[idx["NCAS"]]); ncon = float(cols[idx["NCON"]])
                    fcas = float(cols[idx["FCAS"]]); fcon = float(cols[idx["FCON"]])
                    if ncas + ncon > 0:
                        freq1 = f"{(fcas * ncas + fcon * ncon) / (ncas + ncon):.6f}"
                except (ValueError, IndexError):
                    pass

            rec = [snp_id, rsid, a1, a2, f"{z:.6f}", weight, freq1, str(pos)]
            for key in hit:
                loci_rows[key].append(rec)
            n_kept += 1

    print(f"Read {n_read:,} data rows; kept {n_kept:,} in loci; "
          f"skipped {n_skipped_se:,} with SE<=0")

    summary = []
    for (chrom, start, end), rows in loci_rows.items():
        tag = f"{chrom}_{start}_{end}"
        out_path = os.path.join(args.output_dir, f"{tag}_sumstats.tsv.gz")
        with open_maybe_gz(out_path, "wt") as out:
            out.write("\t".join(out_cols) + "\n")
            for rec in rows:
                out.write("\t".join(rec) + "\n")
        print(f"{chrom}:{start}-{end}  ->  {len(rows)} SNPs  ({out_path})")
        summary.append((f"{chrom}:{start}-{end}", len(rows)))

    summary_path = os.path.join(args.output_dir, "locus_snp_counts.tsv")
    with open_maybe_gz(summary_path, "wt") as sf:
        sf.write("locus\tn_snps\n")
        for loc, n in summary:
            sf.write(f"{loc}\t{n}\n")
    print(f"\nSummary: {summary_path}")


if __name__ == "__main__":
    main()
