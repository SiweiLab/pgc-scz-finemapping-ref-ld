#!/usr/bin/env python3
# ============================================================================
# STEP 0 — RUN INSIDE THE All of Us RESEARCHER WORKBENCH (Hail / Spark env)
# ============================================================================
# Produce the per-locus EAS in-sample LD directly in Hail, straight from the
# CDR v9 srWGS Hail VDS:
#   - {locus}.ld   (n x n signed r, whitespace-delimited -- read by build_eas_ld_npz.py)
#   - {locus}.z    (rsid chromosome position allele1 allele2, same row order as .ld)
#   - {locus}_eas_ids.incl  (EAS research IDs, 1/line; record of exactly who went into R)
#
# Pipeline: v9 ancestry_preds -> EAS sample set -> filter VDS to (EAS x interval)
# -> split multiallelics -> densify to MT -> (optional) restrict to the EUR SNP
# set -> hl.row_correlation(GT dosage) -> .ld + .z.
#
# >>> WHY NOT BGEN/LDstore2 <<<
#   AoU srWGS variant calls are hard genotype calls (GT), not imputed dosages.
#   hl.export_bgen refuses hard-call-only MatrixTables outright ("exporting to
#   BGEN requires a GP ... If you only have hard calls (GT), BGEN is probably
#   not the right format") and LDstore2 v2.0's only input modes are
#   --read-only-bgen / --read-bdose (bdose itself only comes from a BGEN) --
#   there is no hard-call ingestion path in LDstore2 at all. MultiSuSiE's
#   published LDstore2 recipe assumes genuinely probabilistic (imputed)
#   genotypes, e.g. UKB's imputed BGEN -- this repo's own completed EUR pilot
#   (build_ukb_ld_matrix.py) confirms that: it extracts sub-matrices from a
#   *precomputed* public GCTB SBayesR UKB LD reference panel, it never ran
#   LDstore2 itself. AoU's EAS cohort has no equivalent precomputed panel, so
#   in-sample LD has to be computed here -- and since the data is hard calls,
#   hl.row_correlation is the direct fit, not a synthetic one-hot GP forced
#   into BGEN.
#
# >>> BUILD WARNING <<<
#   AoU CDR v9 srWGS is **GRCh38**. Your EUR chr12 pilot locus is **hg19**
#   (chr12_110129434_111489434). Pass this step a GRCh38 interval and, for a
#   mergeable joint locus, a --target-snps file = the EUR .z SNPs lifted to
#   GRCh38 (rsid + hg38 chr/pos/alleles). The EAS R produced downstream is then
#   aligned to EUR **by rsID**, not by hg19 position. Do NOT reuse hg19 coords.
#
# Env vars (verify exact names in the Controlled CDR Directory for your CDR):
#   WGS_VDS_PATH   -> gs:// path to the srWGS Hail VDS   (used if --vds omitted)
#   The v9 ancestry_preds TSV path is NOT a stable env var across CDRs; pass it
#   with --ancestry-file (look it up in the Controlled CDR Directory / auxiliary
#   'ancestry' folder; expected columns include research_id + ancestry_pred).
#   v9 ancestry_preds.tsv has TWO ancestry columns: 'ancestry_pred' (raw hard
#   argmax, no "remaining" bucket) and 'ancestry_pred_other' (adds a 7th "oth"
#   group for samples not confidently assigned). AoU's official CDR v9 QC report
#   (Table G.1 genetic ancestry group counts) matches ancestry_pred_other exactly
#   (eas=12,980), not ancestry_pred (eas=15,893) — this script prefers
#   ancestry_pred_other when present.
#
# Usage (in a Hail-enabled notebook / dataproc cluster on the Workbench):
#   python hail_ld.py \
#       --locus chr12_110129434_111489434 \
#       --interval chr12:109000000-113000000 \        # GRCh38 window (lift the EUR window!)
#       --ancestry-file gs://.../ancestry_preds.tsv \
#       --out <out_dir>/chr12 \
#       [--vds gs://.../srwgs.vds]  [--target-snps eur_locus_hg38.tsv] \
#       [--ancestry-code eas] [--min-call-rate 0.95] [--biallelic-snv-only]
#
# Then: python build_eas_ld_npz.py <out>/{locus}.ld <out>/{locus}.z --locus {locus} --out <dir>
# ----------------------------------------------------------------------------
import os
import argparse
import numpy as np
import hail as hl


def load_eas_ids(ancestry_file, code):
    """v9 ancestry_preds -> set of EAS research IDs (as strings)."""
    ht = hl.import_table(ancestry_file, impute=True, key=None)
    cols = set(ht.row)
    id_col = next((c for c in ("research_id", "s", "sample_id", "person_id") if c in cols), None)
    # Prefer ancestry_pred_other: v9 adds a 7th "oth" (unassigned/admixed) group that
    # ancestry_pred's raw argmax doesn't have. ancestry_pred_other matches AoU's
    # official CDR v9 QC report (Table G.1) exactly; ancestry_pred does not.
    anc_col = next((c for c in ("ancestry_pred_other", "ancestry_pred", "ancestry", "predicted_ancestry")
                    if c in cols), None)
    if id_col is None or anc_col is None:
        raise SystemExit(f"ancestry file missing id/ancestry columns; saw {sorted(cols)}")
    ht = ht.filter(hl.str(ht[anc_col]).lower() == code.lower())
    ids = ht.aggregate(hl.agg.collect_as_set(hl.str(ht[id_col])))
    print(f"[eas] {len(ids)} samples with {anc_col} == '{code}'")
    return ids


def load_exclude_ids(path):
    """One research_id per line, NO header (e.g. a locally-built union of
    relatedness/genomic-QC/sex-concordance exclusion lists, as in
    RUN_ORDER.md step 3 -- with any header row already
    stripped by the caller before upload)."""
    ht = hl.import_table(path, no_header=True, impute=False)
    ids = ht.aggregate(hl.agg.collect_as_set(ht.f0))
    print(f"[exclude] {len(ids)} sample IDs loaded from {path}")
    return ids


def load_target_snps(path):
    """EUR SNP set lifted to GRCh38: expects columns rsid, chromosome(or chrom/chr),
    position(or pos/bp), allele1/allele2 (or a1/a2/ref/alt). Returns two Hail
    tables:
      - locus_ht: keyed by locus only -- cheap coarse restriction of the SPARSE
        VDS before densify (unchanged from the original design, for memory).
      - exact_ht: keyed by (locus, alleles), TWO rows per input SNP (both
        allele1/allele2 orderings, since we don't know a priori which one the
        reference genome calls REF vs ALT) -- lets a genuinely multiallelic
        position (real and non-negligible in AoU's srWGS panel) be resolved to
        exactly the one ALT that matches the GWAS via Hail's native keyed
        variant join (hl.vds.filter_variants) instead of a broadcast struct-
        compare in a row filter, which reproducibly stalled Hail's query
        planner/executor on this VDS at scale (confirmed: 3/3 identical runs of
        that approach hung with near-zero CPU). Without allele-aware matching
        at all, every ALT at a matched locus inherits the same rsID and
        downstream scripts (align_eas_ld_to_sumstats.py / merge_eur_eas_ld.py)
        drop the whole position as ambiguous, discarding a real match along
        with the spurious ones."""
    ht = hl.import_table(path, impute=True)
    cols = {c.lower(): c for c in ht.row}
    def pick(*names):
        for n in names:
            if n in cols:
                return cols[n]
        raise SystemExit(f"--target-snps missing one of {names}; saw {list(cols)}")
    rsid, chrom, pos = pick("rsid", "snp", "id"), pick("chromosome", "chrom", "chr"), pick("position", "pos", "bp")
    a1, a2 = pick("allele1", "a1", "ref"), pick("allele2", "a2", "alt")
    ht = ht.annotate(_contig=hl.str(ht[chrom]))
    ht = ht.annotate(_contig=hl.if_else(ht._contig.startswith("chr"), ht._contig, "chr" + ht._contig))
    ht = ht.annotate(locus=hl.locus(ht._contig, hl.int32(ht[pos]), reference_genome="GRCh38"),
                     rsid=hl.str(ht[rsid]),
                     allele1=hl.str(ht[a1]).upper(),
                     allele2=hl.str(ht[a2]).upper())
    ht = ht.select("locus", "rsid", "allele1", "allele2")
    n = ht.count()
    print(f"[target] {n} EUR SNPs (GRCh38) to restrict the EAS locus to")

    locus_ht = ht.key_by("locus").select()

    fwd = ht.annotate(alleles=[ht.allele1, ht.allele2])
    rev = ht.annotate(alleles=[ht.allele2, ht.allele1])
    exact_ht = fwd.union(rev).key_by("locus", "alleles").select("rsid")

    return locus_ht, exact_ht


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus", required=True, help="tag; keep identical across steps 0-3")
    ap.add_argument("--interval", required=True, help="GRCh38, e.g. chr12:109000000-113000000")
    ap.add_argument("--ancestry-file", required=True, help="v9 ancestry_preds TSV (gs://...)")
    ap.add_argument("--out", required=True,
                     help="MUST be a gs:// path (e.g. gs://<staging-bucket>/ld/<locus>), "
                          "not a local/gcsfuse-mounted path. All writes here go through Hail's "
                          "own FS layer (hl.hadoop_open / Table.export), which resolves a bare "
                          "absolute path against HDFS's default filesystem -- NOT the local disk "
                          "or a gcsfuse mount that happens to have the same-looking path. "
                          "Confirmed the hard way: passing a gcsfuse staging-bucket path here "
                          "silently wrote the .z/.incl files into ephemeral HDFS instead (gone "
                          "when the cluster is deleted) while only the plain-Python-written .ld "
                          "file landed in the real bucket. Always use gs://.")
    ap.add_argument("--vds", default=os.environ.get("WGS_VDS_PATH"),
                    help="srWGS Hail VDS path (defaults to $WGS_VDS_PATH)")
    ap.add_argument("--target-snps", help="EUR .z lifted to GRCh38 (restrict + carry rsID)")
    ap.add_argument("--ancestry-code", default="eas")
    ap.add_argument("--exclude-samples",
                    help="gs:// path, one research_id per line (no header) -- e.g. a "
                         "union of relatedness_flagged_samples.tsv / genomic-QC "
                         "flagged_samples.tsv / sex-concordance exceptions")
    ap.add_argument("--min-call-rate", type=float, default=0.95)
    ap.add_argument("--biallelic-snv-only", action="store_true",
                    help="keep only biallelic SNVs (recommended for a clean LD matrix)")
    ap.add_argument("--dry-run", action="store_true",
                    help="count EAS samples + interval variants on the SPARSE VDS and exit; "
                         "no densify/export (densify is the memory hog — use this to size the cluster)")
    ap.add_argument("--tmp", default="/tmp/hail_ld")
    ap.add_argument("--acaf-mt",
                    help="v2: read genotypes from AoU's precomputed, already split and "
                         "dense ACAF-threshold MatrixTable (e.g. .../acaf_threshold/splitMT/"
                         "hail.mt) instead of the VDS: no split_multi, no densify. Requires "
                         "--target-snps. Without it, the VDS path below runs as in v1.")
    ap.add_argument("--text-ld", action="store_true",
                    help="write R as %%.6g text {locus}.ld (v1's format, for exact v1 "
                         "reproduction) instead of binary {locus}.ld.npy")
    args = ap.parse_args()
    if not args.vds and not args.acaf_mt:
        raise SystemExit("no VDS path: pass --vds or ensure $WGS_VDS_PATH is set")
    if args.acaf_mt and not args.target_snps:
        raise SystemExit("--acaf-mt needs --target-snps")
    os.makedirs(args.out, exist_ok=True)

    hl.init(tmp_dir=args.tmp, default_reference="GRCh38")

    eas_ids = load_eas_ids(args.ancestry_file, args.ancestry_code)
    if args.exclude_samples:
        exclude_ids = load_exclude_ids(args.exclude_samples)
        n_before = len(eas_ids)
        eas_ids = eas_ids - exclude_ids
        print(f"[exclude] {n_before - len(eas_ids)} of {n_before} {args.ancestry_code} "
              f"samples dropped (relatedness/QC/sex-concordance); {len(eas_ids)} remain")

    interval = hl.parse_locus_interval(args.interval, reference_genome="GRCh38")

    if args.acaf_mt:
        mt = load_acaf(args, interval, eas_ids)
    else:
        mt = load_vds_dense(args, interval, eas_ids)
        if mt is None:  # --dry-run
            return
    extract(args, mt, eas_ids)


def load_acaf(args, interval, eas_ids):
    """v2: the ACAF splitMT is already split (one ALT per row, REF first) and dense, so the
    VDS block (filter_samples, both filter_variants, split_multi, to_dense_mt) is replaced by an
    interval filter, a keyed (locus, alleles) semi-join to the targets and a column filter."""
    _, exact_tgt = load_target_snps(args.target_snps)
    mt = hl.read_matrix_table(args.acaf_mt)
    # The MT's own variant_qc/info are over all ~535k samples; the EAS-subset QC is recomputed
    # below exactly as on the VDS path.
    mt = mt.drop("variant_qc", "info")
    mt = hl.filter_intervals(mt, [interval])
    mt = mt.semi_join_rows(exact_tgt)
    eas_ht = hl.Table.parallelize([hl.struct(s=x) for x in eas_ids],
                                  hl.tstruct(s=hl.tstr)).key_by("s")
    mt = mt.semi_join_cols(eas_ht)
    if args.biallelic_snv_only:
        mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    mt = mt.annotate_rows(rsid=exact_tgt[mt.row_key].rsid)
    # the ACAF MT carries filtered entries (AoU's FT entry filter) where the VDS path's
    # to_dense_mt gives missing GT. variant_qc's call_rate already counts them as not called
    # (n_filtered is in its denominator) and row_correlation unfilters internally, so neither
    # changes; but any other GT expression or collect never sees them (an earlier Stage B crash).
    # Unfiltering turns them into missing GT, so both paths agree from here.
    return mt.select_entries("GT").unfilter_entries()


def load_vds_dense(args, interval, eas_ids):
    # --- read VDS, restrict rows to the interval FIRST, THEN subset samples to EAS.
    # remove_dead_alleles=True (below) scans entries to drop alleles with zero
    # carriers among kept samples — narrowing rows first keeps that scan cheap
    # instead of running it over more of the genome than necessary.
    vds = hl.vds.read_vds(args.vds)
    vds = hl.vds.filter_intervals(vds, [interval], keep=True)
    vds = hl.vds.filter_samples(
        vds, hl.Table.parallelize(
            [hl.struct(s=x) for x in eas_ids], hl.tstruct(s=hl.tstr)).key_by("s"),
        keep=True, remove_dead_alleles=True)

    # Restrict to the --target-snps loci on the SPARSE VDS, BEFORE densify. Densify
    # materializes genotypes for every row x sample it's given, so densifying the
    # full interval and only afterward dropping rows outside --target-snps (the old
    # order) was the OOM: it densified ~87k interval variants x ~13k EAS samples
    # when only ~2.3k target loci are ever kept. filter_variants narrows the rows
    # densify has to touch to just those loci first.
    locus_tgt, exact_tgt = load_target_snps(args.target_snps) if args.target_snps else (None, None)
    if locus_tgt is not None:
        vds = hl.vds.filter_variants(vds, locus_tgt, keep=True)

    # Split multiallelics on the SPARSE VDS, BEFORE densify. This VDS's variant_data uses
    # Hail's local-allele sparse encoding (entries are LGT/LA/LAD, not legacy GT/AD/PL) --
    # hl.split_multi_hts only knows how to recode a legacy GT/AD/PL schema; on LGT/LA
    # entries it silently passes them through *unrecoded* per split row instead of
    # erroring (confirmed: entries stayed LGT/LA and every downstream GT-based call, e.g.
    # variant_qc, then failed with "expect an entry field 'GT'"). hl.vds.split_multi is
    # local-allele aware and produces real GT/AD entries directly -- verified against this
    # VDS's schema before relying on it.
    vds = hl.vds.split_multi(vds)

    # Exact (locus, alleles) restriction, still on the SPARSE VDS (before densify,
    # same memory rationale as the coarse locus_tgt filter above). Needs to happen
    # AFTER split_multi since alleles are only resolved to one-ALT-per-row post-split.
    # Using Hail's own keyed variant join (filter_variants) here -- not a broadcast
    # struct-compare in filter_rows -- is what actually fixes the genuinely-
    # multiallelic-position ambiguity without stalling (see load_target_snps).
    if exact_tgt is not None:
        vds = hl.vds.filter_variants(vds, exact_tgt, keep=True)

    # Cheap sanity check on the SPARSE VDS — no densify (densify materializes
    # genotypes for all samples x variants and is what blows up memory).
    if args.dry_run:
        n_var = vds.variant_data.count_rows()
        n_samp = vds.variant_data.count_cols()
        print(f"[dry-run {args.locus}] interval {args.interval}: "
              f"{n_var} variants x {n_samp} EAS samples (sparse; nothing exported)")
        return None

    # densify to a dense MT over the (now target-SNP-restricted, already-split) interval
    mt = hl.vds.to_dense_mt(vds)

    if args.biallelic_snv_only:
        mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

    # rsID: carry EUR rsID when --target-snps given; else fall back to variant id.
    # By this point the row set is ALREADY exactly (locus, alleles)-restricted via
    # the vds-level filter_variants(exact_tgt) call above, using Hail's own keyed
    # variant join -- not a broadcast struct-compare in filter_rows, which
    # reproducibly stalled Hail's query planner/executor on this VDS at scale
    # (confirmed: 3/3 identical runs hung with near-zero CPU across different
    # predicate expressions -- hl.set() equality and an explicit boolean OR both
    # stalled, at different pipeline stages each time; the join itself was the
    # problem, not the comparison expression). So this is now a clean, cheap 1:1
    # lookup (both allele orderings in exact_tgt carry the same rsid, so either
    # matches), not a filter.
    if exact_tgt is not None:
        mt = mt.annotate_rows(rsid=exact_tgt[mt.row_key].rsid)
        mt = mt.filter_rows(hl.is_defined(mt.rsid))
    else:
        mt = mt.annotate_rows(
            rsid=mt.locus.contig.replace("^chr", "") + ":" + hl.str(mt.locus.position)
                 + ":" + mt.alleles[0] + ":" + mt.alleles[1])
    return mt


def extract(args, mt, eas_ids):
    # basic QC over the EAS subset. Checkpoint BEFORE the call-rate and monomorphic filters,
    # so the drops can be counted from it cheaply; the checkpoint also fixes one physical row
    # order (mt's key order) for both row_correlation and the .z export, and stops Hail
    # re-running the whole read/filter plan for each action that follows.
    mt = hl.variant_qc(mt)
    mt = mt.checkpoint(os.path.join(args.tmp, f"{args.locus}_eas_qc.mt"), overwrite=True)
    n_in, n_samp = mt.count()
    cr_ok = mt.variant_qc.call_rate >= args.min_call_rate
    # Monomorphic in the EAS subset (AF 0 or 1): row_correlation gives them NaN rows, which
    # build_eas_ld_npz.py used to drop afterwards. Drop them here instead.
    poly = (mt.variant_qc.AF[1] > 0) & (mt.variant_qc.AF[1] < 1)
    n_cr, n_mono = mt.aggregate_rows((hl.agg.count_where(~cr_ok), hl.agg.count_where(cr_ok & ~poly)))
    mt = mt.filter_rows(cr_ok & poly)
    n_var = n_in - n_cr - n_mono
    print(f"[qc {args.locus}] {n_in} matched rows; dropped call_rate<{args.min_call_rate}: {n_cr}; "
          f"dropped monomorphic: {n_mono}")
    print(f"[locus {args.locus}] {n_var} variants x {n_samp} EAS samples after QC")

    # --- in-sample EAS LD, computed directly in Hail (see BUILD WARNING at top for why
    # not BGEN/LDstore2) ---
    # n_alt_alleles() = ALT dosage, matching allele1=ALT below -- R is already signed to
    # the same effect allele as the .z, no extra sign-flipping needed to align them.
    R = hl.row_correlation(mt.GT.n_alt_alleles()).to_numpy()
    # hl.hadoop_open (not a plain Python open()) -- see --out help: this must go through
    # Hail's own FS layer, same as the .z/.incl writes below, or it can silently land in
    # the wrong place depending on how --out resolves.
    if args.text_ld:
        ld_path = os.path.join(args.out, f"{args.locus}.ld")
        with hl.hadoop_open(ld_path, "w") as f:
            np.savetxt(f, R, fmt="%.6g")
    else:
        ld_path = os.path.join(args.out, f"{args.locus}.ld.npy")
        with hl.hadoop_open(ld_path, "wb") as f:
            np.save(f, R.astype(np.float32))

    # .z in LDstore order: rsid chromosome position allele1 allele2 (allele1 = ALT/effect),
    # same row order as R since both come from the same checkpointed mt; then the EAS-subset
    # QC: eas_af (ALT frequency), call_rate, n_called, filters (PASS if the site's
    # filter set is empty; NA if the source has no filters field).
    rows = mt.rows()
    filt = (hl.if_else(hl.len(rows.filters) == 0, "PASS", hl.delimit(hl.sorted(hl.array(rows.filters)), ";"))
            if "filters" in rows.row else hl.missing(hl.tstr))
    rows = rows.annotate(
        chromosome=rows.locus.contig.replace("^chr", ""),
        position=rows.locus.position,
        allele1=rows.alleles[1],   # ALT as effect allele; align to EAS sumstats downstream
        allele2=rows.alleles[0],   # REF
        eas_af=rows.variant_qc.AF[1],
        call_rate=rows.variant_qc.call_rate,
        n_called=rows.variant_qc.n_called,
        filters=hl.or_else(filt, "NA"))
    z = rows.select("rsid", "chromosome", "position", "allele1", "allele2",
                    "eas_af", "call_rate", "n_called", "filters")
    z_path = os.path.join(args.out, f"{args.locus}.z")
    z.export(z_path, delimiter=" ")

    # EAS include list (1 id/line) — record of exactly which samples went into R
    incl_path = os.path.join(args.out, f"{args.locus}_eas_ids.incl")
    with hl.hadoop_open(incl_path, "w") as f:
        for x in sorted(eas_ids):
            f.write(x + "\n")

    print("DONE step 0. Wrote:")
    print(f"  {ld_path}  ({R.shape[0]}x{R.shape[0]} signed r)")
    print(f"  {z_path}")
    print(f"  {incl_path}")
    print(f"  n_samples: {n_samp}")
    print("Next: python build_eas_ld_npz.py "
          f"{ld_path} {z_path} --locus {args.locus} --out <out_dir>")


if __name__ == "__main__":
    main()