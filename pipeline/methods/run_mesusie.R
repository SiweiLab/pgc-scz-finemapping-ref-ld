#!/usr/bin/env Rscript
# Run MESuSiE on a harmonized locus.
#
# Reads the per-ancestry {locus}.{anc}.ma + {locus}.{anc}.ld.gz files produced by
# ../harmonize.py. MESuSiE classifies each signal as shared vs ancestry-specific.
#
# Install: library(devtools); install_github("borangao/MESuSiE")
# Tutorial: https://borangao.github.io/meSuSie_Analysis/
#
# Outputs (in --out):
#   {locus}_mesusie_pip.tsv         rsID + PIP
#   {locus}_mesusie_config.tsv      per-CS shared/ancestry-specific classification
#
# Usage:
#   Rscript run_mesusie.R chr1_5495783_8495783 EUR,AFR,EAS \
#       /scratch/.../ma_harmonized /scratch/.../mesusie_out [coverage] [L]
#
# L defaults to 10 (matching every previously published run of this wrapper,
# which had no L argument at all -- meSuSie_core(..., L = 10) was hardcoded
# until 2026-09-08's L-sweep pilot needed it exposed).


args <- commandArgs(trailingOnly = TRUE)
locus      <- args[1]

# No defaults: a parameter with a default can silently disagree with run_params.py, which
# is how the manuscript's SuSiE-merged ended up at min_abs_corr 0.1.
need <- function(i, name) {
  if (length(args) < i || is.na(args[i]) || args[i] == "") {
    stop(sprintf(paste0("missing required argument %d (%s). Pass every analysis parameter ",
                        "explicitly from pipeline/run_params.py; this script has no defaults."),
                 i, name), call. = FALSE)
  }
  args[i]
}
ancestries <- strsplit(args[2], ",")[[1]]
harm_dir   <- args[3]
out_dir    <- args[4]
coverage   <- as.numeric(need(5, "coverage"))
L_arg      <- as.integer(need(6, "L"))

# Packages load only after every required argument is present, so a missing parameter
# fails on the guard even where MESuSiE is not installed.
suppressMessages({ library(MESuSiE); library(data.table) })
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Fixed re-implementation of MESuSiE::meSuSie_get_cs ──────────────────────
# The installed package's meSuSie_get_cs(res, Xcorr, coverage=..., ...) accepts
# a `coverage` argument but never forwards it: its body hardcodes
# `status = in_CS(alpha_either, coverage = 0.95)` regardless of the caller's
# value (confirmed by printing the installed function body directly) -- a
# real upstream bug, not a usage error. `in_CS()` itself IS exported and DOES
# take a working `coverage` arg, so this re-implementation is just the
# original function body with that one call fixed to use the real argument.
get_cs_fixed <- function(res, Xcorr, coverage = 0.95, prior_tol = 1e-9,
                         cor_method = "min.abs.corr", cor_threshold = 0.5) {
  include_idx <- unlist(lapply(res$V, function(x) max(diag(x)) > prior_tol))
  alpha_either <- t(Reduce(cbind, lapply(res$alpha, function(x) apply(x, 1, sum))))
  status <- MESuSiE::in_CS(alpha_either, coverage = coverage)
  cs <- lapply(1:nrow(status), function(i) which(status[i, ] != 0))
  claimed_coverage <- sapply(1:length(cs), function(i) sum(alpha_either[i, ][cs[[i]]]))
  include_idx <- include_idx * (lapply(cs, length) > 0)
  include_idx <- include_idx * (!(duplicated(cs)))
  include_idx <- as.logical(include_idx)
  if (sum(include_idx) == 0) {
    return(list(cs = NULL, cs_category = NULL, purity = NULL,
                cs_index = NULL, coverage = NULL, requested_coverage = coverage))
  }
  cs <- cs[include_idx]; claimed_coverage <- claimed_coverage[include_idx]
  purity <- data.frame(do.call(rbind, lapply(1:length(cs), function(i) {
    MESuSiE::get_purity(cs[[i]], Xcorr)
  })))
  colnames(purity) <- c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
  is_pure <- which(purity[, colnames(purity) %in% cor_method] >= cor_threshold)
  if (length(is_pure) == 0) {
    return(list(cs = NULL, cs_category = NULL, purity = NULL,
                cs_index = NULL, coverage = NULL, requested_coverage = coverage))
  }
  cs <- cs[is_pure]; purity <- purity[is_pure, ]; claimed_coverage <- claimed_coverage[is_pure]
  row_names <- paste0("L", which(include_idx)[is_pure])
  names(cs) <- row_names; rownames(purity) <- row_names
  cs_index <- which(include_idx)[is_pure]
  cs_category <- unlist(lapply(res$alpha[cs_index], function(x) res$name_config[which.max(colSums(x))]))
  names(cs_category) <- row_names
  list(cs = cs, cs_category = cs_category, purity = purity,
       cs_index = cs_index, coverage = claimed_coverage, requested_coverage = coverage)
}

# ── Load per-ancestry summary stats + LD ─────────────────────────────────────
# MESuSiE expects a list of summary-stat data.frames and a list of LD matrices,
# all on the SAME, identically-ordered SNP set (guaranteed by harmonize_loci.py).
summ_list <- list(); R_list <- list()
for (anc in ancestries) {
  ma <- fread(file.path(harm_dir, paste0(locus, ".", anc, ".ma")))
  # MESuSiE summary format: SNP, Z (and N); derive Z from b/se.
  ma[, Z := b / se]
  summ_list[[anc]] <- data.frame(SNP = ma$SNP, Z = ma$Z, N = ma$N,
                                 Beta = ma$b, Se = ma$se)
  R_mat <- as.matrix(fread(file.path(harm_dir,
                           paste0(locus, ".", anc, ".ld.gz")), header = FALSE))
  # meSuSie_core's own bundled example (data(LD_list)) carries SNP-name row/col
  # names matching summary_stat_list's SNP order exactly -- our .ld.gz has none
  # (plain fread(header=FALSE)); set them explicitly to match that convention.
  dimnames(R_mat) <- list(ma$SNP, ma$SNP)
  R_list[[anc]] <- R_mat
}
names(summ_list) <- ancestries
names(R_list)    <- ancestries

# ── Run MESuSiE ──────────────────────────────────────────────────────────────
# Verified against the installed package (2026-07-07): meSuSie_core(R_mat_list,
# summary_stat_list, L, ...) and summary_stat_list columns SNP/Beta/Se/Z/N all
# match ?meSuSie_core exactly -- this call itself was already correct.
fit <- MESuSiE::meSuSie_core(R_list, summ_list, L = L_arg)

# ── Write outputs ────────────────────────────────────────────────────────────
pip <- fit$pip
snps <- summ_list[[1]]$SNP
fwrite(data.frame(rsid = snps, pip = pip),
       file.path(out_dir, paste0(locus, "_mesusie_pip.tsv")), sep = "\t")

# Credible sets with shared/ancestry-specific config label, re-extracted at the
# REQUESTED coverage via get_cs_fixed() (fit$cs itself is always the package's
# hardcoded-0.95 extraction -- see get_cs_fixed's comment above). `fit` here is
# the meSuSieObject_obj itself (meSuSie_core returns it after mutating it in
# place), so it already carries the $V/$alpha/$name_config fields get_cs_fixed
# needs -- no need to refit.
cs <- get_cs_fixed(fit, R_list, coverage = coverage)
if (!is.null(cs) && length(cs$cs) > 0) {
  rows <- do.call(rbind, lapply(names(cs$cs), function(cs_name) {
    idx <- cs$cs[[cs_name]]
    data.frame(cs_id = cs_name, rsid = snps[idx],
               config = cs$cs_category[[cs_name]])
  }))
  fwrite(rows, file.path(out_dir, paste0(locus, "_mesusie_config.tsv")), sep = "\t")
  cat(sprintf("MESuSiE done: %d credible sets, %d SNPs PIP>0.9\n",
              length(cs$cs), sum(pip > 0.9, na.rm = TRUE)))
} else {
  cat(sprintf("MESuSiE done: 0 credible sets, %d SNPs PIP>0.9\n",
              sum(pip > 0.9, na.rm = TRUE)))
}
