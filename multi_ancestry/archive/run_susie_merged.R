#!/usr/bin/env Rscript
# Run SuSiE-RSS (susieR::susie_rss, matching susie/run_susie_gctb.R's recipe:
# L=10, coverage=0.95, min_abs_corr=0.1) plus the lambda/kriging diagnostic
# (estimate_s_rss/kriging_rss), directly on FINEMAP's own .z/.ld inputs for
# the 9-loci EUR+EAS merged-LD rerun -- same z/R FINEMAP already used
# (finemap/build_finemap_input_gctb.py's output), so this is an exact
# apples-to-apples second method on identical data, no separate LD-loading
# path needed.
#
# Skips run_susie_gctb.R's locuszoomr/EnsDb.Hsapiens.v86/patchwork/reticulate
# dependencies (not installed on this cluster, no "carma" conda env either)
# -- just susieR + ggplot2 (both already available). Plotting is handled
# separately in Python (plot_locus_finemap_merged.py) to match its existing
# panel styling; this script only writes the kriging diagnostic PNG (which
# comes straight from susieR's own kriging_rss()$plot, not locuszoomr) plus
# tsv summary/PIP tables.
#
# Usage: Rscript run_susie_merged.R <locus> <finemap_dir> <out_dir> <coverage> <min_abs_corr> <L>
# All three analysis parameters are required (B30 coverage/min_abs_corr, B31 L); the defaults
# described below are historical and no longer apply.
# coverage defaults to 0.95 (susie_rss's own default); pass 0.99 to match the
# published PGC3-SCZ SuSiEx credible-set coverage for an apples-to-apples
# comparison against ST11's multi-ancestry results.
# min_abs_corr defaults to 0.1 (this script's original value) for backward
# compatibility with existing outputs -- but the per-ancestry methods in this
# same comparison (MultiSuSiE --min-abs-corr, MESuSiE cor_threshold, SuSiEx
# --min_purity) all use 0.5. Pass 0.5 explicitly to make the merged-LD negative
# control an apples-to-apples comparison on purity filtering, not just coverage.
suppressMessages({
  library(susieR)
  library(ggplot2)
})

args         <- commandArgs(trailingOnly = TRUE)

# No defaults: a parameter with a default can silently disagree with run_params.py, which
# is how the manuscript's SuSiE-merged ended up at min_abs_corr 0.1 (ORIGINAL_RUN_SPEC 3.6).
need <- function(i, name) {
  if (length(args) < i || is.na(args[i]) || args[i] == "") {
    stop(sprintf(paste0("missing required argument %d (%s). Pass every analysis parameter ",
                        "explicitly from aou_ld/run_params.py; this script has no defaults."),
                 i, name), call. = FALSE)
  }
  args[i]
}
locus        <- args[1]
finemap_dir  <- args[2]
out_dir      <- args[3]
coverage     <- as.numeric(need(4, "coverage"))
min_abs_corr <- as.numeric(need(5, "min_abs_corr"))
L            <- as.integer(need(6, "L"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

zf       <- file.path(finemap_dir, paste0(locus, ".z"))
ldf      <- file.path(finemap_dir, paste0(locus, ".ld"))
master_f <- file.path(finemap_dir, paste0(locus, ".master"))

z_df   <- read.table(zf, header = TRUE, stringsAsFactors = FALSE)
R      <- as.matrix(read.table(ldf))
master <- read.table(master_f, header = TRUE, sep = ";", stringsAsFactors = FALSE)
N      <- master$n_samples[1]

z <- z_df$beta  # se=1 throughout (build_finemap_input_gctb.py's convention), so beta==z
stopifnot(nrow(R) == length(z))

cat("=== ", locus, " ===\n", sep = "")
cat("  N SNPs:", length(z), "| n_samples:", N, "\n")

lambda <- tryCatch(
  estimate_s_rss(z, R, n = N),
  error = function(e) { cat("  WARNING estimate_s_rss:", conditionMessage(e), "\n"); NA }
)
cat(sprintf("  lambda: %.6g\n", lambda))

diag_f <- file.path(out_dir, paste0(locus, "_kriging.png"))
tryCatch({
  condz <- kriging_rss(z, R, n = N)
  png(diag_f, width = 1200, height = 900, res = 120)
  print(condz$plot + coord_equal() +
          ggtitle(paste(locus, sprintf("| EUR+EAS merged LD | lambda = %.4g", lambda))))
  dev.off()
  cat("  saved kriging:", diag_f, "\n")
}, error = function(e) cat("  WARNING kriging_rss:", conditionMessage(e), "\n"))

fit <- tryCatch(
  susie_rss(z = z, R = R, n = N, L = L, coverage = coverage, min_abs_corr = min_abs_corr),
  error = function(e) {
    if (grepl("prior variance is unreasonably large", conditionMessage(e))) {
      cat("  WARNING: prior variance error -- retrying with check_prior=FALSE\n")
      tryCatch(
        susie_rss(z = z, R = R, n = N, L = L, coverage = coverage,
                  min_abs_corr = min_abs_corr, check_prior = FALSE),
        error = function(e2) { cat("  ERROR retry:", conditionMessage(e2), "\n"); NULL }
      )
    } else { cat("  ERROR susie_rss:", conditionMessage(e), "\n"); NULL }
  }
)
if (is.null(fit)) { cat("  ABORT: no SuSiE fit\n\n"); quit(status = 1) }

n_cs      <- length(fit$sets$cs)
converged <- fit$converged
cat("  CS:", n_cs, "| converged:", converged, "\n")

cs_names <- names(fit$sets$cs)
cs_label <- setNames(paste0("CS", seq_along(cs_names)), cs_names)
in_cs    <- rep(".", length(z))
if (n_cs > 0)
  for (cs_name in cs_names) in_cs[fit$sets$cs[[cs_name]]] <- cs_label[cs_name]

summary_df <- data.frame(locus = locus, n_snps = length(z), N = round(N),
                          lambda = round(lambda, 6), coverage = coverage, L = L,
                          n_cs = n_cs, converged = converged)
write.table(summary_df, file.path(out_dir, paste0(locus, "_susie_summary.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

pip_df <- data.frame(locus = locus, rsid = z_df$rsid, chromosome = z_df$chromosome,
                      position = z_df$position, Zscore = round(z, 6),
                      PIP = round(fit$pip, 6), in_cs = in_cs)
write.table(pip_df, file.path(out_dir, paste0(locus, "_susie_pip.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("  saved:", file.path(out_dir, paste0(locus, "_susie_summary.tsv")), "\n")
cat("  saved:", file.path(out_dir, paste0(locus, "_susie_pip.tsv")), "\n\n")
