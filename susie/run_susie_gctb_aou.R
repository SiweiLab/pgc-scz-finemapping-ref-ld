# AoU-perimeter minimal-dependency variant of run_susie_gctb.R.
#
# Why this exists: run_susie_gctb.R needs reticulate+numpy (to read the .npz LD
# matrix) plus locuszoomr/EnsDb.Hsapiens.v86/patchwork/ggrepel for the annotated
# locus plot. On a fresh AoU cluster only base R + a user-writable lib are
# available, and EnsDb.Hsapiens.v86 is a heavy Bioconductor annotation package
# not needed to confirm the PIP result. This script instead reads the
# space-delimited .z/.ld/.master triplet that finemap/build_finemap_input_gctb.py
# already produces from the SAME {locus}_ukb_R.npz + {locus}_sumstats.tsv.gz
# (so it is bit-for-bit the same z/R FINEMAP used), needing only susieR +
# ggplot2 (both easy to install/already present). Produces the same
# summary/credible-set/pip TSVs as run_susie_gctb.R, plus a simple (non-genomic-
# track) PIP scatter plot.
#
# Usage:
#   Rscript --vanilla run_susie_gctb_aou.R <analysis> <locus> <finemap_in_dir> \
#       <ukb_snps.txt.gz> <output_dir>
.libPaths(c(path.expand("~/R/library"), .libPaths()))
suppressMessages({
  library(susieR)
  library(ggplot2)
})

set.seed(1)

args           <- commandArgs(trailingOnly = TRUE)
analysis       <- args[1]
locus          <- args[2]
finemap_in_dir <- args[3]
ukb_snps_f     <- args[4]
output_dir     <- args[5]
# No defaults: a parameter with a default can silently disagree with run_params.py.
need <- function(i, name) {
  if (length(args) < i || is.na(args[i]) || args[i] == "") {
    stop(sprintf(paste0("missing required argument %d (%s). Pass every analysis parameter ",
                        "explicitly from aou_ld/run_params.py; this script has no defaults."),
                 i, name), call. = FALSE)
  }
  args[i]
}
L_arg          <- as.integer(need(6, "L"))
coverage       <- as.numeric(need(7, "coverage"))
min_abs_corr   <- as.numeric(need(8, "min_abs_corr"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

z_f      <- file.path(finemap_in_dir, paste0(locus, ".z"))
ld_f     <- file.path(finemap_in_dir, paste0(locus, ".ld"))
master_f <- file.path(finemap_in_dir, paste0(locus, ".master"))

z_df   <- read.table(z_f, header = TRUE, stringsAsFactors = FALSE)
R      <- as.matrix(read.table(ld_f, header = FALSE))
master <- read.table(master_f, header = TRUE, sep = ";", stringsAsFactors = FALSE)
N      <- master$n_samples[1]

snps <- read.table(gzfile(ukb_snps_f), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mn   <- setNames(snps$MarkerName, snps$SNP_ID)
MarkerName <- unname(mn[z_df$rsid])
MarkerName[is.na(MarkerName)] <- z_df$rsid[is.na(MarkerName)]

stopifnot(nrow(R) == nrow(z_df), ncol(R) == nrow(z_df))
z <- as.double(z_df$beta)   # se == 1 by construction, so beta == Zscore

cat("===", locus, "===\n")
cat("  SNPs:", length(z), "| N:", N, "\n")

off <- R[upper.tri(R)]
cat(sprintf("  R: min=%.3f  max=%.3f  frac_negative=%.1f%%\n",
            min(off), max(off), 100 * mean(off < 0)))
cat("  Min eigenvalue:", min(eigen(R, only.values = TRUE, symmetric = TRUE)$values), "\n")

lambda <- tryCatch(
  estimate_s_rss(z, R, n = N),
  error = function(e) { cat("  WARNING: estimate_s_rss failed:", conditionMessage(e), "\n"); NA }
)
cat(sprintf("  Lambda: %.4g\n", lambda))

diag_f <- file.path(output_dir, paste0(analysis, "_", locus, "_aou_kriging.png"))
tryCatch({
  condz <- kriging_rss(z, R, n = N)
  png(diag_f, width = 1200, height = 900, res = 120)
  print(condz$plot + coord_equal() +
          ggtitle(paste(analysis, locus, sprintf("| GCTB UKB (AoU) | lambda = %.4g", lambda))))
  dev.off()
  cat("  Saved kriging:", diag_f, "\n")
}, error = function(e) { cat("  WARNING: kriging_rss failed:", conditionMessage(e), "\n") })

fit <- tryCatch(
  susie_rss(z = z, R = R, n = N, L = L_arg, coverage = coverage, min_abs_corr = min_abs_corr),
  error = function(e) {
    if (grepl("prior variance is unreasonably large", conditionMessage(e))) {
      cat("  WARNING: prior variance error - retrying with check_prior=FALSE\n")
      tryCatch(
        susie_rss(z = z, R = R, n = N, L = L_arg, coverage = coverage,
                  min_abs_corr = min_abs_corr, check_prior = FALSE),
        error = function(e2) { cat("  ERROR on retry:", conditionMessage(e2), "\n"); NULL }
      )
    } else { cat("  ERROR in susie_rss:", conditionMessage(e), "\n"); NULL }
  }
)
if (is.null(fit)) stop("susie_rss failed")

n_cs      <- length(fit$sets$cs)
converged <- fit$converged
cat("  CS:", n_cs, "| Converged:", converged, "\n")

cs_names <- names(fit$sets$cs)
cs_label <- setNames(paste0("CS", seq_along(cs_names)), cs_names)
in_cs    <- rep(".", length(z))
if (n_cs > 0)
  for (cs_name in cs_names)
    in_cs[fit$sets$cs[[cs_name]]] <- cs_label[cs_name]

summary_f <- file.path(output_dir, paste0(analysis, "_susie_gctb_aou_summary.tsv"))
cs_f      <- file.path(output_dir, paste0(analysis, "_susie_gctb_aou_credible_sets.tsv"))
pip_f     <- file.path(output_dir, paste0(analysis, "_susie_gctb_aou_pip_all.tsv"))

write(paste(c("analysis", "locus", "n_snps", "N", "lambda", "n_cs", "converged"),
            collapse = "\t"), summary_f)
write(paste(c(analysis, locus, length(z), round(N), round(lambda, 6), n_cs, converged),
            collapse = "\t"), summary_f, append = TRUE)

write(paste(c("analysis", "locus", "cs", "MarkerName", "PIP"), collapse = "\t"), cs_f)
if (n_cs > 0)
  for (cs_name in cs_names) {
    idx <- fit$sets$cs[[cs_name]]
    for (j in seq_along(idx))
      write(paste(c(analysis, locus, cs_label[cs_name], MarkerName[idx[j]],
                    round(fit$pip[idx[j]], 6)), collapse = "\t"),
            cs_f, append = TRUE)
  }

pip_df <- data.frame(
  analysis = analysis, locus = locus,
  MarkerName = MarkerName, SNP_ID = z_df$rsid, position = z_df$position,
  Zscore = round(z, 6), PIP = round(fit$pip, 6), in_cs = in_cs
)
write.table(pip_df, pip_f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

plot_f <- file.path(output_dir, paste0(analysis, "_", locus, "_aou_pip.png"))
p <- ggplot(pip_df, aes(x = position, y = PIP, colour = in_cs)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = 0.1, linetype = "dashed", colour = "grey50") +
  theme_bw() +
  ggtitle(sprintf("%s %s | GCTB UKB (AoU) | lambda=%.4g | CS=%d", analysis, locus, lambda, n_cs)) +
  labs(x = "position (hg19)", y = "PIP (SuSiE-RSS)")
png(plot_f, width = 2000, height = 1200, res = 200)
print(p)
dev.off()

cat("  Saved:", plot_f, "\n")
cat("Done. Outputs:\n", summary_f, "\n", cs_f, "\n", pip_f, "\n")
