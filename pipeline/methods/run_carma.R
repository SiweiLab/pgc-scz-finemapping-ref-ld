#!/usr/bin/env Rscript
# CARMA pilot on one merged-LD locus, at CARMA's own installed defaults.
# Merged-LD arm, so z.list and ld.list both have length 1. Inputs are FINEMAP-format
# finemap_in_rebuilt/{locus}.z (se = 1, so beta is the Z score) and {locus}.ld.
#
# Usage: Rscript run_carma.R <locus> <in_dir> <out_dir> <rho.index> <num.causal> <BF.index> [<seed>]
#   The three parameters are REQUIRED; this script has no defaults.
# Run under the conda env at ~/envs/carma (R 4.5.3, CARMA 1.0 @ 856f0fb).
suppressPackageStartupMessages(library(CARMA))
args <- commandArgs(trailingOnly = TRUE)

# Every analysis parameter must be passed explicitly. A parameter with a default is a
# parameter that can silently disagree with run_params.py -- which is exactly how the
# manuscript's SuSiE-merged ended up at min_abs_corr 0.1.
need <- function(i, name) {
  if (length(args) < i || is.na(args[i]) || args[i] == "") {
    stop(sprintf(paste0("missing required argument %d (%s). Pass every analysis parameter ",
                        "explicitly from pipeline/run_params.py; this script has no defaults."),
                 i, name), call. = FALSE)
  }
  args[i]
}
locus <- args[1]; in_dir <- args[2]; out_dir <- args[3]
rho_index   <- as.numeric(need(4, "rho.index (CARMA's coverage)"))
num_causal  <- as.integer(need(5, "num.causal"))
bf_index    <- as.numeric(need(6, "BF.index"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

z <- read.table(file.path(in_dir, paste0(locus, ".z")), header = TRUE, stringsAsFactors = FALSE)
ld <- as.matrix(read.table(file.path(in_dir, paste0(locus, ".ld"))))
stopifnot(nrow(ld) == nrow(z), ncol(ld) == nrow(z), all(z$se == 1))
zscore <- z$beta / z$se

cat("CARMA", as.character(packageVersion("CARMA")), "| locus", locus, "| n_snps", nrow(z), "\n")
# output.labels is where CARMA writes its own side files; point it at out_dir.
# lambda.list is NOT optional despite formals() defaulting it to NULL: the installed body
# calls log(lambda.list[[i]]) unconditionally ("non-numeric argument to mathematical
# function"), and the man page gives no default. 1/sqrt(p) is the value the package's own
# example passes -- the only documented one. Every other argument is left at formals().
lambda <- 1 / sqrt(nrow(z))
cat(sprintf("lambda (1/sqrt(p)) %.6f\n", lambda))
t0 <- proc.time()[["elapsed"]]
# rho.index is CARMA's coverage (paper: "probability rho or greater of containing at least
# one causal variant"). Passed explicitly even though its default already equals ours -- a
# default that silently matches is indistinguishable from one that silently does not.
# Argument 7, the RNG seed (run_params.SEED_CARMA). CARMA's model search is stochastic; without a seed, identical
# runs differ (chr5_86721950: max |dPIP| 0.023) and so does their run time. It is optional ONLY so that a driver
# already running with the old 6-argument argv keeps working; that case is logged as unseeded.
seed <- if (length(args) >= 7 && !is.na(args[7]) && args[7] != "") as.integer(args[7]) else NA_integer_
if (is.na(seed)) cat("WARNING: no seed (argument 7); CARMA's stochastic search is UNSEEDED\n") else set.seed(seed)
cat(sprintf("seed %s\n", if (is.na(seed)) "none" else seed))
res <- CARMA(z.list = list(zscore), ld.list = list(ld), lambda.list = list(lambda),
             output.labels = out_dir,
             rho.index = rho_index, num.causal = num_causal, BF.index = bf_index,
             outlier.switch = TRUE)
elapsed <- proc.time()[["elapsed"]] - t0
cat(sprintf("elapsed_sec %.1f\n", elapsed))

r <- res[[1]]
cat("result components:", paste(names(r), collapse = ", "), "\n")
saveRDS(res, file.path(out_dir, paste0(locus, "_carma.rds")))

pip <- data.frame(rsid = z$rsid, bp = z$position, z = zscore, pip = r$PIPs)
write.table(pip, file.path(out_dir, paste0(locus, "_carma_pip.tsv")), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("top PIPs:\n"); print(head(pip[order(-pip$pip), ], 5), row.names = FALSE)

cs <- r$`Credible set`
cat("Credible set component:\n"); str(cs, max.level = 2)
if (!is.null(cs[[2]]) && length(cs[[2]]) > 0) {
  for (i in seq_along(cs[[2]])) {
    idx <- cs[[2]][[i]]
    cat(sprintf("  CS%d: %d variants; top %s (pip %.4f)\n", i, length(idx),
                z$rsid[idx][which.max(r$PIPs[idx])], max(r$PIPs[idx])))
  }
}
out <- r$Outliers
cat("Outliers component:\n"); str(out)
if (is.data.frame(out) && nrow(out) > 0) {
  idx <- out$Index
  write.table(data.frame(out, rsid = z$rsid[idx], bp = z$position[idx], z = zscore[idx]),
              file.path(out_dir, paste0(locus, "_carma_outliers.tsv")), sep = "\t",
              row.names = FALSE, quote = FALSE)
}
