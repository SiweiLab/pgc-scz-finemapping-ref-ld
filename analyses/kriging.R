#!/usr/bin/env Rscript
# susieR's LD-consistency diagnostics (estimate_s_rss + kriging_rss) on the exact SuSiEx inputs of the v2 primary run
# (no-palindromic): per ancestry, EUR z against GCTB-10k EUR LD and EAS z against AoU EAS LD. Diagnostic only, so nothing is
# removed or refitted.
#
# Usage: Rscript kriging.R <locus> <susiex_dir> <out_dir> <n_eur> <n_eas>
#   susiex_dir = ~/work/v2_harmonized/noambig/union/<locus>_susiex   ({EUR,EAS}.sst.txt, {EUR,EAS}.ld.bin)
# Writes, per ancestry:
#   <out>/<locus>_<ANC>_kriging.tsv  snp, bp, z, condmean, condvar, z_std_diff, logLR, outlier (susieR's rule: logLR > 2 and |z| > 2)
#   <out>/<locus>_<ANC>_kriging.png  susieR's own kriging_rss plot
#   <out>/<locus>_<ANC>_s.txt        s from estimate_s_rss (null-mle), n, P
suppressMessages(library(susieR))
a <- commandArgs(trailingOnly = TRUE)
locus <- a[1]; d <- a[2]; out <- a[3]
n_by <- c(EUR = as.numeric(a[4]), EAS = as.numeric(a[5]))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
for (anc in c("EUR", "EAS")) {
  sst <- read.table(file.path(d, paste0(anc, ".sst.txt")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  P <- nrow(sst)
  R <- matrix(readBin(file.path(d, paste0(anc, ".ld.bin")), "numeric", size = 4, n = P * P), P, P, byrow = TRUE)
  z <- sst$stat
  n <- n_by[[anc]]
  s <- estimate_s_rss(z, R, n = n)
  k <- kriging_rss(z, R, n = n, s = s)
  cd <- k$conditional_dist
  cd$outlier <- cd$logLR > 2 & abs(cd$z) > 2
  write.table(cbind(snp = sst$snp, bp = sst$bp, cd), file.path(out, paste0(locus, "_", anc, "_kriging.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  ggplot2::ggsave(file.path(out, paste0(locus, "_", anc, "_kriging.png")),
                  k$plot + ggplot2::ggtitle(sprintf("%s %s  s = %.4g", locus, anc, s)), width = 5, height = 5, dpi = 150)
  writeLines(sprintf("locus\tanc\ts\tn\tP\tn_outlier\n%s\t%s\t%.6g\t%d\t%d\t%d", locus, anc, s, as.integer(n), P,
                     sum(cd$outlier)), file.path(out, paste0(locus, "_", anc, "_s.txt")))
  cat(locus, anc, "P", P, "s", signif(s, 4), "outliers", sum(cd$outlier), "\n")
}
