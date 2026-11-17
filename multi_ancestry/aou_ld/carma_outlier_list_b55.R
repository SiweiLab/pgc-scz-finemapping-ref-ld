# B55: list CARMA's outlier variants (result$Outliers$Index, 1-based into the .z rows) for every locus of a run.
# Usage: Rscript carma_outlier_list_b55.R <carma_out_dir> <merged_in_dir> <out.tsv>
a <- commandArgs(trailingOnly = TRUE)
fs <- list.files(a[1], pattern = "_carma\\.rds$", full.names = TRUE)
rows <- list()
for (f in fs) {
  L <- sub("_carma\\.rds$", "", basename(f))
  r <- readRDS(f)[[1]]
  o <- r$Outliers
  z <- read.table(file.path(a[2], paste0(L, ".z")), header = TRUE, stringsAsFactors = FALSE)
  if (!is.null(o) && nrow(o) > 0) rows[[L]] <- data.frame(locus = L, index = o$Index, rsid = z$rsid[o$Index], z = o$Z)
  else rows[[L]] <- data.frame(locus = L, index = NA, rsid = NA, z = NA)
}
d <- do.call(rbind, rows)
write.table(d, a[3], sep = "\t", row.names = FALSE, quote = FALSE)
cat("loci", length(fs), "outliers", sum(!is.na(d$index)), "\n")
