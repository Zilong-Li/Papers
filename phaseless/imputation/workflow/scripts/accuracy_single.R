saveRDS(snakemake, snakemake@output[["rds"]])
q()

snakemake@source("common.R")

acc_r2_by_af <- function(d0, d1, af, bins) {
  id <- intersect(intersect(rownames(d0), rownames(d1)), names(af))
  res <- r2_by_freq(bins, af, d0, d1, which_snps = id, flip = TRUE)
  as.data.frame(cbind(bin = bins[-1], single = res[, "simple"], orphan = res[, "simple"]))
}
