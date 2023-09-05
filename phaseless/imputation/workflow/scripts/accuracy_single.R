## saveRDS(snakemake, snakemake@output[["rds"]])
## q()


## snakemake <- readRDS("/maps/projects/alab/people/rlk420/phaseless/imputation/ceu/results/stitch/accuracy.down2.0x.chr20.rds")
## setwd("/maps/projects/alab/people/rlk420/phaseless/imputation/ceu/")

snakemake@source("common.R")

samples <- as.character(snakemake@params[["samples"]])

af <- fread(snakemake@input[["af"]], header=FALSE, data.table=FALSE)

id <- sapply(strsplit(af[,1], ":"), function(l) paste(l[1], l[2], sep=":"))
af <- as.numeric(af[,2])
names(af) <- id

## dat <- apply(truth[,-1], 2, function(a) {
##   o <- sapply(str_split(a, fixed("|")), as.numeric)
##   colSums(o)
## } )
## rownames(dat) <- truth[,1]
## truth <- dat
## dim(truth)
## rm(dat)

library(vcfppR)
vcffile <- snakemake@input[["vcf"]]
truthfile <- snakemake@input[["truth"]]

## truthfile <- "/maps/projects/alab/people/rlk420/truth/CEU/everyone.UKBB_GEL.chr20.ligate.phased.vcf.gz"
## Sys.setenv(PATH="/usr/bin:/bin:/sbin:/usr/sbin:/usr/local/bin:/usr/local/sbin:/home/rlk420/.local/bin")
## system(paste("bcftools", "index", "-f", vcffile))

vcf <- tableDS(vcffile, snakemake@wildcards[["chrom"]], samples)
vcfid <- paste(vcf$chr,vcf$pos, sep=":")
ds <- vcf[["ds"]]
names(ds) <- vcfid

truth <- tableGT(truthfile, snakemake@wildcards[["chrom"]], samples)
truthid <- paste(truth$chr,truth$pos, sep=":")
gt <- truth[["gt"]]
names(gt) <- truthid
samples <- unlist(strsplit(samples, ","))
ord <- match(samples, truth[["samples"]])  ## fix order
gt <- lapply(gt, function(g2) {
  colSums(matrix(g2, nrow=2))
})

id <- intersect(intersect(names(gt), names(ds)), names(af))

bins <- sort(unique(c(
  c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
  c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
  c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
  seq(0.1, 0.5, length.out = 5)
)))


acc_r2_by_af <- function(gt, ds, af, id, bins) {
  truthG <- do.call(rbind, gt[id])
  truthG <- truthG[,ord]
  testDS <- do.call(rbind, ds[id])
  res <- r2_by_freq(bins, af, truthG, testDS, which_snps = id, flip = TRUE)
  as.data.frame(cbind(bin = bins[-1], r2 = res[, "simple"], n = res[, "n"], nA = res[, "nA"]))
}


out <- acc_r2_by_af(gt, ds, af, id, bins)
saveRDS(out, snakemake@output[["rds"]])

