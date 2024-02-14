snakemake@source("common.R")

## saveRDS(snakemake, snakemake@output[["rds"]])
## q()
## setwd("/maps/projects/alab/people/rlk420/quilt2/human/HRC_NIPT")
## snakemake <- readRDS("results/summary/quilt.nipt.accuracy.regular.refsize0.chr20.down1.0x.0.0f.rds")

r2persam <- function(imputed, truth, samples, sites){
  ord <- match(imputed$samples, samples$FamilyID)
  mid <- samples[ord,"Mid"]
  ord2 <- match(mid, truth$samples)
  gt <- truth[[10]]
  ds <- imputed[[10]]
  true <- gt[match(sites, truth$id), ord2]
  test <- ds[match(sites, imputed$id), ]
  r2mid <- diag(cor(true, test)**2)
  list("mother"=r2mid, "kid"=rep(NA, length(r2mid)))
}

r2perfam <- function(imputed, truth, samples, sites){
  ord <- match(imputed$samples, samples$FamilyID)
  mid <- samples[ord,"Mid"]
  ord2 <- match(mid, truth$samples)
  gt <- truth[[10]]
  ds <- imputed$mds
  true <- gt[match(sites, truth$id), ord2]
  test <- ds[match(sites, imputed$id), ]
  r2mid <- diag(cor(true, test)**2)
  ## get kid genotype dosage
  kid <- samples[ord,"Kid"]
  ord2 <- match(kid, truth$samples)
  ds <- imputed$fds
  true <- gt[match(sites, truth$id), ord2]
  test <- ds[match(sites, imputed$id), ]
  r2kid <- diag(cor(true, test)**2)
  list("mother"=r2mid, "kid"=r2kid)
}

r2famperaf <- function(imputed, truth, af, bins, sites) {
  ord <- match(imputed$samples, samples$FamilyID)
  mid <- samples[ord,"Mid"]
  ord2 <- match(mid, truth$samples)
  gt <- truth[[10]]
  ds <- imputed$mds
  true <- gt[match(sites, truth$id), ord2]
  test <- ds[match(sites, imputed$id), ]
  rownames(true) <- sites
  rownames(test) <- sites
  res <- r2_by_freq(bins, af, true, test, which_snps = sites, flip = FALSE)
  kid <- samples[ord,"Kid"]
  ord2 <- match(kid, truth$samples)
  ds <- imputed$fds
  true <- gt[match(sites, truth$id), ord2]
  test <- ds[match(sites, imputed$id), ]
  rownames(true) <- sites
  rownames(test) <- sites
  res2 <- r2_by_freq(bins, af, true, test, which_snps = sites, flip = FALSE)
  as.data.frame(cbind(bin = bins[-1], r2mds = res[, "simple"], r2fds = res2[,"simple"], n = res[, "n"], nA = res[, "nA"]))
}

r2samperaf <- function(imputed, truth, af, bins, sites) {
  ord <- match(imputed$samples, samples$FamilyID)
  mid <- samples[ord,"Mid"]
  ord2 <- match(mid, truth$samples)
  gt <- truth[[10]]
  ds <- imputed[[10]]
  true <- gt[match(sites, truth$id), ord2]
  test <- ds[match(sites, imputed$id), ]
  rownames(true) <- sites
  rownames(test) <- sites
  res <- r2_by_freq(bins, af, true, test, which_snps = sites, flip = FALSE)
  as.data.frame(cbind(bin = bins[-1], r2mds = res[, "simple"], r2fds = rep(NA, nrow(res)), n = res[, "n"], nA = res[, "nA"]))
}

ff <- as.numeric(snakemake@wildcards[["ff"]])
samplefile <- snakemake@config[["samples_nipt"]]
samples <- read.table(samplefile,h=T)
niptvcf <- snakemake@input[["vcf"]]
imputed <- vcftable(niptvcf, collapse = T, format = ifelse(ff==0, "DS", "GT"))
imputed$id <- paste(imputed$chr, imputed$pos, imputed$ref, imputed$alt, sep = ":")
chr <- unique(imputed$chr)
region <- paste0(chr,":", paste0(range(imputed$pos), collapse = "-"))
truthvcf <- snakemake@params[["truth"]]
truth <- vcftable(truthvcf, region = region)
truth$id <- paste(truth$chr, truth$pos, truth$ref, truth$alt, sep = ":")

af <- fread(snakemake@input[["af"]], header=FALSE, data.table=FALSE)
sites <-  intersect(imputed$id, truth$id)
sites <-  intersect(af[,1], sites)
af <- af[match(sites, af[,1]), 2]
names(af) <- sites

bins <- sort(unique(c(
  ## c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
  c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
  c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
  seq(0.1, 0.5, length.out = 5)
)))

if(ff > 0) {
  TMP <- vcftable(niptvcf, format = "MDS")
  imputed[["mds"]] <- TMP[[10]]
  TMP <- vcftable(niptvcf, format = "FDS")
  imputed[["fds"]] <- TMP[[10]]
  r2samples <- r2perfam(imputed, truth, samples, sites)
  r2bins <- r2famperaf(imputed, truth, af, bins, sites)
} else {
  r2samples <- r2persam(imputed, truth, samples, sites)
  r2bins <- r2samperaf(imputed, truth, af, bins, sites)
}

rds <- list(r2bins, r2samples)

saveRDS(rds,snakemake@output[["rds"]])
## q()

## source("workflow/scripts/common.R")

## r <- readRDS("nipt.rds")
## acc.mat <- r[[1]]
## acc.fet <- r[[2]]
## truth <- r[[3]]
## mds <- r[[4]]
## fds <- r[[5]]
## af <- r[[6]]

mycols <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")
palette(mycols)

## pdf(paste0(snakemake@output[["rds"]], ".pdf"), w=12, h=6)
png(paste0(snakemake@output[["rds"]], ".png"), w=12, h=6, res = 300, units = "in")

par(mfrow=c(1,2))

r2bins <- rds[[1]]
out <- as.data.frame(cbind(bin = r2bins[, "bin"],
                           Mat = r2bins[, "r2mds"],
                           Fet = r2bins[, "r2fds"]))
groups <- colnames(out)[-1]

x <- out$bin[!sapply(out$bin, is.na)] # remove AF bin with NA results
x <- log10(as.numeric(x))
labels <- 100 * out$bin
labels <- labels[!sapply(out$bin, is.na)]

title <- paste(snakemake@wildcards[["chrom"]],
               ", refsize:",
               snakemake@wildcards[["size"]],
               ", depth:",
               snakemake@wildcards[["depth"]],
               ", ff:",
               snakemake@wildcards[["ff"]] )

plot(1, col = "transparent", axes = FALSE,
     xlim = c(min(x), max(x)), ylim = c(0, 1.0),
     ylab = "Aggregated R2 within each MAF bin", xlab = "Allele Frequency (%)",
     main = title)
nd <- length(groups)
for (i in 1:nd) {
  y <- as.numeric(out[, groups[i]])
  points(x, y, pch = 17, col = mycols[i], cex = 1.3, type = "b")
}
axis(side = 1, at = x, labels = labels)
axis(side = 2)
legend("bottomright", legend = groups, fill=1:nd, bty = "n")

boxplot(r2samples, ylab = "r2 per sample", main = paste("#sites=",length(sites)), col=1:2, cex.axis=2, cex.lab = 2, cex.main = 2)

dev.off()

q()
