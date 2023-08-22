

## saveRDS(snakemake, snakemake@output[["rds"]])
## q()

## snakemake <- readRDS("/maps/projects/alab/people/rlk420/quilt2/human/HRC_NIPT/results/summary/quilt.nipt.accuracy.regular.refsize0.chr20.down1.0x.0.2f.rds")
## setwd("/maps/projects/alab/people/rlk420/quilt2/human/HRC_NIPT")

snakemake@source("common.R")
truth <- fread(snakemake@input[["truth"]], header=TRUE, data.table=FALSE)
af <- fread(snakemake@input[["af"]], header=FALSE, data.table=FALSE)
vcf<- fread(snakemake@input[["vcf"]], header=TRUE, data.table=FALSE)

id <- af[,1]
af <- as.numeric(af[,2])
names(af) <- id

vcf$ID <- paste(vcf[,1],vcf$POS, vcf$REF, vcf$ALT, sep=":")
vcfid <- vcf$ID
vcf <- vcf[,-seq(9)]

## mds <- sapply(sapply(dat[1,], strsplit, ":"),"[[", 3 )
## fds <- sapply(sapply(dat[1,], strsplit, ":"),"[[", 5 )

fds <- apply(vcf, 2, function(d) {
  ## as.numeric(sapply(sapply(d, strsplit, ":"),"[[", 5))
  as.numeric(sapply(str_split(d, ":"), "[[", 5))
})
rownames(fds) <- vcfid

mds <- apply(vcf, 2, function(d) {
  as.numeric(sapply(str_split(d, ":"), "[[", 3))
})
rownames(mds) <- vcfid

dat <- sapply(seq_along(nrow(truth)), function(i) {
  ## t <- sapply(truth[,i],strsplit,"|",fixed=T)
})

dat <- apply(truth, 1, function(r) {
  t <- str_split(r[-1], fixed("|"))
  as.numeric(unlist(t))
})
dat <- t(dat)
rownames(dat) <- truth[,1]
truth <- dat
rm(dat)
# save.image(file="nipt.RData")

bins <- sort(unique(c(
  c(0, 0.01 / 100, 0.02 / 100, 0.05 / 100),
  c(0, 0.01 / 10, 0.02 / 10, 0.05 / 10),
  c(0, 0.01 / 1, 0.02 / 1, 0.05 / 1),
  seq(0.1, 0.5, length.out = 5)
)))

acc_r2_by_af <- function(d0, d1, af, bins) {
  id <- intersect(intersect(rownames(d0), rownames(d1)), names(af))
  res <- r2_by_freq(bins, af, d0, d1, which_snps = id, flip = TRUE)
  as.data.frame(cbind(bin = bins[-1], r2 = res[, "simple"], n = res[, "n"], nA = res[, "nA"]))
}


## truth matal ds
dat <- truth[,seq(3, ncol(truth), 6)] + truth[,seq(4, ncol(truth), 6)]
acc.mat <- acc_r2_by_af(dat, mds, af, bins)

## truth fetal ds
dat <- truth[,seq(1, ncol(truth), 6)] + truth[,seq(2, ncol(truth), 6)]
acc.fet <- acc_r2_by_af(dat, fds, af, bins)

## 3 imputed haps
gts3 <- apply(vcf, 1, function(d) {
  o <- str_split(sapply(str_split(d, ":"), "[[", 1),
                 fixed("|"))
  as.numeric(unlist(o))
})
gts3 <- t(gts3)
rownames(gts3) <- vcfid


## truth paternal gts
## accuracy of paternal transimitted
ord <- sort(c(seq(5,ncol(truth),6), seq(6,ncol(truth),6)))
dat <- truth[,ord]
acc.pt1 <- acc_r2_by_af(dat[,seq(1,ncol(dat),2)], gts3[,seq(3,ncol(gts3),3)], af, bins)
acc.pt2 <- acc_r2_by_af(dat[,seq(2,ncol(dat),2)], gts3[,seq(3,ncol(gts3),3)], af, bins)

rds <- list(Mat=acc.mat, Fet=acc.fet, Pt1=acc.pt1, Pt2=acc.pt2)
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

pdf(paste0(snakemake@output[["rds"]], ".pdf"))

out <- as.data.frame(cbind(bin = acc.mat[, "bin"],
                           Mat = acc.mat[, "r2"],
                           Fet = acc.fet[, "r2"],
                           Pt1 = acc.pt1[, "r2"],
                           Pt2 = acc.pt2[, "r2"]))
groups <- colnames(out)[-1]
## fed(out)

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
     ylab = "Aggregated R2 within each MAF bin", xlab = "Minor Allele Frequency",
     main = title)
nd <- length(groups)
for (i in 1:nd) {
  y <- as.numeric(out[, groups[i]])
  points(x, y, pch = 17, col = mycols[i], cex = 1.3, type = "b")
}
axis(side = 1, at = x, labels = labels)
axis(side = 2)
legend("bottomright", legend = groups, fill=1:nd, bty = "n")

dev.off()
