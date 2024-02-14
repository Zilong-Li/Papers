
## setwd("/maps/projects/alab/people/rlk420/quilt2/human/HRC_NIPT")
## saveRDS(snakemake, snakemake@output[["rds"]])
## q()

## snakemake <- readRDS("/maps/projects/alab/people/rlk420/quilt2/human/HRC_NIPT/results/summary/quilt.nipt.accuracy.regular.refsize0.chr20.rds")
## source("workflow/scripts/common.R")

snakemake@source("common.R")
depths <- as.character(unlist(snakemake@config["coverage"]))
fetalfrac <- as.character(unlist(snakemake@config["fetalfrac"]))

input <- lapply(snakemake@input, readRDS)

lsp <- list()
i <- 1
for(dep in seq_along(depths)){
  lff <- list()
  for(ff in seq_along(fetalfrac)){
    rds <- input[[i]]
    out <- as.data.frame(cbind(bin = rds[[1]][,"bin"],
                               Mat = rds[[1]][,"r2mds"],
                               Fet = rds[[1]][,"r2fds"]))
    lff[[ff]] <- out
    i <- i + 1
  }
  names(lff) <- fetalfrac
  lsp[[dep]] <- lff
}
names(lsp) <- depths

dsp <- list()
i <- 1
for(dep in seq_along(depths)){
  lff <- list()
  for(ff in seq_along(fetalfrac)){
    rds <- input[[i]]
    lff[[ff]] <- do.call(cbind, rds[[2]])
    i <- i + 1
  }
  nsamples <- unique(sapply(lff, nrow))
  d <- as.data.frame(do.call(rbind, lff))
  d$ff <- rep(as.numeric(fetalfrac), each=nsamples)
  d <- reshape2::melt(d, id.vars = "ff")
  dsp[[dep]] <- d
}
names(dsp) <- depths

saveRDS(list(lsp, dsp), snakemake@output[["rds"]])

## lsp <- readRDS("results/summary/quilt.nipt.accuracy.regular.refsize0.chr20.rds")

bins <- sapply(lsp, function(lff) {
  o <- sapply(lff, function(d) {
    d <- d[!sapply(d[, "Mat"], is.na), ] # rm na
    d$bin
  })
  unlist(o)
})

bins <- sort(unique(as.vector(unlist(bins))))
x <- log10(as.numeric(bins))
labels <- 100 * bins

ffpchs <- c(10,15, 16, 17, 18, 19)
ffpchs <- ffpchs[seq_along(fetalfrac)]
names(ffpchs) <- fetalfrac

mycols <- c("#e69f00", "#d55e00", "#56b4e9", "#cc79a7", "#009e73", "#0072b2", "#f0e442")
palette(mycols)

## pdf(paste0(snakemake@output[["rds"]], ".pdf"), h=6,w=12)
png(paste0(snakemake@output[["rds"]], ".png"), h=12,w=6*length(depths), res = 300, units = "in")
par(mfrow=c(2,length(depths)))

for(dep in depths){
  plot(1, col = "transparent", axes = FALSE,
     xlim = c(min(x), max(x)), ylim = c(0, 1.0), cex.lab = 1.5, cex.axis = 1.5,
     ylab = "Aggregated R2 within each MAF bin", xlab = "Minor Allele Frequency %",
     main = paste0("depth:", dep))
  for(ff in fetalfrac){
    points(x, lsp[[dep]][[ff]][,"Mat"], pch = ffpchs[ff], col = mycols[1], cex = 1.3, type = "b")
    points(x, lsp[[dep]][[ff]][,"Fet"], pch = ffpchs[ff], col = mycols[2], cex = 1.3, type = "b")
  }
  axis(side = 1, at = x, labels = labels)
  axis(side = 2, at = seq(0, 1, 0.2))
  legend("bottomright", legend = names(ffpchs), pch = ffpchs, bty = "n", cex = 1.2, inset = c(0.2, 0), title = "Fetal Fraction")
  legend("bottomright", legend = c("Mat", "Fet"), fill = mycols, bty = "n", cex = 1.2, title="Haps")
}

nff <- length(fetalfrac)
a <- setdiff(seq(3*nff), seq(3,3*nff, 3))
b <- seq(1.5, 3*nff, 3)
for(dep in depths){
  d <- dsp[[dep]]
  boxplot(value ~ variable + ff, d , xaxt = "n", col = 1:2, at = a, border = NULL,
          cex.lab = 1.5, cex.axis = 1.5, xlab = "Fetal Fraction", ylab = "r2 per sample")
  axis(side = 1, at = b, labels = fetalfrac,)
  legend("bottomright", legend = c("Mother", "Kid"), fill = 1:2, bty = "n", cex=2)

}

dev.off()


