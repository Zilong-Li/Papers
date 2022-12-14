
snakemake@source("common.R")

library(data.table)
library(grid)
library(gridExtra)
library(gtable)
library(cowplot)

inputs <- names(snakemake@input)
if ("vecplink1" %in% inputs) {
  truth <- as.matrix(fread(snakemake@input[["vecplink1"]][1])[, -c(1, 2)])
} else {
  truth <- as.matrix(fread(snakemake@input[["vecflashpca"]][1], header = T)[, -c(1, 2)])
}

l.pc <- list(
  flashpca = as.matrix(fread(snakemake@input[["vecflashpca"]][1], header = T)[, -c(1, 2)]),
  fastpca = readfastpca(snakemake@input[["vecfastpca"]][1]),
  plink2 = as.matrix(fread(snakemake@input[["vecplink2"]], header = T)[, -c(1, 2)]),
  terapca = as.matrix(fread(snakemake@input[["vecterapca"]][1], header = T)[, -1]),
  propca = as.matrix(fread(snakemake@input[["vecpropca"]][1])),
  pcaone.a = as.matrix(fread(snakemake@input[["vecpcaonea"]][1])),
  pcaone.h = as.matrix(fread(snakemake@input[["vecpcaoneh"]][1])),
  pcaone.f = as.matrix(fread(snakemake@input[["vecpcaonef"]][1]))
)

plot_pca_asia <- function(x, y, p, m, bs = 1, ...) {
  par(oma = rep(0, 4), mar = c(bs, 1, 1, 0.5), mgp = c(0, 0, 0), xaxt = "n", yaxt = "n")
  plot(x, y, bty = "l", col = p, pch = 19, ylab = "PC2", cex = 0.8, ...)
  mtext(m, side = 3, line = -0.3, font = 2)
}

data <- snakemake@wildcards[["data"]]
pops <- read.table(snakemake@config[[data]]$fam, header = F)[, 1:2]
coljco <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF", "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#4A6990FF")
if (length(unique(pops[, 1])) > 10) {
  coljco <- palette(colorAnd("large"))
}
pctruth <- function() {
  plot_pca_asia(truth[, 1], -truth[, 2], coljco[as.numeric(factor(pops[, 1]))], "Full SVD", bs = 0.1, xlab = "")
  legend("topleft", legend = sort(unique(pops[, 1])), col = coljco[as.numeric(factor(sort(unique(pops[, 1]))))], pch = 19, bty = "n", cex = 0.8, y.intersp = 0.8, inset = c(0, -0.03))
}
pcpcaone <- function() {
  plot_pca_asia(-l.pc$pcaone.f[, 1], -l.pc$pcaone.f[, 2], coljco[as.numeric(factor(pops[, 1]))], "PCAone", bs = 0.1, xlab = "")
}
pcflashpca <- function() {
  plot_pca_asia(-l.pc$flashpca[, 1], l.pc$flashpca[, 2], coljco[as.numeric(factor(pops[, 1]))], "FlashPCA2", bs = 0.1, xlab = "")
}
pcfastpca <- function() {
  plot_pca_asia(l.pc$fastpca[, 1], -l.pc$fastpca[, 2], coljco[as.numeric(factor(pops[, 1]))], "FastPCA", bs = 0.1, xlab = "")
}
pcplink2 <- function() {
  plot_pca_asia(-l.pc$plink2[, 1], -l.pc$plink2[, 2], coljco[as.numeric(factor(pops[, 1]))], "Plink2", bs = 0.1, xlab = "")
}
pcpropca <- function() {
  plot_pca_asia(l.pc$propca[, 1], l.pc$propca[, 2], coljco[as.numeric(factor(pops[, 1]))], "ProPCA", xlab = "PC1")
}
pcterapca <- function() {
  plot_pca_asia(-l.pc$terapca[, 1], l.pc$terapca[, 2], coljco[as.numeric(factor(pops[, 1]))], "TeraPCA", xlab = "PC1")
}
a0 <- function() {
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, col = "transparent", axes = F, xlab = "", ylab = "")
  text(0, 0, "Two estimated PCs (K=2)", font = 1, cex = 1.3)
}

pdf(snakemake@output[["pdf"]], h = 8, w = 6)
w <- 0.32
h <- 1 - 3 * w
ggdraw() + draw_label("B", 0.0, 0.98, 0, 0, fontface = "bold") + draw_plot(a0, 0, 0, 1, h) + draw_plot(pcterapca, 0, h, 0.5, w) + draw_plot(pcpropca, 0.5, h, 0.5, w) + draw_plot(pcflashpca, 0, w + h, 0.5, w) + draw_plot(pcplink2, 0.5, w + h, 0.5, w) + draw_plot(pctruth, 0, 2 * w + h, 0.5, w) + draw_plot(pcpcaone, 0.5, 2 * w + h, 0.5, w)
## ggdraw() + draw_label("B", 0.0, 0.98, 0, 0, fontface = "bold") + draw_plot(a0, 0, 0, 1, h) + draw_plot(pcterapca, 0, h, 0.5, w) + draw_plot(pcpropca, 0.5, h, 0.5, w) + draw_plot(pcflashpca, 0, w + h, 0.5, w) + draw_plot(pcfastpca, 0.5, w + h, 0.5, w) + draw_plot(pctruth, 0, 2 * w + h, 0.5, w) + draw_plot(pcpcaone, 0.5, 2 * w + h, 0.5, w)
dev.off()
