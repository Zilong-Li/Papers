
snakemake@source("common.R")

library(data.table)
library(grid)
library(gridExtra)
library(gtable)
library(cowplot)

KK <- snakemake@params[["pcs"]]

pcaone.a <- sapply(snakemake@input[["pcaonea"]], function(f) {
  lines <- readLines(f)
  r <- strsplit(lines[grep("Matrix Operation", lines)], " ")
  as.numeric(rev(r[[length(r)]])[1])
})

flashpca <- sapply(snakemake@input[["flashpca"]], function(f) {
  lines <- readLines(f)
  r <- strsplit(lines[grep("Matrix operation", lines)], " ")
  as.numeric(rev(r[[length(r)]])[1])
})

pcaone.h <- sapply(snakemake@input[["pcaoneh"]], function(f) {
  lines <- readLines(f)
  r <- strsplit(lines[grep("stops", lines)], " ")[[1]]
  as.numeric(gsub("epoch=", "", r[length(r)]))
})

pcaone.f <- sapply(snakemake@input[["pcaonef"]], function(f) {
  lines <- readLines(f)
  r <- strsplit(lines[grep("stops", lines)], " ")[[1]]
  as.numeric(gsub("epoch=", "", r[length(r)]))
})

terapca <- sapply(snakemake@input[["terapca"]], function(f) {
  lines <- readLines(f)
  r <- strsplit(lines[grep("converged", lines)], " ")
  n <- as.numeric(gsub(",", "", r[[length(r)]][3])) + 1
  n * 2
  ## paste0(n, "*")
})

propca <- sapply(snakemake@input[["propca"]], function(f) {
  lines <- readLines(f)
  r <- strsplit(lines[grep("^Iteration", lines)], " ")
  n <- as.numeric(r[[length(r)]][2])
  n
  ## paste0(n, "*")
})

## plink2 <- paste0(KK, "*")
plink2 <- KK

epochs <- data.frame(pcaone.f, pcaone.h, terapca, plink2, propca, flashpca, pcaone.a)
rownames(epochs) <- KK

acc <- lapply(snakemake@input[["rds"]], function(f) {
  readRDS(f)[, c("minSSE", "MEV")]
})
names(acc) <- KK

log <- lapply(snakemake@input[["rds"]], function(f) {
  readRDS(f)[, c("Time.s", "RAM.kb")]
})
names(log) <- KK


# now show epochs vs mev for pcaone.f and pcaone.h for K=10
## K <- "10"
K <- "40"
iK <- which(KK == K)
## file.plink1 <- snakemake@input[["vecplink1"]][[iK]]
## truth <- as.matrix(fread(file.plink1)[, -c(1, 2)])
file.pcaone.full <- snakemake@input[["vecfull"]]
truth <- as.matrix(fread(file.pcaone.full)[,1:K])
file.pcaone.h <- snakemake@input[["vecpcaoneh"]][[iK]]
file.pcaone.f <- snakemake@input[["vecpcaonef"]][[iK]]
prefix <- gsub(".eigvecs", "", file.pcaone.h)
C.pcaone.h <- sapply(1:epochs[iK, "pcaone.h"], function(i) {
  test <- as.matrix(fread(paste0(prefix, ".epoch.", i - 1, ".eigvecs")))
  mev(test, truth)
})
prefix <- gsub(".eigvecs", "", file.pcaone.f)
C.pcaone.f <- sapply(1:epochs[iK, "pcaone.f"], function(i) {
  test <- as.matrix(fread(paste0(prefix, ".epoch.", i - 1, ".eigvecs")))
  mev(test, truth)
})
C <- list("pcaone.f" = C.pcaone.f, "pcaone.h" = C.pcaone.h)

# save all summary data
saveRDS(list(acc = acc, log = log, epochs = epochs, C = C), snakemake@output[["rds"]])

# right-aligned
tt2 <- ttheme_default(
  core = list(fg_params = list(hjust = 1, x = 0.9)),
  rowhead = list(fg_params = list(hjust = 1, x = 0.95))
)

pdf(snakemake@output[["pdf2"]], h = 6, w = 12)
gtbls <- lapply(log, function(d) tableGrob(round(d, digits = 6), theme = tt2))
plot_grid(plotlist = gtbls, labels = paste0("K", names(gtbls)))
dev.off()


## setwd("~/Project/zll/shangban/pcaone/paper/workflow/scripts/")
## source("common.R")
## rds <- readRDS("~/Downloads/fig1A.asian.rds")
## C <- rds$C
## C$pcaone.f
## epochs <- rds$epochs
## KK <- c(2, 5, 10, 20, 30, 40)
## acc <- rds$acc
## log <- rds$log

## ram.max <- apply(m.ram, 1, max)
## barplot(ram.max/1024^2,  horiz = T, names.arg = PROGS['name', names(ram.max)])

# https://stackoverflow.com/questions/4227223/convert-a-list-to-a-data-frame
acc <- do.call(cbind.data.frame, acc) # to data.frame. make wider
m.sse <- as.matrix(acc[, seq(1, ncol(acc), 2), drop = FALSE])
m.mev <- as.matrix(acc[, seq(2, ncol(acc), 2), drop = FALSE])
m.sse <- ifelse(m.sse <= 1e-10, 1e-10, m.sse)
progs <- rev(rownames(m.mev))
log <- do.call(cbind.data.frame, log)
m.ram <- as.matrix(log[, seq(2, ncol(log), 2), drop = FALSE])/1024^2
colnames(m.ram) <- KK

plotacc <- function(mat, xlabels, mapprog, ...) {
  par(mar = c(3.5, 3.5, 2, 1), mgp = c(2, 0.6, 0))
  plot(1, bty = "l", type = "n", xaxt = "n", col = "transparent", cex.main = 2.0, cex.axis = 1.2, cex.lab = 1.2, xlim = c(1, nrow(mat)), ...)
  for (n in colnames(mat)) {
    points(mat[, n], type = "b", col = mapprog["col", n], lwd = 2)
  }
  axis(1, at = 1:nrow(mat), labels = xlabels, cex.axis = 1.2)
}

plotacc2 <- function(h, f, ...) {
  par(mar = c(3.5, 3.5, 2, 1), mgp = c(2, 0.6, 0))
  plot(1, bty = "l", type = "n", xaxt = "n", col = "transparent", cex.main = 2.0, cex.axis = 1.2, cex.lab = 1.2, xlim = c(1, length(h)), ...)
  abline(h = 1, col = "gray", lwd = 2, lty = 2)
  points(h, type = "b", col = PROGS["col", "pcaone.h"], lwd = 2)
  points(f, type = "b", col = PROGS["col", "pcaone.f"], lwd = 2)
  axis(1, at = seq(1, length(h) + 1, 4), cex.axis = 1.2)
}

figE <- function() {
  plotacc2(C$pcaone.h, C$pcaone.f, ylim = c(0.2, 1.0), xlab = "Epochs", ylab = paste0("Accuracy (MEV) of K=",K," PCs"))
  legend("bottomright", legend = fancyname(c("pcaone.f", "pcaone.h")), col = as.character(PROGS["col", c("pcaone.f", "pcaone.h")]), pch = 21, lwd = 1.3, cex = 1.3, bty = "n")
}

figD <- function() {
    plotacc(t(m.mev), xlabels = KK, mapprog = PROGS, ylim = range(m.mev, na.rm = T), xlab = "Number of estimated PCs (K)", ylab = "Accuracy (MEV) of PCs")
}

figC <- function() {
  plotacc(t(m.ram), xlabels = KK, mapprog = PROGS, ylim = range(m.ram, na.rm = T), xlab = "Number of estimated PCs (K)", ylab = "Memory in Gigabytes")
  legend("topleft", legend =  fancyname(progs), col = as.character(PROGS["col", progs]), pch = 21, lwd = 1.3, cex = 1.2, bty = "n", inset = c(0, -0.05))
}

figA <- function() {
  name <- colnames(epochs)
  epochs <- cbind("K" = KK, epochs)
  colnames(epochs) <- c(expression(bold("K(PCs)")) , fancyname(name, bold = T))
  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
  g <- tableGrob(epochs, rows = NULL, theme=tt)
  g$widths <- unit(rep(1 / ncol(g) * 0.95, ncol(g)), "npc")
  g <- gtable_add_grob(g,
    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
    t = 2, b = nrow(g), l = 1, r = ncol(g)
  )
  g <- gtable_add_grob(g,
    grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
    t = 1, l = 1, r = ncol(g)
  )
  g
}

figAtitle <- function() {
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, col = "transparent", axes = F, xlab = "", ylab = "")
  text(0, 0, "Number of passes over data (epochs) of each method for top K PCs", font = 2, cex = 1.3)
}

## ggdraw() + draw_plot(figC, 0, 0, 0.33, 0.65) + draw_plot(figD, 0.33, 0, 0.33, 0.65) + draw_plot(figE, 0.66, 0, 0.33, 0.65)+ draw_plot(figA(), 0, 0.61, 1, 0.35) + draw_plot(figAtitle, 0, 0.9, 1, 0.15) + draw_label("A", 0.01, 0.97, 0, 0, fontface = "bold") + draw_label("C", 0.01, 0.58, 0, 0, fontface = "bold") + draw_label("D", 0.34, 0.58, 0, 0, fontface = "bold")+ draw_label("E", 0.67, 0.58, 0, 0, fontface = "bold")

## pdf("t.pdf", w = 12, h = 8)
## ggdraw() + draw_plot(figC, 0, 0, 0.5, 0.7) + draw_plot(figD, 0.5, 0, 0.5, 0.7) + draw_plot(figA(), 0, 0.65, 1, 0.35) + draw_plot(figAtitle, 0, 0.89, 1, 0.15) + draw_label("A", 0.01, 0.97, 0, 0, fontface = "bold") + draw_label("C", 0.01, 0.65, 0, 0, fontface = "bold") + draw_label("D", 0.51, 0.65, 0, 0, fontface = "bold")

## pdf("~/Downloads/t.pdf", h = 6, w = 12)
pdf(snakemake@output[["pdf"]], h = 6, w = 12)
ggdraw() + draw_plot(figC, 0, 0, 0.33, 0.65) + draw_plot(figD, 0.33, 0, 0.33, 0.65) + draw_plot(figE, 0.66, 0, 0.33, 0.65)+ draw_plot(figA(), 0, 0.61, 1, 0.35) + draw_plot(figAtitle, 0, 0.9, 1, 0.15) + draw_label("A", 0.01, 0.97, 0, 0, fontface = "bold") + draw_label("C", 0.01, 0.58, 0, 0, fontface = "bold") + draw_label("D", 0.34, 0.58, 0, 0, fontface = "bold")+ draw_label("E", 0.67, 0.58, 0, 0, fontface = "bold")
dev.off()
