
snakemake@source("common.R")


A <- readRDS(snakemake@input[["A"]])
B <- readRDS(snakemake@input[["B"]])


A.acc <- do.call(cbind.data.frame, A$acc)
A.log <- do.call(cbind.data.frame, A$log)
A.names <- gsub("data","",names(A$log))
A.time <- A.log[, seq(1, ncol(A.log), 2)]
colnames(A.time) <- A.names
A.ram <- A.log[, seq(2, ncol(A.log), 2)]
colnames(A.ram) <- A.names

B.acc <- do.call(cbind.data.frame, B$acc)
B.log <- do.call(cbind.data.frame, B$log)
B.names <- gsub("data","",names(B$log))
B.time <- B.log[, seq(1, ncol(B.log), 2)]
colnames(B.time) <- B.names
B.ram <- B.log[, seq(2, ncol(B.log), 2)]
colnames(B.ram) <- B.names

fig2.names <- c(A.names, B.names)

fig2.time <- cbind(A.time, B.time)
fig2.time.sp <- 1/sweep(fig2.time, 2, as.numeric(fig2.time[1,]), "/")
## fig2.time.min <- fig2.time / 60

fig2.ram <- cbind(A.ram, B.ram)
fig2.ram.sp <- sweep(fig2.ram, 2, as.numeric(fig2.ram[1,]), "/")
## fig2.ram.gb <- fig2.ram / 1024^2

fig2.acc <- cbind(A.acc, B.acc)
fig2.acc.minsse <- fig2.acc[, seq(1, ncol(fig2.acc), 2)]
colnames(fig2.acc.minsse) <- fig2.names
fig2.acc.mev <- fig2.acc[, seq(2, ncol(fig2.acc), 2)]
colnames(fig2.acc.mev) <- fig2.names

res <- list("time" = fig2.time, "ram" = fig2.ram, "minSSE" = fig2.acc.minsse, "MEV" = fig2.acc.mev)
progs <- rownames(res$time)
saveRDS(res, snakemake@output[["rds"]])

plotmat <- function(mat, xlabels, ...) {
  par(mar = c(4.5, 4.5, 0.5, 0))
  plot(1, bty = "l", type = "n", xaxt = "n", col = "transparent", cex.main = 2.0, cex.axis = 1.5, cex.lab = 2.0, xlim = c(1, nrow(mat)), ...)
  mapprog <- PROGS # inject global variable
  for (n in colnames(mat)) {
    points(1:6, mat[1:6, n], type = "b", col = mapprog['col', n], lwd = 2)
    points(7:13, mat[7:13, n], type = "b", col = mapprog['col', n], lwd = 2)
  }
  axis(1, at = seq_len(nrow(mat)), labels = xlabels, cex.axis = 1.5)
}
plotfig2 <- function(fig2,...) {
    plotmat(t(fig2), xlabels = colnames(fig2), ylim = range(fig2, na.rm = T),  ...)
}

mplot <- function(mat, x, labels = TRUE, ylim = range(mat, na.rm = T), ...) {
  plot(x, mat[,1], bty = "l",  xaxt = "n",col = "transparent", ylim = ylim, cex.main = 2.0, cex.axis = 1.5, cex.lab = 1.5, ...)
  ## plot(1, bty = "l",  col = "transparent", ylim = ylim, cex.main = 2.0, cex.axis = 1.5, cex.lab = 2.0, ...)
  mapprog <- PROGS # inject global variable
  for (n in colnames(mat)) {
    points(x, mat[, n], type = "b", col = mapprog['col', n], lwd = 2)
  }
  axis(1, at = x, labels = labels, cex.axis = 1.5)
}

mplot.t <- function(mat, x, ...) {
    mplot(t(mat), x, ...)
}

plegend <- function() {
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, col = "transparent", axes = F, xlab = "", ylab = "")
  name <- rev(rownames(res$ram))
  legend("bottomleft", legend = fancyname(name), fill = as.character(PROGS["col", name]), horiz = T, cex = 1.3, xjust = 0, yjust = 0, xpd = T, inset = c(-0.01, 0), bty = "n",text.font = 2, x.intersp = 0.5, adj = 0.02)
## legend("topleft", legend = as.character(PROGS["name", rev(rownames(res$ram))]), col = as.character(PROGS["col", rev(rownames(res$ram))]), pch = 21, lwd = 1.5, cex = 1.5, bty = "n", y.intersp = 1, inset = c(0, -0.03))
}


## setwd('/Users/zilong/Project/zll/shangban/pcaone/paper/workflow/scripts/')
## res <- readRDS("~/Downloads/fig2.k20.rds")
## source("common.R")

pdf(snakemake@output[["pdf"]], w = 9, h = 6)
## pdf("~/Downloads/t.pdf", w = 9, h = 6)
layout(matrix(c(1, 1, 2, 3, 4, 5), nrow = 3, byrow = T), widths = c(1, 1), heights = c(1, 3, 5))
## layout.show(5)
plegend()
par(mar = c(1, 4.5, 0.5, 0))
mplot.t(res$ram[,7:13]/1024^2, x = c(2,4,6,8,10,20,30), labels  = FALSE, xlab = "",ylab = "Memory in GBs", log = 'x')
legend("topleft", "#Samples=50,000", bty="n", cex=1.5)
par(mar = c(1, 4.5, 0.5, 0))
mplot.t(res$ram[,1:6]/1024^2, x = c(1,2,4,6,8,10), labels  = FALSE, xlab = "",ylab = "", log = 'x')
legend("topleft", "#SNPs=3,000,000", bty="n", cex=1.5)
par(mar = c(4.5, 4.5, 0.5, 0))
mplot.t(res$time[,7:13]/60, x = c(2,4,6,8,10,20,30), xlab = "N (x100000) SNPs",ylab = "Time in Minutes", log = 'xy')
par(mar = c(4.5, 4.5, 0.5, 0))
mplot.t(res$time[,1:6]/60, x = c(1,2,4,6,8,10), xlab = "N (x10000) Samples",ylab = "", log = 'xy')
dev.off()

## pdf("figure2.pdf", w = 14, h = 12)
## pdf(snakemake@output[["pdf"]], w = 12, h = 6)
## par(mfrow = c(1, 2))
## ## plotfig2(fig2.time.sp, xlab = "", ylab = "Time Speedup (x)")
## progs <- rev(progs)
## plotfig2(res$time/60, xlab = "",ylab = "Time in Minutes scaled by log10", log = 'y')
## legend("topleft", legend = as.character(PROGS["name", progs]), col = as.character(PROGS["col", progs]), pch = 21, lwd = 1.5, cex = 2, bty = "n", y.intersp = 1, inset = c(0, -0.03))
## ## plotfig2(res$MEV, xlab = "",ylab = "MEV of top K=20" )
## ## plotfig2(fig2.ram.sp, xlab = "Subset of UK Biobank datat",ylab = "Relative RAM allocation (x)")
## plotfig2(res$ram/1024^2,xlab = "Subset of UK Biobank datat", ylab = "Memory in GBs scaled by log10", log = 'y')
## ## plotfig2(res$minSSE,xlab = "Subset of UK Biobank datat", ylab = "minSSE of top K=20", log = 'y')
## warnings()
## dev.off()
