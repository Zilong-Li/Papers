
saveRDS(snakemake, snakemake@output[["rds"]])
snakemake@source("common.R")

## quit()
## source("workflow/scripts/common.R")
## snakemake <- readRDS("results/supplementary/stop_criteria//asian.acc.k10.bfile.rds")
## snakemake <- readRDS("results/supplementary/stop_criteria//1000G.acc.k10.bfile.rds")

full = as.matrix(fread(snakemake@input[["full"]]))
## full <- as.matrix(fread(cmd = "cut -d' ' -f3- results/asian/plink1.k40.eigenvec"))

ff <- snakemake@input[["pcaonef"]]
fh <- snakemake@input[["pcaoneh"]]

ff <- sub(pattern = ".eigvecs", replacement = "", x = ff)
fh <- sub(pattern = ".eigvecs", replacement = "", x = fh)

maxp <- 40

ffw <- lapply(ff, function(f) {
  paste0(f, ".epoch.", 1:maxp-1, ".eigvecs")
})
names(ffw) <- paste0("w=",snakemake@params$windows)


fhw <- lapply(fh, function(f) {
  paste0(f, ".epoch.", 1:maxp-1, ".eigvecs")
})


get_mev_true <- function(fh, full) {
  h_mev_true <- sapply(fh, function(f) {
    t <- as.matrix(fread(f))
    mev(t, full[,1:ncol(t)])
  })
  h_mev_true
}


get_mev_succ <- function(fh, maxp) {
  lfh <- list(fh[1:(maxp-1)], fh[2:maxp])
  h_mev_succ <- sapply(1:(maxp-1), function(i) {
    t1 <- as.matrix(fread(lfh[[1]][i]))
    t2 <- as.matrix(fread(lfh[[2]][i]))
    mev(t1, t2)
  })
}

f_mev_true <- get_mev_true(ffw[["w=64"]], full)
f_mev_succ <- get_mev_succ(ffw[["w=64"]], maxp)

(f_mev <- cbind(log10(1-f_mev_succ),  log10(1-f_mev_true[-1])))

h_mev_true <- get_mev_true(fhw[[1]], full)
h_mev_succ <- get_mev_succ(fhw[[1]], maxp)

(h_mev <- cbind(log10(1-h_mev_succ),  log10(1-h_mev_true[-1])))

plot_mev_per_epoch <- function(d, legend = TRUE, ...) {
  d[is.infinite(d)] <- -6
  d[d[,1] < -6.0, 1] <- -6
  d[d[,2] < -6.0, 2] <- -6
  plot(d[,2], col = "blue", type = "b", ylab = "log10(1-MEV)", xlab = "Epoch", ylim = c(-6, 0), xaxt = "n", cex = 1.5, ...)
  lines(d[,1], col = "orange", type = "b", cex = 1.5)
  axis(1, at = 1:nrow(d), labels = 1:nrow(d)+1)
  abline(h = -4, lty = 2, col = "gray60")
  x <- which(d[,1] < -4)[1]
  y <- d[x, 1]
  text(x = x, y = y, labels = "X", col = "red", cex = 1.2)
  if(legend)
    legend("topright", legend = c("MEV between current and truth", "MEV between successive epochs"), col = c("blue", "orange"), pch = 21, cex = 1.3)
}

pdf(paste0(snakemake@output[["rds"]], ".pdf"), h = 6, w = 12)
par(mfrow = c(1, 2))
plot_mev_per_epoch(f_mev, main = expression(bold("PCAone")))
plot_mev_per_epoch(h_mev, legend = FALSE, main = expression(bold("PCAone"[H+Y])))
dev.off()
