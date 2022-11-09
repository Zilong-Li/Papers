
snakemake@source("common.R")

library(grid)
library(gridExtra)
library(cowplot)

datas <- snakemake@params[["dataset"]]
K <- snakemake@wildcards[["k"]]
## res <- lapply(snakemake@input, readRDS)
res <- lapply(snakemake@input, function(f) {
    readRDS(f)[, c("minSSE","MEV")]
})
names(res) <- datas

log <- lapply(snakemake@input, function(f) {
    readRDS(f)[, c("Time.s","RAM.kb")]
})
names(log) <- datas

saveRDS(list(acc = res, log = log), snakemake@output[["rds"]])

# https://stackoverflow.com/questions/4227223/convert-a-list-to-a-data-frame
res <- do.call(cbind.data.frame, res) # to data.frame. make wider
m.sse <- as.matrix(res[,seq(1, ncol(res), 2), drop = FALSE])
m.mev <- as.matrix(res[,seq(2, ncol(res), 2), drop = FALSE])

# right-aligned
tt2 <- ttheme_default(core=list(fg_params=list(hjust=1, x=0.9)),
                      rowhead=list(fg_params=list(hjust=1, x=0.95)))

pdf(snakemake@output[["pdf"]], h = 6, w = 12)
palette(wong)
par(mfrow = c(1, 2))
progs <- rev(rownames(m.mev))
plotacc(t(m.mev), xlabels = datas, ylim=range(m.mev, na.rm = T), xlab = "Program setting K", ylab = paste0("MEV of top ", K, " PCs"))
legend("bottomright", legend = as.character(PROGS["name", progs]), col = as.character(PROGS["col", progs]), pch = 21, lwd = 1.3)
plotacc(t(m.sse), xlabels = datas, ylim=range(m.sse, na.rm = T), xlab = "Program setting K", ylab = paste0("minSSE of top ", K, " PCs"))
dev.off()
## legend("topleft", rownames(m.sse), col = 1:nrow(m.sse), pch = 21, lwd = 1.3)

pdf(snakemake@output[["pdf2"]], h = 12, w = 12)
gtbls <- lapply(log, function(d) tableGrob(round(d, digits=6),theme = tt2))

plot_grid(plotlist=gtbls, labels = names(gtbls), ncol = 3)


dev.off()

