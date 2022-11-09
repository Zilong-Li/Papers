
snakemake@source("common.R")

library(data.table)
library(gridExtra)

truth <- as.matrix(fread(snakemake@input[["pcaonea"]]))

pcaone.h <- as.matrix(fread(snakemake@input[["pcaoneh"]]))
pcaone.f <- as.matrix(fread(snakemake@input[["pcaonef"]]))
terapca <- as.matrix(fread(snakemake@input[["terapca"]]))

ll <- list(pcaone.f, pcaone.h, terapca)

res.sse <- sapply(ll, minSSE, truth = truth)
res.mev <- sapply(ll, mev, truth = truth)
res <- data.frame(minSSE = res.sse, MEV = res.mev)

saveRDS(res, snakemake@output[["rds"]])

pdf(snakemake@output[["pdf"]])
grid.table(as.data.frame(res))
dev.off()
