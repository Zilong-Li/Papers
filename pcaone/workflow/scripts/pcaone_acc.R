
snakemake@source("common.R")

library(data.table)
library(gridExtra)

windows <- snakemake@params[["windows"]]

## truth <- as.matrix(fread(snakemake@input[["plink1"]])[, -c(1,2)])
truth = as.matrix(fread(snakemake@input[["pcaonea"]]))

pcaone.h = as.matrix(fread(snakemake@input[["pcaoneh"]]))
## pcaone.f = as.matrix(fread(snakemake@input[["pcaonef"]]))
pcaone.f = lapply(snakemake@input[["pcaonef"]], function(fn) as.matrix(fread(fn)))

res.sse <- sapply(pcaone.f, minSSE, truth = truth)
res.mev <- sapply(pcaone.f, mev, truth = truth)
res.pcaonef <- data.frame(minSSE = res.sse, MEV = res.mev, row.names = paste0("pcaonef.w",windows))
res.sse <- minSSE(pcaone.h, truth = truth)
res.mev <- mev(pcaone.h, truth = truth)
res.pcaoneh <- data.frame(minSSE = res.sse, MEV = res.mev, row.names = "pcaoneh")

d <- rbind(res.pcaonef, res.pcaoneh)
saveRDS(d, snakemake@output[["rds"]])

pdf(snakemake@output[["pdf"]])
grid.table(as.data.frame(d))
dev.off()
