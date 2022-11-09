
snakemake@source("common.R")

library(gridExtra)

windows <- snakemake@params[["windows"]]

pcaone.h = list(pcaoneh = readLines(snakemake@input[["pcaoneh"]]))
pcaone.f = lapply(snakemake@input[["pcaonef"]], readLines)

time.pcaonef <- logtimes(pcaone.f)
ram.pcaonef <- logram(pcaone.f)
time.pcaoneh <- logtimes(pcaone.h)
ram.pcaoneh <- logram(pcaone.h)
epochs.pcaonef <- logepochs(pcaone.f)
epochs.pcaoneh <- logepochs(pcaone.h)

d <- data.frame(time.pcaonef, time.pcaoneh, epochs.pcaonef, epochs.pcaoneh, ram.pcaonef, ram.pcaoneh, row.names=paste0("windows=",windows))
saveRDS(d, snakemake@output[["rds"]])
pdf(snakemake@output[["pdf"]], w = 10)
grid.table(t(d))
dev.off()

## save.image()
