snakemake@source("common.R")

library(data.table)
library(gridExtra)
library(parallel)

windows <- snakemake@params[["windows"]]


truth <- as.matrix(fread(snakemake@input[["pcaonea"]]))
pcaone.h <- as.matrix(fread(snakemake@input[["pcaoneh"]]))
pcaone.f <- as.matrix(fread(snakemake@input[["pcaonef"]][[1]]))
onlinepca <- as.matrix(fread(snakemake@input[["onlinepca"]]))

# order matters!!!!
ll <- list(pcaone.f, pcaone.h, onlinepca, "pcaone.a" = truth)
res.sse <- unlist(mclapply(ll, minSSE, truth = truth, mc.cores = 4))
res.mev <- unlist(mclapply(ll, mev, truth = truth, mc.cores = 4))

llog <- list(
  pcaone.f = readLines(snakemake@input[["logpcaonef"]][[1]]),
  pcaone.h = readLines(snakemake@input[["logpcaoneh"]]),
  onlinepca = readLines(snakemake@input[["logonlinepca"]]),
  pcaone.a = readLines(snakemake@input[["logpcaonea"]])
)

times <- round(logtimes(llog)) # in secs
rams <- logram(llog) # in kbytes


res <- data.frame("minSSE" = res.sse, "MEV" = res.mev, "Time.s" = times, "RAM.kb" = rams)
rownames(res) <- c("PCAone", "PCAone_H+Y", "OnlinePCA_Halko", "PCAone_Arnoldi")
## save.image()

onlinepca <- flip(onlinepca, truth)
pcaone.f <- flip(pcaone.f, truth)
pcaone.h <- flip(pcaone.h, truth)

Af <- matchMin(x = pcaone.f, y = truth)
Ah <- matchMin(x = pcaone.h, y = truth)
Ao <- matchMin(x = onlinepca, y = truth)

byk.sse <- cbind(
  "PCAone" = colSums((pcaone.f - Af)^2),
  "PCAone_H+Y" = colSums((pcaone.h - Ah)^2),
  "OnlinePCA_Halko" = colSums((onlinepca - Ao)^2)
)
saveRDS(list("all" = res, "byk" = byk.sse), snakemake@output[["rds"]])

pdf(snakemake@output[["pdf"]])
grid.table(as.data.frame(res))
dev.off()

## ccol <- c("darkred","darkblue","darkgreen")
## barplot(t(rmsd),beside=T,col=ccol,xlab="PCs",names=1:40,ylab="minSSE per PC",main="Brain scRNAs Accuracy")
