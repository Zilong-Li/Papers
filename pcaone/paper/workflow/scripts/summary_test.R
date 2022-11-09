
library(gridExtra)

l <- readRDS(snakemake@input[["log"]])
colnames(l) <- gsub("pcaonef","PCAone",colnames(l))
colnames(l) <- gsub("pcaoneh","PCAone_H+Y",colnames(l))
a <- readRDS(snakemake@input[["acc"]])

h <- a[nrow(a), ]
colnames(h) <- paste0(colnames(h), ".PCAone_H+Y")
rownames(h) <- "PCAone_H+Y"
f <- a[-nrow(a), ]
colnames(f) <- paste0(colnames(f), ".PCAone")

pdf(snakemake@output[["pdf"]], w = 10)
o <- t(cbind(l, f, h))
o <- o[order(rownames(o)),]
grid.table(round(o, digits = 6))
dev.off()
