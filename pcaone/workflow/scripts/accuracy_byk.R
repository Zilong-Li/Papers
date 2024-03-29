
snakemake@source("common.R")

library(data.table)
library(grid)
library(gridExtra)
library(gtable)
library(cowplot)


inputs <- names(snakemake@input)
K <- snakemake@wildcards[["k"]]
truth <- NULL
if ("vecfull" %in% inputs)
  truth <- as.matrix(fread(snakemake@input[["vecfull"]]))[,1:K]
## if ("vecplink1" %in% inputs) {
##   truth <- readplink1(snakemake@input[["vecplink1"]])
## }
if (is.null(truth)) {
  truth <- as.matrix(fread(snakemake@input[["vecpcaonea"]]))
}
## truth <- as.matrix(fread(snakemake@input[["vecflashpca"]], header = T)[, -c(1, 2)])

l.pc <- list(
  flashpca = readflashpca(snakemake@input[["vecflashpca"]]),
  ## fastpca = readfastpca(snakemake@input[["fastpca"]]),
  plink2 = readplink2(snakemake@input[["vecplink2"]]),
  terapca = as.matrix(fread(snakemake@input[["vecterapca"]], header = T)[, -1]),
  propca = as.matrix(fread(snakemake@input[["vecpropca"]])),
  pcaone.a = as.matrix(fread(snakemake@input[["vecpcaonea"]])),
  pcaone.h = as.matrix(fread(snakemake@input[["vecpcaoneh"]])),
  pcaone.f = as.matrix(fread(snakemake@input[["vecpcaonef"]]))
)

l.log <- list(
  flashpca = readLines(snakemake@input[["logflashpca"]]),
  plink2 = readLines(snakemake@input[["logplink2"]]),
  terapca = readLines(snakemake@input[["logterapca"]]),
  propca = readLines(snakemake@input[["logpropca"]]),
  pcaone.a = readLines(snakemake@input[["logpcaonea"]]),
  pcaone.h = readLines(snakemake@input[["logpcaoneh"]]),
  pcaone.f = readLines(snakemake@input[["logpcaonef"]])
)

times <- round(logtimes(l.log)) # in secs
rams <- logram(l.log) # in kbytes

res.sse <- sapply(l.pc, minSSE, truth = truth)
res.mev <- sapply(l.pc, mev, truth = truth)
d <- as.data.frame(cbind("minSSE" = res.sse, "MEV" = res.mev, "Time.s" = times, "RAM.kb" = rams))
saveRDS(d, snakemake@output[["rds"]])


# https://cran.r-project.org/web/packages/gridExtra/vignettes/tableGrob.html

# right-aligned
tt2 <- ttheme_default(
  core = list(fg_params = list(hjust = 1, x = 0.9)),
  rowhead = list(fg_params = list(hjust = 1, x = 0.95))
)
# left-aligned
tt3 <- ttheme_default(
  core = list(fg_params = list(hjust = 0, x = 0.1)),
  rowhead = list(fg_params = list(hjust = 0, x = 0))
)

pdf(snakemake@output[["pdf"]])
grid.table(round(d, digits = 6), theme = tt2) # right-aligned
dev.off()
