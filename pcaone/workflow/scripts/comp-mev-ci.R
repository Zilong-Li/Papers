snakemake@source("common.R")

source("workflow/scripts/common.R")

rule <- snakemake@rule

if (rule == "calc_mev_ci") {
  inputs <- names(snakemake@input)
  K <- snakemake@wildcards[["k"]]
  truth <- NULL
  if ("vecfull" %in% inputs)
    truth <- as.matrix(fread(snakemake@input[["vecfull"]]))[,1:K]
  if (is.null(truth))
    truth <- as.matrix(fread(snakemake@input[["vecpcaonea"]]))

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

  res.sse <- sapply(l.pc, minSSE, truth = truth)
  res.mev <- sapply(l.pc, mev, truth = truth)
  d <- as.data.frame(cbind("minSSE" = res.sse, "MEV" = res.mev))
  saveRDS(d, snakemake@output[["rds"]])
  quit()
}

if(rule == "collapse_mev_ci") {

  ## snakemake <- readRDS("results/supplementary/suppl_mev_ci/mev.ci.k10.rds")

  res <- lapply(snakemake@input[["rds"]], readRDS)

  prog <- rownames(res[[1]])

  d <- sapply(res, "[[", "MEV")
  ## d <- sapply(res, "[[", "minSSE")

  pdf(paste0(snakemake@output[["rds"]], ".pdf"), w = 9, h = 6)
  boxplot(t(d), names = fancyname2(prog), ylab = "MEV(K=10)", main = "10 random subsets of Asian data with 1 millions SNPs")
  dev.off()

  saveRDS(snakemake, snakemake@output[["rds"]])
  quit()
}

## iram <- as.matrix(fread("results/asian/pcaone.a.k40.eigvecs", h = F, data.table = F))
## pcaone <- as.matrix(fread("results/asian/pcaone.f.k40.eigvecs", h = F, data.table = F))
## halko <- as.matrix(fread("results/asian/pcaone.h.k40.eigvecs", h = F, data.table = F))

## norm_vec <- function(x) {
##   sqrt(sum(x^2))
## }

## mev <- function(Y, truth) {
##   if (!isTRUE(all.equal(dim(Y), dim(truth)))) {
##     return(NA)
##   }
##   mean(apply(truth, 2, function(x) {
##     norm_vec(t(Y) %*% x)
##   }))
## }

## (t(Y[sam,1:2]) %*% truth[sam,1:2])

## new <- iram

## mean(Y[,1])

## rmse(Y, truth)


## sam <- 1:100
## Y <- pcaone[sam,1:2]
## truth <- new[sam,1:2]


## new[,2] <- -new[,2]

## mev(pcaone[sam,1:2], iram[sam,1:2])

## mev(pcaone[sam,1:2], new[sam,1:2])


## dim(pcaone)


## mev.ci <- function(x, truth, iter = 100, k = NULL) {
##   n <- dim(x)[1]
##   o <- sapply(1:iter, function(i) {
##     sam <- sample(1:n, replace = T)
##     if (is.null(k)) {
##       mev(x[sam, ], truth[sam, ])
##     } else {
##       mev(x[sam, 1:k], truth[sam, 1:k])
##     }
##   })
##   o
## }

## (res <- mev.ci(pcaone, truth = iram, iter = 1000, k = 2))
## (SE <- sd(res))

## (CI95 <- est + c(-1, 1) * qnorm(0.975) * SE)

## c(est = est, CI = CI95)

## (CI95 <- quantile(res,c(0.025,0.975)))

## c(est=mev,CI=CI95)

## dev.off()
