## snakemake@source("common.R")

## source("workflow/scripts/common.R")

rule <- snakemake@rule
scenario <- snakemake@config["scenario"]

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442")
palette(cbPalette)

plot_mev_byk <- function(d.mev, ymin = NULL, legend = TRUE) {
  if(is.null(ymin))
    ymin <- min(d.mev)
  n <- ncol(d.mev)
  data <- rownames(d.mev)
  plot(1, bty = "l", type = "n", xaxt = "n", col = "transparent", cex.main = 2.0, cex.axis = 1.2, cex.lab = 1.2, xlim = c(1, n), ylim = c(ymin, 1.0), ylab = "MEV", xlab = "")
  for (d in data) {
    points(d.mev[d, ], type = "b", lwd = 2, col = which(factor(data) == d))
  }
  axis(1, at = 1:n, labels = paste0("K=",colnames(d.mev)), cex.axis = 1.2)
  if(legend)
    legend("bottomright", legend = data, fill = factor(data), cex = 1.2, bty = "n")
}

plot_eigvals <- function(d.mev, ymin = NULL, legend = TRUE) {
  if(is.null(ymin))
    ymin <- min(d.mev)
  n <- ncol(d.mev)
  data <- rownames(d.mev)
  plot(1, ylim = c(min(d.mev), max(d.mev)), log = "y", bty = "l", type = "n", xaxt = "n", col = "transparent", cex.main = 2.0, cex.axis = 1.2, cex.lab = 1.2, xlim = c(1, n), ylab = "Eigenvalues", xlab = "PCs")
  for (d in data) {
    points(d.mev[d, ], type = "b", lwd = 2, col = which(factor(data) == d))
  }
  ## axis(1, at = 1:n, labels = paste0("K=",colnames(d.mev)), cex.axis = 1.2)
  if(legend)
    legend("topright", legend = data, fill = factor(data), cex = 1.2, bty = "n")
}


## this is by PCs
if (rule %in% c("collapse_bgen_summary",  "collapse_binary_summary",  "collapse_csv_summary")) {
  log <- lapply(snakemake@input[["log"]], readRDS)
  names(log) <- snakemake@params[["k"]]
  acc <- lapply(snakemake@input[["acc"]], readRDS)
  names(acc) <- snakemake@params[["k"]]
  res <- list(log = log, acc = acc)
  saveRDS(res, snakemake@output[["rds"]])

  ## res <- readRDS("results/supplementary/suppl_rnaseq/summary.gtex_rnaseq.binary.rds")
  ## acc <- res[["acc"]]
  ## log <- res[["log"]]

  d.epochs <- sapply(log, function(d) {
    as.numeric(d[,3:4])
  })
  rownames(d.epochs) <- c("PCAone", "PCAone_H+Y")
  colnames(d.epochs) <- paste0("K=", colnames(d.epochs))

  d.mev <- sapply(acc, function(d) {
    d[,2]
  })
  rownames(d.mev) <- c("PCAone", "PCAone_H+Y")

  pdf(snakemake@output[["pdf"]], w = 14, h = 7)
  par(mfrow = c(1, 2))
  if (scenario == "suppl_rnaseq")
    barplot(d.epochs, beside = T, col = factor(rownames(d.mev)), ylab = "Epochs",  main = "GTEx bulk RNA-seq data (M=56200, N=17382)")
  if (scenario == "suppl_meta")
    barplot(d.epochs, ylim = c(0, max(d.epochs)+3),beside = T, col = factor(rownames(d.mev)), ylab = "Epochs",  main = "MetaHit data (M=192140, N=528)")
  legend("topleft",bty = "n", cex = 1.5, legend = c(expression("PCAone"), expression(paste("PCAone"[H + Y]))), fill = factor(rownames(d.epochs)))
  plot_mev_byk(d.mev, ymin = min(d.mev), legend = FALSE)
  dev.off()

}

if (rule == "collapse_bfile_summary") {

  ## saveRDS(snakemake, snakemake@output[["rds"]])
  ## pdf(snakemake@output[["pdf1"]])
  ## plot(1)
  ## dev.off()
  ## pdf(snakemake@output[["pdf2"]])
  ## plot(1)
  ## dev.off()
  ## quit()
  ## snakemake <- readRDS("results/supplementary/suppl_bfile/summary.bfile.rds")

  res <- lapply(snakemake@input[["rds"]], readRDS)
  names(res) <- snakemake@params[["k"]]

  saveRDS(res, snakemake@output[["rds"]])

  log <- sapply(res, "[[", "log")
  acc <- sapply(res, "[[", "acc")
  eigvals <- sapply(res, "[[", "eigvals")
  rownames(acc)

  eigvals.a <- t(apply(eigvals, 1, function(fn) {
    e <- read.table(fn[length(fn)])[,1]
    as.numeric(e)
  }))


  if( scenario == "suppl_ukb" ) {
    ukbnames <- list("datas1" = "ukb_nsamples_10K",
                     "datas2" = "ukb_nsamples_20K",
                     "datas3" = "ukb_nsamples_40K",
                     "datas4" = "ukb_nsamples_60K",
                     "datas5" = "ukb_nsamples_80K",
                     "datas6" = "ukb_nsamples_100K",
                     "datas7" = "ukb_nsnps_200K",
                     "datas8" = "ukb_nsnps_400K",
                     "datas9" = "ukb_nsnps_600K",
                     "datas10" = "ukb_nsnps_800K",
                     "datas11" = "ukb_nsnps_1M",
                     "datas12" = "ukb_nsnps_2M",
                     "datas13" = "ukb_nsnps_3M"
                     )
    rownames(acc) <- sapply(rownames(acc), function(i) ukbnames[[i]])
    rownames(log) <- sapply(rownames(log), function(i) ukbnames[[i]])
  }


  if( scenario == "suppl_bfile" ) {
    ukbnames <- list("1000G_nsnps200k" = "1000G_nsnps_200K",
                      "1000G_nsnps400k" = "1000G_nsnps_400K",
                      "1000G_nsnps600k" = "1000G_nsnps_600K",
                      "1000G_nsnps800k" = "1000G_nsnps_800K",
                      "1000G_nsnps1000k" =  "1000G_nsnps_1M",
                      "1000G_nsnps2000k" =  "1000G_nsnps_2M",
                      "ukb_array_nsamples10k" = "ukb_nsamples_10K",
                      "ukb_array_nsamples20k" = "ukb_nsamples_20K",
                      "ukb_array_nsamples40k" = "ukb_nsamples_40K",
                      "ukb_array_nsamples60k" = "ukb_nsamples_60K",
                      "ukb_array_nsamples80k" = "ukb_nsamples_80K",
                      "ukb_array_nsamples100k" = "ukb_nsamples_100K",
                      "ukb_array_nsamples200k" = "ukb_nsamples_200K" )
    rownames(acc) <- sapply(rownames(acc), function(i) ukbnames[[i]])
    rownames(log) <- sapply(rownames(log), function(i) ukbnames[[i]])
  }

  epochs.h <- apply(log, 2, function(l) {
    sapply(l, function(d) {
      d[rownames(d) %in% c("windows=64"), 4] # choose default w=64 only
    })
  })


  epochs.f <- apply(log, 2, function(l) {
    sapply(l, function(d) {
      d[rownames(d) %in% c("windows=64"), 3] # choose default w=64 only
    })
  })

  mev.f <- apply(acc, 2, function(l) {
    sapply(l, function(d) {
      d[rownames(d) %in% c("pcaonef.w64"), 2]
    })
  })

  mev.h <- apply(acc, 2, function(l) {
    sapply(l, function(d) {
      d[7, 2]
      d[rownames(d) %in% c("pcaoneh"), 2]
    })
  })

  pdf(snakemake@output[["pdf1"]], w = 12, h = 9)
  par(mfrow = c(2, 2))
  emax <- max(cbind(epochs.h[1:6,], epochs.f[1:6,])) + 2
  ymin <- min(cbind(mev.h[1:6,], mev.f[1:6,]))
  barplot(epochs.f[1:6, ], beside = T, ylim = c(0, emax), col = factor(rownames(mev.f)[1:6]), ylab = "Epochs", cex.axis = 1.2, cex.lab = 1.2, names.arg = paste0("K=", colnames(epochs.f)), main = expression("PCAone"), cex.main = 2)
  plot_mev_byk(mev.f[1:6, ], ymin)

  barplot(epochs.h[1:6, ], beside = T, ylim = c(0, emax), col = factor(rownames(mev.h)[1:6]), ylab = "Epochs", cex.axis = 1.2, cex.lab = 1.2, names.arg = paste0("K=", colnames(epochs.h)), main = expression("PCAone"[H+Y]), cex.main = 2)
  plot_mev_byk(mev.h[1:6, ], ymin, legend = FALSE)
  dev.off()

  pdf(snakemake@output[["pdf2"]], w = 12, h = 9)
  par(mfrow = c(2, 2))
  emax <- max(cbind(epochs.h[1:6+6,], epochs.f[1:6+6,])) + 2
  ymin <- min(cbind(mev.h[1:6+6,], mev.f[1:6+6,]))
  barplot(epochs.f[1:6+6, ], beside = T, ylim = c(0, emax), col = factor(rownames(mev.f)[1:6+6]), ylab = "Epochs", cex.axis = 1.2, cex.lab = 1.2, names.arg = paste0("K=", colnames(epochs.f)), main = expression("PCAone"), cex.main = 2)
  plot_mev_byk(mev.f[1:6+6, ], ymin)
  barplot(epochs.h[1:6+6, ], beside = T, ylim = c(0, emax), col = factor(rownames(mev.h)[1:6+6]), ylab = "Epochs", cex.axis = 1.2, cex.lab = 1.2, names.arg = paste0("K=", colnames(epochs.h)), main = expression("PCAone"[H+Y]), cex.main = 2)
  plot_mev_byk(mev.h[1:6+6, ], ymin, legend = FALSE)
  dev.off()
}

if (rule == "collect_bfile_summary") {
  log <- lapply(snakemake@input[["log"]], readRDS)
  names(log) <- snakemake@params[["data"]]
  acc <- lapply(snakemake@input[["acc"]], readRDS)
  names(acc) <- snakemake@params[["data"]]
  eigvals <- snakemake@input[["val"]]
  names(eigvals) <- snakemake@params[["data"]]

  res <- list(log = log, acc = acc,eigvals = eigvals)

  saveRDS(res, snakemake@output[["rds"]])

  mw <- min(sapply(lapply(log, rownames), length))
  if( length(mw) == 1)
    quit()

  windows <- rownames(log[[1]])
  windows <- gsub("windows", replacement = "w", x = windows)

  ## res <- readRDS("results/supplementary/suppl_ukb_bysamples/summary.k30.bfile.rds")
  ## res <- readRDS("results/supplementary/suppl_bfile/summary.k40.bfile.rds")
  ## res <- readRDS("results/supplementary/suppl_bfile/summary.k10.bfile.rds")

  K <- snakemake@wildcards[["k"]]

  acc <- res[["acc"]]
  log <- res[["log"]]
  data <- names(acc)

  res.pcaone <- lapply(data, function(n) {
    f <- acc[[n]]["MEV"]
    f <- f[!rownames(f) %in% c("pcaoneh"),]
    o <- list(mev = f, epochs = log[[n]][["epochs.pcaonef"]])
    o
  })
  names(res.pcaone) <- data

  res.halko <- lapply(data, function(n) {
    o <- list(mev = acc[[n]]["pcaoneh","MEV"], epochs = log[[n]][["epochs.pcaoneh"]])
    o
  })
  names(res.halko) <- data

  # Epochs

  pcaone.mev <- lapply(res.pcaone, "[[", "mev")

  (pcaone.epochs <- t(sapply(res.pcaone, "[[", "epochs")))
  (halko.epochs <- t(sapply(res.halko, "[[", "epochs")))

  plot_mev_byw <- function(pcaone.mev, data) {
    l.mev <- pcaone.mev[data]
    ymin <- min(do.call(rbind.data.frame, l.mev))
    plot(1, bty = "l", type = "n", xaxt = "n", col = "transparent", cex.main = 2.0, cex.axis = 1.2, cex.lab = 1.2, xlim = c(1, 6), ylim = c(ymin, 1.0), ylab = paste0("MEV(PC=", K, ")"), xlab = "")
    for (d in data) {
      points(l.mev[[d]], type = "b", lwd = 2, col = which(factor(data) == d))
    }
    axis(1, at = 1:6, labels = windows, cex.axis = 1.2)
    legend("bottomright", legend = data, fill = factor(data), cex = 1.2, bty = "n")
  }

  ## out <- paste0(snakemake@output[["rds"]],".pdf")
  pdf(paste0(snakemake@output[["rds"]], ".pdf"), w = 12, h = 9)
  par(mfrow = c(2, 2))
  barplot(pcaone.epochs[1:6, ], beside = T, col = factor(data[1:6]), ylab = "Epochs", cex.axis = 1.2, cex.lab = 1.2, names.arg = windows)
  plot_mev_byw(pcaone.mev, data[1:6])
  barplot(pcaone.epochs[1:6 + 6, ], beside = T, col = factor(data[1:6+6]), ylab = "Epochs", cex.axis = 1.2, cex.lab = 1.2, names.arg = windows)
  plot_mev_byw(pcaone.mev, data[1:6 + 6])
  dev.off()

  quit()
}

## barplot(t(cbind(halko.epochs[,1],pcaone.epochs[,4])), beside = T)


## do.call(rbind.data.frame, res.pcaone)

# MEV
