
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
library(gridGraphics)
library(gridExtra)


norm_vec <- function(x) {
  sqrt(sum(x^2))
}

mev <- function(Y, truth) {
  if (!isTRUE(all.equal(dim(Y), dim(truth)))) {
    return(NA)
  }
  mean(apply(truth, 2, function(y) {
    norm_vec(t(Y) %*% y)
  }))
}

mev_byk <- function(X, Y) {
  stopifnot(all.equal(dim(X), dim(Y)))
  o <- sapply(2:ncol(X), function(n) {
    mev(X[, 1:n], Y[, 1:n])
  })
}

flipmat <- function(X, Y) {
  switched <- colSums(abs(X - Y)) > 2 * colSums(abs(X + Y))
  if (any(switched)) {
    Y[, switched] <- -Y[, switched]
  }
  Y
}

rmse <- function(X, Y) {
  Y <- flipmat(X, Y)
  sqrt(mean((X - Y)^2))
}

## meanSSE measurement
meanSSE <- function(y, truth = x, returnSSE = TRUE) {
  y2 <- y
  K <- ncol(truth)
  acc <- function(x, y) {
    sum((x - y)^2)
  }
  for (i in 1:K) {
    w <- which.min(sapply(1:K, function(j) acc(truth[, i], y[, j])))
    w2 <- which.min(sapply(1:K, function(j) acc(truth[, i], -y[, j])))
    if (acc(truth[, i], -y[, w2]) < acc(truth[, i], y[, w])) {
      y2[, i] <- -y[, w2]
    } else {
      y2[, i] <- y[, w]
    }
  }
  if (returnSSE) {
    return(acc(truth, y2) / K)
  }
  return(y2)
}

minSSE <- function(y, truth = x) {
  if (!isTRUE(all.equal(dim(y), dim(truth)))) {
    return(NA)
  }
  K <- ncol(truth)
  res <- rep(0, K)
  acc <- function(x, y) {
    sum((x - y)^2)
  }
  for (i in 1:K) {
    w <- which.min(sapply(1:K, function(j) acc(truth[, i], y[, j])))
    w2 <- which.min(sapply(1:K, function(j) acc(truth[, i], -y[, j])))
    s1 <- acc(truth[, i], y[, w])
    s2 <- acc(truth[, i], -y[, w2])
    if (s2 < s1) {
      res[i] <- s2
    } else {
      res[i] <- s1
    }
  }
  res <- ifelse(sum(res) <= 1e-10, 1e-10, sum(res))
  res
}

# matrix x and y
# return a columns of y for min dist to x
matchMin <- function(x, y, returnRMSE = FALSE) {
  if (!isTRUE(all.equal(dim(x), dim(y)))) {
    return(NA)
  }
  y2 <- y
  K <- ncol(x)
  acc <- function(x, y) {
    sqrt(sum((x - y)^2))
  }
  for (i in 1:K) {
    w <- which.min(sapply(1:K, function(j) acc(x[, i], y[, j])))
    w2 <- which.min(sapply(1:K, function(j) acc(x[, i], -y[, j])))
    if (acc(x[, i], -y[, w2]) < acc(x[, i], y[, w])) {
      y2[, i] <- -y[, w2]
    } else {
      y2[, i] <- y[, w]
    }
  }
  if (returnRMSE) {
    return(acc(x, y2))
  }
  return(y2)
}

flip <- function(x, y) {
  if (!isTRUE(all.equal(dim(x), dim(y)))) {
    return(NULL)
  }
  if (length(dim(x)) == 1) {
    if (sum(abs(x - y)) > sum(abs(y + x))) {
      return(-x)
    } else {
      return(x)
    }
  }
  w <- colSums(abs(x - y)) > colSums(abs(x + y))
  x[, w] <- -x[, w]
  return(x)
}

## parse the output of /usr/bin/time -v
gettimes <- function(ss) {
  sapply(strsplit(ss, ":"), function(s) {
    s <- as.numeric(s)
    n <- length(s)
    sum(sapply(1:n, function(i) {
      s[i] * 60^(n - i)
    }))
  })
}

# check file empty and Exit status: 0
checklog <- function(lines) {
  if (identical(lines, character(0))) {
    return(NULL)
  }
  if (length(grep("Exit status: 0", lines)) == 0) {
    return(NULL)
  }
  lines
}

logram <- function(l) {
  sapply(l, function(lines) {
    lines <- checklog(lines)
    if (is.null(lines)) {
      return(NA)
    }
    a <- strsplit(lines[grep("Maximum", lines)], " ")[[1]]
    as.numeric(a[length(a)]) # in kbytes
  })
}

logtimes <- function(l) {
  sapply(l, function(lines) {
    lines <- checklog(lines)
    if (is.null(lines)) {
      return(NA)
    }
    a <- strsplit(lines[grep("Elapsed", lines)], " ")[[1]]
    gettimes(a[length(a)]) # in seconds
  })
}

logepochs <- function(l) {
  sapply(l, function(lines) {
    lines <- checklog(lines)
    if (is.null(lines)) {
      return(NA)
    }
    a <- strsplit(lines[grep("stops", lines)], " ")[[1]]
    as.numeric(gsub("epoch=", "", a[length(a)]))
  })
}

gnutime <- function(dl) {
  sapply(dl, function(d) {
    sum(gettimes(d[, 1]))
  })
}

gunram <- function(dl) {
  sapply(dl, function(d) {
    max(d[, 2]) / 1024 # MB units
  })
}

plotacc <- function(mat, xlabels, ...) {
  par(mar = c(3.5, 3.5, 2, 1), mgp = c(2, 0.6, 0))
  plot(1, bty = "l", type = "n", xaxt = "n", col = "transparent", cex.main = 2.0, cex.axis = 1.2, cex.lab = 1.2, xlim = c(1, nrow(mat)), ...)
  ## abline(h = 1, col = "gray", lwd = 2, lty = 2)
  for (i in seq_len(ncol(mat))) {
    points(mat[, i], type = "b", col = i, lwd = 2)
  }
  axis(1, at = seq_len(nrow(mat)), labels = xlabels, cex.axis = 1.2)
}

mplot <- function(mat, x, ylim = range(mat), ...) {
  ## plot(x, mat[,1], bty = "l",  col = "transparent", ylim = range(mat), ...)
  plot(x, mat[, 1], bty = "l", col = "transparent", ylim = ylim, ...)
  for (i in seq_len(ncol(mat))) {
    points(x, mat[, i], type = "b", col = i, lwd = 2)
  }
}

## wong <- c("black", "#f0e442", "#e69f00", "#d55e00","#009e73","#56b4e9","#cc79a7", "#0072b2")

colorAnd <- function(x = "colorblind6b") {
  if (x == "line") {
    palette(c("#b2182b", "#2166ac", "#4DAF4A", "#FF7F00", "#F781BF", "#984EA3"))
  } else if (x == "rasmus") {
    palette(c("mistyrose", "lavender", "lightyellow", "lightblue", "lightgreen", "seashell", "lightcyan"))
  } else if (x == "colorblind6b") {
    palette(c("#b2182b", "#ef8a62", "#fddbc7", "#d1e5f0", "#67a9cf", "#2166ac"))
  } else if (x == "anders") {
    palette(c("darkgreen", "#00A600FF", "yellow", "#E9BD3AFF", "orange", "coral4", "red4", "black"))
  } else if (x == "large") {
    palette(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"))
  } else if (x == "ffs") {
    palette(c("#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#444444"))
  }
}

wong <- c("#000000", "#f0e442", "#e69f00", "#d55e00", "#cc79a7", "#56b4e9", "#009e73")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

PROGS <- data.frame("flashpca" = c("#999999", "FlashPCA2"), "plink2" = c("#009E73", "Plink2(FastPCA)"), "terapca" = c("#E69F00", "TeraPCA"), "propca" = c("#D55E00", "ProPCA"), "pcaone.a" = c("#CC79A7", "PCAone_Arnoldi"), "pcaone.h" = c("#0072B2", "PCAone_H+Y"), "pcaone.f" = c("#56B4E9", "PCAone"))
rownames(PROGS) <- c("col", "name")

fancyname <- list("pcaone.h" = expression(bold("PCAone"[H + Y])), "pcaone.f" = expression(bold("PCAone")))

fancyname <- function(c, bold = FALSE) {
  nb <- list(
    "pcaone.h" = expression(paste(bold("PCAone"[H + Y]))),
    "pcaone.f" = expression(paste(bold("PCAone"))),
    "pcaone.a" = expression(paste(bold("PCAone"[Arnoldi]))),
    "terapca" = expression(paste(bold("TeraPCA"))),
    "propca" = expression(paste(bold("ProPCA*"))),
    "flashpca" = expression(paste(bold("FlashPCA2"))),
    "plink2" = expression(paste(bold("Plink2(FastPCA)*")))
  )
  np <- list(
    "pcaone.h" = expression(paste("PCAone"[H + Y])),
    "pcaone.f" = expression(paste("PCAone")),
    "pcaone.a" = expression(paste("PCAone"[Arnoldi])),
    "terapca" = expression(paste("TeraPCA")),
    "propca" = expression(paste("ProPCA")),
    "flashpca" = expression(paste("FlashPCA2")),
    "plink2" = expression(paste("Plink2(FastPCA)"))
  )
  if (bold) {
    return(sapply(c, function(i) nb[[i]]))
  } else {
    return(sapply(c, function(i) np[[i]]))
  }
}

fancyname2 <- function(c, bold = FALSE) {
  nb <- list(
    "pcaone.h" = expression(paste(bold("PCAone"[H + Y]))),
    "pcaone.f" = expression(paste(bold("PCAone"))),
    "pcaone.a" = expression(paste(bold("PCAone"[Arnoldi]))),
    "terapca" = expression(paste(bold("TeraPCA"))),
    "propca" = expression(paste(bold("ProPCA*"))),
    "flashpca" = expression(paste(bold("FlashPCA2"))),
    "plink2" = expression(paste(bold("Plink2*")))
  )
  np <- list(
    "pcaone.h" = expression(paste("PCAone"[H + Y])),
    "pcaone.f" = expression(paste("PCAone")),
    "pcaone.a" = expression(paste("PCAone"[Arnoldi])),
    "terapca" = expression(paste("TeraPCA")),
    "propca" = expression(paste("ProPCA")),
    "flashpca" = expression(paste("FlashPCA2")),
    "plink2" = expression(paste("Plink2"))
  )
  if (bold) {
    return(sapply(c, function(i) nb[[i]]))
  } else {
    return(sapply(c, function(i) np[[i]]))
  }
}

readfastpca <- function(fn) {
  if (file.size(fn) == 0L) {
    return(NULL)
  } else {
    a <- fread(cmd = paste("tail +2", fn))
    return(as.matrix(a[, 1:(ncol(a) - 1)][, -1]))
  }
}

readflashpca <- function(fn) {
  if (file.size(fn) == 0L) {
    return(NULL)
  } else {
    a <- fread(fn, header = T)
    return(as.matrix(a[, -c(1, 2)]))
  }
}

readplink1 <- function(fn) {
  if (file.size(fn) == 0L) {
    return(NULL)
  } else {
    return(as.matrix(fread(fn)[, -c(1, 2)]))
  }
}

readplink2 <- function(fn) {
  if (file.size(fn) == 0L) {
    return(NULL)
  } else {
    a <- fread(fn, header = T)
    return(as.matrix(a[, -c(1, 2)]))
  }
}

## make eigenvalues barplot

ploteigvals <- function(ff, name = NULL) {
  # input eigvals files
  n <- length(ff)
  eigvals <- lapply(ff, function(f) {
    v <- read.table(f)[, 1]
    v
  })
  if (is.null(name)) {
    name <- names(eigvals)
  }
  nc <- ceiling(n / 2)
  par(mfrow = c(2, nc))
  for (i in seq_len(length(eigvals))) {
    barplot(eigvals[[i]], names.arg = seq_len(length(eigvals[[i]])), main = name[i], xlab = "PCs", ylab = "Eigenvalues", cex.lab = 1.5)
  }
}


make_ukb_pca <- function(pc, fam, label){
  d1 = fread(pc, select = 1:4, h=F, data.table=F)
  d2 = fread(fam, h=F, select = 1:2, col.names = c("fid", "iid"), data.table=F)
  d3 = fread(label, col.names = c("iid","ethnic"), h=F, data.table=F)
  d2 <- d2 %>% left_join(d3, by="iid")
  d  <- cbind(d2$ethnic, d1)
  colnames(d)  <- c("ethnic", paste0("PC", 1:(dim(d)[2] - 1)))
  dd <- d %>% mutate(ethnic = case_when(ethnic == 1 ~ 'Other Asian backgroud',
                                        ethnic == 1001 ~ 'British',
                                        ethnic == 2001 ~ 'White and Black Caribbean',
                                        ethnic == 3001 ~ 'Indian',
                                        ethnic == 4001 ~ 'Caribbean',
                                        ethnic == 2 ~ 'Other mixed background',
                                        ethnic == 1002 ~ 'Irish',
                                        ethnic == 2002 ~ 'White and Black Afrian',
                                        ethnic == 3002 ~ 'Pakistani',
                                        ethnic == 4002 ~ 'Afrian',
                                        ethnic == 3 ~ 'Other Asian backgroud',
                                        ethnic == 1003 ~ 'Other White backgroud',
                                        ethnic == 2003 ~ 'White and Asian',
                                        ethnic == 3003 ~ 'Bangladeshi',
                                        ethnic == 4003 ~ 'Other Black backgroud',
                                        ethnic == 4 ~ 'Other Black backgroud',
                                        ethnic == 2004 ~ 'Other mixed background',
                                        ethnic == 3004 ~ 'Other Asian backgroud',
                                        ethnic == 5 ~ 'Chinese',
                                        ethnic == 6 ~ 'Other ethnic',
                                        TRUE ~ as.character(NA)))
  dd <- d %>% mutate(ethnic = case_when(ethnic == 1 ~ 'Other',
                                        ethnic == 1001 ~ 'British',
                                        ethnic == 2001 ~ 'Other',
                                        ethnic == 3001 ~ 'Indian',
                                        ethnic == 4001 ~ 'Caribbean',
                                        ethnic == 2 ~ 'Other',
                                        ethnic == 1002 ~ 'Irish',
                                        ethnic == 2002 ~ 'Other',
                                        ethnic == 3002 ~ 'Pakistani',
                                        ethnic == 4002 ~ 'Afrian',
                                        ethnic == 3 ~ 'Other',
                                        ethnic == 1003 ~ 'Other',
                                        ethnic == 2003 ~ 'Other',
                                        ethnic == 3003 ~ 'Bangladeshi',
                                        ethnic == 4003 ~ 'Other',
                                        ethnic == 4 ~ 'Other',
                                        ethnic == 2004 ~ 'Other',
                                        ethnic == 3004 ~ 'Other',
                                        ethnic == 5 ~ 'Chinese',
                                        ethnic == 6 ~ 'Other',
                                        TRUE ~ "Other"))
  dd
}
