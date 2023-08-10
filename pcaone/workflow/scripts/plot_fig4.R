library(cowplot)
library(ggplot2)
library(data.table)
library(caTools)

Sys.setenv(DISPLAY = "localhost:15.0" )
barplot(1:10)

palette(c("#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BEBADA","#FFED6F","#8DD3C7"))
## df = data.table::fread("/emc/zilong/pcaone/manuscript/ukb/all.unrelated.h.manh", h=F, data.table=F)
df = data.table::fread("/maps/projects/alab/scratch/zilong/emc/pcaone/manuscript/ukb/all.unrelated.h.manh", h=F, data.table=F)
npc <- dim(df)[2] - 3
colnames(df) <- c("chr", "snp", "bp", paste0("pc", 1:npc))
chrtable <- table(df$chr)
window <- 1000
mat <- apply(df[,-(1:3)], 2, function(x ) abs(caTools::runmax(x,  k = window )))
groups <- list("Population structure" = c(1, 2, 4),
            "Centromere" = c(3, 5, 8, 9, 10, 11, 13, 14, 19, 20, 22, 29, 32, 34, 35, 37),
            "HLA" =  c(6, 12, 15, 16, 27, 30, 38),
            "Other structure" =  c(7, 17, 18, 23:26, 28, 31, 33, 36, 39, 40),
            "Inversion" =  c(21))

# anno : a vecotor of pc to be colored
plotLegend <- function(anno = c(1, 2, 4)) {
    par(mar = c(0,0,0,0), mgp = c(2, 0.6, 0))
    plot(0, 0, col = "transparent", axes = F, xlab = "", ylab = "")
    cats <- names(groups)
    letext <- c(cats,  paste0("PC", anno))
    legend("bottomleft", legend = letext, col = 1:length(letext) , pch = c(rep(16,length(cats)), rep(NA, 3)), lty = c(rep(NA, length(cats)), rep(1, 3)), adj = 0, lwd = 4, pt.cex = 3.0, cex = 1.2, text.font = 2, horiz = T, xjust = 0, yjust = 0, xpd = T, inset = c(0, 0), bty = "n", x.intersp = 0.1)
}

plotLoadings <- function(window = 1000, anno = c(1, 2, 4 )) {
    if (!is.list(groups)) stop("groups must be a list")
    cats <- names(groups)
    ymax <- c(0, max(mat)) * 1.05
    nchr <- length(chrtable)
    par(mar = c(3.5, 3.5, 4, 1), mgp = c(2, 0.6, 0))
    plot(1:nrow(mat), mat[,1], type = 'l', col = "gray60", lwd = 2, cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5,  ylim = ymax, axes = F, xlab = "Chromosomes", ylab = "SNP loadings", main = "Top 40 PC loadings of 339,582 individuals and 3,841,676 SNPs\n( Cost: 5 hours 44 mins, 6GB, 20 threads )")
    axis(2, cex.lab = 1.5, cex.axis = 1.5)
    ## text(cumsum(as.vector(chrtable)) - chrtable/2, rep(-0.001,22), 1:22, pch = "|", xpd = T, ...)
    text(cumsum(as.vector(chrtable)) - chrtable/2, rep(-ymax[2]/30, nchr), 1:nchr, pch = "|", xpd = T, cex = 1.1)
    apply(mat, 2, function(x) lines(1:nrow(mat), x, lwd = 2, col = "gray60") )
    if (is.vector(anno)) {
        sapply(1:length(anno), function(i) lines(1:nrow(mat), mat[, anno[i]], lwd = 2, col = i+length(cats)) )
    }
    abline(v = cumsum(as.vector(chrtable)), col = "white", lwd = 3)
    w <- apply(mat, 2, which.max)
    for(i in 1:ncol(mat)) {
        c <- which(cats == names(which(sapply(groups, function(x) i %in% x ))))
        if (i %in% c(10, 12, 27)) {
            points(w[i], mat[w[i], i] - 0.001, col = c, pch = 16, cex = 3)
            text(w[i], mat[w[i], i] - 0.001, i)
        } else {
            points(w[i], mat[w[i], i], col = c, pch = 16, cex = 3)
            text(w[i], mat[w[i], i], i)
        }
    }
}

png("ukb-loadings.png", height = 9, width = 18, res = 300, units = "in" )
ggdraw() + draw_plot(plotLegend, 0, 0.88, 1, 0.1 ) + draw_plot(plotLoadings, 0, 0, 1, 0.9)
dev.off()

########################### PCA and LocusZoom #####################

refGenes <- "data/pcaone/refGene.hg38.txt.gz"
d <- read.table("data/pcaone/ukb.all.loadings1-4.chr5.gz", h = T)
colnames(d) <- c("chr", "pos", paste0("pc", 1:4))

plotlocus <- function(){
  locusZoomNoLD <- function(pval, pos, chr, refGenes = NULL, w = 1, plog = FALSE, main = "") {
    xforrs <- 0.03
    regsz <- 1.2
    width <- 22
    rsqplus <- 0.045
    rightylabxplus <- 0.05
    xmargin <- 0.005
    cexaxis <- 1.2
    cexlab <- 1.2
    blue <- "dodgerblue4"

    pval[is.na(pval)] <- 1
    N <- length(pval)
    rsqwithn <- rep(0, N - 1)
    ## palette(brewer.pal(10, "RdBu"))[1:10]
    cols <- c(10, 8, 6, 4, 2, 2)
    getcol <- function(x) {
      cuts <- c(.2, .4, .6, .8, 1)
      cols[6 - sum(as.numeric(x) < cuts)]
    }

    min.pos <- min(pos) / 1e6
    max.pos <- max(pos) / 1e6

    adj <- 0
    par(mar = c(0, 0, 0, 0))
    nf <- layout(matrix(c(1, 2), 2, 1, byrow = TRUE), heights = c(4.0, 2.5), widths = c(9))
    par(mar = c(0, 5.5, 2, 1))

    bgcols <- sapply(rsqwithn, getcol)
    bgcols[w] <- 1
    if (plog) {
      ylim <- c(0, max(10.5, max(-log10(pval))))
      pval <- -log10(pval)
    } else {
      ylim <- c(0, max(pval)*1.05)
    }
    plot(pos / 1e6, pval, col = "black", bg = bgcols, xaxs = "i", xlim = c(min.pos - adj, max.pos + adj), xaxt = "n", xlab = "", ylab = "", ylim = ylim, pch = 21, cex = 2, las = 1, cex.axis = cexaxis, main = main)
    ## axis(2, ylim = c(0, 15), col = "black", label = FALSE, cex.axis = cexaxis)
    ## mtext(2, text = expression(paste("-l", og[10], "[", italic(P), "]")), line = 3, cex = cexlab)
    mtext(2, text = "SNP Loadings of PC1", line = 4, cex = cexlab)
    xstart <- max.pos - (max.pos - min.pos) * 0.15
    scalewidth <- (max.pos - min.pos) * 0.0125
    xx <- c(xstart, xstart, xstart + scalewidth, xstart + scalewidth)
    #  legend("topleft",legend=c(rs),pt.bg=1,col="black",cex=1.6,pch=21)

    ybot <- ylim[2] * 0.65
    ytop <- ylim[2] * 0.95
    ysz <- ytop - ybot
    cuts <- c(.2, .4, .6, .8, 1)
    txtcuts <- c("0.2", "0.4", "0.6", "0.8", "1.0")
    txtposs <- c(0.2, 0.4, 0.6, 0.8, 1.00)
    scalexplus <- (max.pos - min.pos) * 0.06

    dat <- read.table(refGenes, as.is = T, head = T, comment.char = "")
    xx2 <- dat[dat[, "chrom"] == paste("chr", chr, sep = "") & dat[, "cdsStart"] < max.pos * 1e6 & dat[, "cdsEnd"] > min.pos * 1e6, ]

    start <- xx2$txStart # casia. this column is needed in the refgene file
    end <- xx2$txEnd # casia. this column is needed in the refgene file
    nams <- xx2$name2 # casia. this column is needed in the refgene file
    cnts <- xx2$exonCount # casia. this column is needed in the refgene file but you can have all as NA

    ## abline(v = pos, lty = 2)
    ## par(mar = c(5.2, 6.2, -0.1, 6.3) + 0.1)
    par(mar = c(4, 5.5, 1, 1))

    plot(c(0, 0), c(0, 0), type = "n", xlim = c(min.pos - adj, max.pos + adj), ylim = c(-0.8, 0.1), xlab = "", xaxs = "i", yaxt = "n", ylab = "", main = "", cex.lab = 2., cex.axis = cexaxis, tck = -0.02)
    mtext(1, text = paste("Position on chromosome ", chr, " (Mb)", sep = ""), line = 2.5, cex = cexlab)

    ord <- order(start)
    start <- start[ord]
    end <- end[ord]
    exoncnts <- cnts[ord]
    nams <- nams[ord]
    keep <- !duplicated(nams)
    start <- start[keep]
    end <- end[keep]
    exoncnts <- cnts[keep]
    nams <- nams[keep]
    ord <- ord[keep]
    he <- rep(c(0, -0.18, -0.36, -0.54, -0.72), 100)[1:length(nams)] - 0.05
    if (length(start) > 0) {
      segments(start / 1e6, he, end / 1e6, he)
      keep <- !duplicated(nams)
      sapply(1:sum(keep), function(x) {
        text((end[keep][x] + start[keep][x]) / 2e6, he[keep][x] + 0.08, bquote(italic(.(nams[keep][x]))), cex = cexlab - 0.4)
      })
      estart <- as.numeric(unlist(sapply(xx2$exonStarts[ord], function(y) {
        strsplit(y, ",")[[1]]
      }))) / 1e6
      eend <- as.numeric(unlist(sapply(xx2$exonEnds[ord], function(y) {
        strsplit(y, ",")[[1]]
      }))) / 1e6
      rect(estart, rep(he, xx2$exonCount[ord]) - 0.01, eend, rep(he, xx2$exonCount[ord]) + 0.01, col = "black")
    }
  }

  t <- subset(d, pos > 33800000 & pos < 34100000 & chr == 5)
  locusZoomNoLD(t$pc1, t$pos, 5, refGenes = refGenes)
}

## plotlocus()

## png("locuszoom-pc1.png", w = 8,  h = 6, res = 300, units = "in")
## locusZoomNoLD(t$pc1, t$pos, 5, refGenes = refGenes)
## dev.off()

###### pca plot
bigpal=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")


## do <- dd[order(dd$ethnic),]

do <- readRDS("data/pcaone/ukb.all.pca.rds")

catt <- c("Other","Afrian","British","Chinese","Indian","Caribbean","Irish", "Pakistani","Bangladeshi")
ord <- unlist(sapply(catt, function(x,df){which(df$ethnic == x)}, df=do))
dd <- do[ord, ]


bigpal=c("#1B9E77","#FF7F00","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","gray60","#8DD3C7","#FFFFB3","#BEBADA","#FB8072")
p1 <- ggplot(dd, aes(PC1, PC2, colour=ethnic)) + theme_classic() + geom_point(alpha = 0.8)+scale_color_manual(values=bigpal)+theme(legend.position="none", axis.text = element_text(size = 14), axis.title = element_text(size = 16))+labs(x="PC1", y="PC2") + ggtitle("")
p2 <- ggplot(dd, aes(PC1, PC4, colour=ethnic)) + theme_classic() + geom_point(alpha = 0.8)+scale_color_manual(values=bigpal)+theme(legend.position="right", legend.title.align = 0.1, legend.text = element_text(size = 20), legend.title = element_text(size = 26), axis.text = element_text(size = 14), axis.title = element_text(size = 16))+labs(x="PC1", y="PC4", colour="") + ggtitle("") + guides(color = guide_legend(override.aes = list(size = 5) ))

## plot_grid(p1, p2)

w1 <- 8
w2 <- 6
w3 <- 8.5
ww <- w1+w2+w3
png("/Users/zilong/Dropbox/Zll/300-Work/301-PCAone/figures/ukb-pca.png", height = w2, width = ww, res = 300, units = "in" )
ggdraw() + draw_plot(plotlocus, 0, 0, w1/ww, 1) + draw_plot(p1, w1/ww, 0, w2/ww, 1) + draw_plot(p2, (w1+w2)/ww, 0, w3/ww, 1) + draw_plot_label(c("A", "B", "C"), c(0, w1/ww, (w1+w2)/ww), 1)
dev.off()
print("done")
