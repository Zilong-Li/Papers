## library(caTools)
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)

Sys.setenv(DISPLAY = "localhost:15.0" )
barplot(1:10)

palette(c("#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BEBADA","#FFED6F","#8DD3C7"))
df = data.table::fread("/emc/zilong/pcaone/manuscript/ukb/all.unrelated.h.manh", h=F, data.table=F)
npc <- dim(df)[2] - 3
colnames(df) <- c("chr", "snp", "bp", paste0("pc", 1:npc))
chrtable <- table(df$chr)
window <- 1000
mat <- apply(df[,-(1:3)], 2, function(x ) abs(runmax(x,  k = window )))

groups <- list("Population structure" = c(1, 2, 4),
            "Centromere" = c(3, 5, 8, 9, 10, 11, 13, 14, 19, 20, 22, 29, 32, 34, 35, 37),
            "HLA" =  c(6, 12, 15, 16, 27, 30, 38),
            "Other structure" =  c(7, 17, 18, 23:26, 28, 31, 33, 36, 39, 40),
            "Inversion" =  c(21))
# anno : a vecotor of pc to be colored
plotLegend <- function(anno = c(1, 2, 4 )) {
    par(mar = c(0,0,0,0), mgp = c(2, 0.6, 0))
    plot(0, 0, col = "transparent", axes = F, xlab = "", ylab = "")
    cats <- names(groups)
    letext <- c(cats,  paste0("PC", anno))
    legend("bottomleft", legend = letext, col = 1:length(letext) , pch = c(rep(16,length(cats)), rep(NA, 3)), lty = c(rep(NA, length(cats)), rep(1, 3)), adj = 0, lwd = 4, pt.cex = 3.0, cex = 1.2, text.font = 2, horiz = T, xjust = 0, yjust = 1, xpd = T, inset = c(0, 0.5), bty = "n", x.intersp = 0.15)
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

# setwd("/emc/zilong/pcaone/manuscript")
# d = fread("ukb/all.unrelated.h.pc.ethnic", h=T, data.table=F)

dat1 <- "/emc/zilong/ukb/imputed/merge.maf0.05"
loadings <- "/emc/zilong/pcaone/output/ukb/big/merge.maf0.05/pcaone.H.k40.n20.20208.loadings.gz"
d1 = fread(loadings, select = 1:4,col.names = paste0("PC", 1:4), h=F, data.table=F)
d2 = fread(paste0(dat1, ".bim"), h=F, select = c(1, 4), col.names = c("#chr", "pos"), data.table=F)
write.table(cbind(d2, d1), gzfile("ukb.all.loadings1-4.gz"), quote=F, row.names=F, sep="\t")

dat1 <- "/emc/zilong/ukb/imputed/merge.maf0.05"
pc <- "/emc/zilong/pcaone/output/ukb/big/merge.maf0.05/pcaone.H.k40.n20.20208.eigvecs"
label = "/emc/zilong/pcaone/manuscript/ukb.ethnic.merged"
d1 = fread(pc, select = 1:4, h=F, data.table=F)
d2 = fread(paste0(dat1, ".fam"), h=F, select = 1:2, col.names = c("fid", "iid"), data.table=F)
d3 = fread(label, col.names = c("iid","ethnic"), h=F, data.table=F)
# merge(setDT(d2), setDT(d3), all.x=T, by = c("iid")) # data.table
d2  <- d2 %>% left_join(d3, by="iid")
d  <- cbind(d2$ethnic, d1)
colnames(d)  <- c("ethnic", paste0("PC", 1:(dim(d)[2] - 1)))
# replace values

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
saveRDS(dd, "ukb.all.pca.rds")

bigpal=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")
# nolegend
p1=ggplot(dd, aes(pc1, pc2, colour=ethnic)) + theme_classic() + geom_point()+scale_color_manual(values=bigpal)+theme(legend.position="none")+labs(x="PC1", y="PC2")
p2=ggplot(dd, aes(pc1, pc4, colour=ethnic)) + theme_classic() + geom_point()+scale_color_manual(values=bigpal)+theme(legend.position="right", legend.title.align = 0.1)+labs(x="PC1", y="PC4", colour="Ethnic Group")
p <- plot_grid(p1, p2)
ggsave("pc124.png", p, w=14, h=6)
ggsave("pc1-2.png", p1, w=6, h=6)
