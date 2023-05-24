source("workflow/scripts/common.R")

label <- "data/ukb_info/ukb.ethnic.merged"

bigpal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

pcaone.gt <- make_ukb_pca("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/pcaone.f.gt.eigvecs", "/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/ukb_snp_chip.bgen.vars.fam", label)

pcaone.imp <- make_ukb_pca("/maps/projects/alab/people/rlk420/pcaone/data/ukb/pcaone.f.eigvecs", "/maps/projects/alab/people/rlk420/pcaone/data/ukb/merge.maf0.05.fam", label)


p1 <- ggplot(pcaone.gt, aes(PC1, -PC2, colour = ethnic)) +
  theme_classic() +
  geom_point() +
  scale_color_manual(values = bigpal) +
  theme(legend.position = "none") +
  labs(x = "PC1", y = "PC2", title = "PCA Based on UKBB Array Genotypes")

p2 <- ggplot(pcaone.imp, aes(PC1, PC2, colour = ethnic)) +
  theme_classic() +
  geom_point() +
  scale_color_manual(values = bigpal) +
  theme(legend.position = "right", legend.title.align = 0.1) +
  labs(x = "PC1", y = "PC2", colour = "Ethnic Group", title = "PCA Based on All UKBB imputed Genotypes")

p <- plot_grid(p1, p2, rel_widths = c(6, 7))
ggsave("results/supplementary/ukb-array-vs-imputed.png", p, w = 13, h = 6)


library(corrplot)

pcaone.imp <- fread("/maps/projects/alab/people/rlk420/pcaone/data/ukb/pcaone.f.eigvecs", h = F, data.table = F)
pcaone.gt <- fread("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/pcaone.f.gt.eigvecs", h = F, data.table = F)
pcaone.ds <- fread("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/pcaone.f.ds.eigvecs", h = F, data.table = F)

cor.gtds <- abs(cor(pcaone.gt, pcaone.ds))

fam.imp <- fread("/maps/projects/alab/people/rlk420/pcaone/data/ukb/pcaone.f.perm.fam", h = F, data.table = F)
fam.gt <- fread("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/pcaone.f.gt.perm.fam", h = F, data.table = F)

ord.gt <- order(fam.gt[, 2])
ord.imp <- order(fam.imp[, 2])


png("results/supplementary/ukb-pc-corr.png", units = "in", res = 200, w = 12, h = 6)
par(mfrow = c(1, 2))
rownames(cor.gtds) <- paste0("PC",  1:nrow(cor.gtds))
colnames(cor.gtds) <- paste0("PC",  1:nrow(cor.gtds))
corrplot(cor.gtds, method = "color", is.corr = F, tl.col = "black", tl.cex = 0.8, title = "Array genotypes VS imputed dosages", mar = c(0, 0, 2, 0))
cor.gtimp <- abs(cor(pcaone.gt[ord.gt, ], pcaone.imp[ord.imp, ]))
rownames(cor.gtimp) <- paste0("PC",  1:nrow(cor.gtimp))
colnames(cor.gtimp) <- paste0("PC",  1:nrow(cor.gtimp))
corrplot(cor.gtimp, method = "color", is.corr = F, tl.col = "black", tl.cex = 0.8, title = "Array genotypes VS all imputed genotypes", mar = c(0, 0, 2, 0))
dev.off()
