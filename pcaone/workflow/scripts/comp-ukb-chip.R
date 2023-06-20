source("workflow/scripts/common.R")

label <- "data/ukb_info/ukb.ethnic.merged"

bigpal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

## pcaone.gt <- make_ukb_pca("results/ukb_array/pcaone.f.w64.k40.eigvecs", "data/ukb_array/ukb_snp_chip_in_phasing.fam", label)
pcaone.gt <- make_ukb_pca("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/pcaone.f.gt.eigvecs", "/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/ukb_snp_chip.bgen.vars.fam", label)
pcaone.ds <- make_ukb_pca("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/imputed/pcaone.f.ds.eigvecs", "/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/ukb_snp_chip.bgen.vars.fam", label)
pcaone.imp <- make_ukb_pca("/maps/projects/alab/people/rlk420/pcaone/data/ukb/pcaone.f.eigvecs", "/maps/projects/alab/people/rlk420/pcaone/data/ukb/merge.maf0.05.fam", label)

p1 <- ggplot(pcaone.gt, aes(PC1, -PC2, colour = ethnic)) +
  theme_classic() +
  geom_point() +
  scale_color_manual(values = bigpal) +
  theme(legend.position = "none") +
  labs(x = "PC1", y = "PC2", title = "UKB SNP chip genotypes")

p2 <- ggplot(pcaone.imp, aes(PC1, PC2, colour = ethnic)) +
  theme_classic() +
  geom_point() +
  scale_color_manual(values = bigpal) +
  theme(legend.position = "none") +
  labs(x = "PC1", y = "PC2", colour = "Ethnic Group", title = "UKB imputed genotypes")

p3 <- ggplot(pcaone.ds, aes(-PC1, -PC2, colour = ethnic)) +
  theme_classic() +
  geom_point() +
  scale_color_manual(values = bigpal) +
  theme(legend.position = "right", legend.title.align = 0.1) +
  labs(x = "PC1", y = "PC2", colour = "Ethnic Group", title = "UKB imputed dosages")

p <- plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(4, 4, 5))
ggsave("results/supplementary/ukb-array-vs-imputed.png", p, w = 13, h = 4)


library(corrplot)

pcaone.ds <- fread("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/imputed/pcaone.f.ds.eigvecs", h = F, data.table = F)
fam.ds <- fread("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/ukb_snp_chip.bgen.vars.fam", h = F, data.table = F)
ord.ds <- order(fam.ds[, 2])

arnoldi.gt <- fread("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/pcaone.a.gt.eigvecs", h = F, data.table = F)
halko.gt <- fread("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/pcaone.h.gt.eigvecs", h = F, data.table = F)
mev(halko.gt, arnoldi.gt)
## pcaone.gt <- fread("results/ukb_array/pcaone.f.w64.k40.eigvecs", h = F, data.table = F)
## fam.gt <- fread("data/ukb_array/ukb_snp_chip_in_phasing.fam", h = F, data.table = F)
pcaone.gt <- fread("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/pcaone.f.gt.eigvecs", h = F, data.table = F)
fam.gt <- fread("/maps/projects/mono-AUDIT/people/rlk420/datahub/ukbb/snp_chip/ukb_snp_chip.bgen.vars.fam", h = F, data.table = F)
ord.gt <- order(fam.gt[, 2])
pcaone.imp <- fread("/maps/projects/alab/people/rlk420/pcaone/data/ukb/pcaone.f.eigvecs", h = F, data.table = F)
fam.imp <- fread("/maps/projects/alab/people/rlk420/pcaone/data/ukb/pcaone.f.perm.fam", h = F, data.table = F)
ord.imp <- order(fam.imp[, 2])

dim(pcaone.imp)
dim(pcaone.gt)

cor.gtds <- abs(cor(pcaone.gt[ord.gt,], pcaone.ds[ord.ds,]))
rownames(cor.gtds) <- paste0("",  1:nrow(cor.gtds))
colnames(cor.gtds) <- paste0("",  1:nrow(cor.gtds))
cor.gtimp <- abs(cor(pcaone.gt[ord.gt, ], pcaone.imp[ord.imp, ]))
rownames(cor.gtimp) <- paste0("",  1:nrow(cor.gtimp))
colnames(cor.gtimp) <- paste0("",  1:nrow(cor.gtimp))

cor.dsimp <- abs(cor(pcaone.ds[ord.ds,], pcaone.imp[ord.imp,]))
rownames(cor.dsimp) <- paste0("",  1:nrow(cor.dsimp))
colnames(cor.dsimp) <- paste0("",  1:nrow(cor.dsimp))

png("results/supplementary/ukb-pc-corr.png", units = "in", res = 200, w = 12, h = 4)
par(mfrow = c(1, 3))
corrplot(cor.gtimp, method = "color", is.corr = F, tl.col = "black", tl.cex = 0.8, title = "PC Correlation", mar = c(0, 3, 2, 0))
mtext(text = "SNP Chip genotypes (M=498,444)", side = 2, las = 0, line = 2, cex = 1.2)
mtext(text = "Imputed genotypes (M=6,133,304)", side = 1, las = 1, line = 4, cex = 1.2)
corrplot(cor.gtds, method = "color", is.corr = F, tl.col = "black", tl.cex = 0.8, title = "PC Correlation", mar = c(0, 2, 2, 0))
mtext(text = "SNP Chip genotypes (M=498,444)", side = 2, las = 0, line = 3, cex = 1.2)
mtext(text = "Imputed dosages (M=6,133,304)", side = 1, las = 1, line = 4, cex = 1.2)
corrplot(cor.dsimp, method = "color", is.corr = F, tl.col = "black", tl.cex = 0.8, title = "PC Correlation", mar = c(0, 2, 2, 0))
mtext(text = "Imputed dosages (M=6,133,304)", side = 2, las = 0, line = 3, cex = 1.2)
mtext(text = "Imputed genotypes (M=6,133,304)", side = 1, las = 1, line = 4, cex = 1.2)
dev.off()

