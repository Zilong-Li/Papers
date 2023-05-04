
source("workflow/scripts/common.R")

meta_label_file <- "/maps/projects/greenland-AUDIT/people/lmj775/Grouop_share/metagenome.species.kraken_numdf.color.tsv"
meta_label <- fread(meta_label_file, h = T, sep = "\t", data.table = F)

meta_label

pca_mat_label_file <- "/maps/projects/greenland-AUDIT/people/rlk420/datahub/metagenome/jgi.sample.ids"
pca_label <- fread(pca_mat_label_file, h = F, data.table = F)[,1]

pca_label


pcaone_f_eigvecs <- "/maps/projects/greenland-AUDIT/people/rlk420/datahub/metagenome/pcaone.f.eigvecs"
pcaone.f <- fread(pcaone_f_eigvecs, h = F, data.table = F)
colnames(pcaone.f) <- paste0("PC", 1:dim(pcaone.f)[2])
pcaone.f$ind <- pca_label


res <- pcaone.f %>% left_join(meta_label, by="ind")
## res[is.na(res)] <- "black"

## p1 <- ggplot(res, aes(PC1, PC2)) +
##   geom_point(aes( colour=`Escherichia coli` )) +
##   scale_color_identity() + # this uses color from given column
##   theme_classic() +
##   labs(title = "Escherichia coli")
## p2 <- ggplot(res, aes(PC1, PC2)) +
##   geom_point(aes( colour=`Prevotella copri` )) +
##   scale_color_identity() + # this uses color from given column
##   theme_classic() +
##   labs(title = "Prevotella copri")
## p3 <- ggplot(res, aes(PC1, PC2)) +
##   geom_point(aes( colour=`Blautia massiliensis` )) +
##   scale_color_identity() + # this uses color from given column
##   theme_classic() +
##   labs(title = "Blautia massiliensis")
## p4 <- ggplot(res, aes(PC1, PC2)) +
##   geom_point(aes( colour=`Bifidobacterium longum` )) +
##   scale_color_identity() + # this uses color from given column
##   theme_classic() +
##   labs(title = "Bifidobacterium longum")
## plot_grid(p1, p2, p3, p4)

png("results/metabinning/Escherichia coli.png", units = "in", res = 200, w = 12, h = 12)
pairs(res[,1:6],  col = res[,"Escherichia coli"], diag.panel = NULL, upper.panel = NULL, main="Escherichia coli")
dev.off()

png("results/metabinning/Prevotella copri.png", units = "in", res = 200, w = 12, h = 12)
pairs(res[,1:6],  col = res[,"Prevotella copri"], diag.panel = NULL, upper.panel = NULL, main="Prevotella copri")
dev.off()

png("results/metabinning/Blautia massiliensis.png", units = "in", res = 200, w = 12, h = 12)
pairs(res[,1:6],  col = res[,"Blautia massiliensis"], diag.panel = NULL, upper.panel = NULL, main="Blautia massiliensis")
dev.off()

png("results/metabinning/Bifidobacterium longum.png", units = "in", res = 200, w = 12, h = 12)
pairs(res[,1:6],  col = res[,"Bifidobacterium longum"], diag.panel = NULL, upper.panel = NULL, main="Bifidobacterium longum")
dev.off()
