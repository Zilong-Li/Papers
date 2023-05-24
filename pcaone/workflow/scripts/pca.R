
saveRDS(snakemake, snakemake@output[["rds"]] )
## q()

snakemake@source("common.R")

## snakemake <- readRDS("results/supplementary/suppl_rnaseq/gtex_rnaseq.pca.rds")
## source("workflow/scripts/common.R")

rule <- snakemake@rule
scenario <- snakemake@config$scenario

vec <- snakemake@input[["vec"]]
val <- snakemake@input[["val"]]
pcs <- as.matrix(fread(vec[length(vec)], data.table = F))
eigs <- read.table(val[length(val)])[,1]

## readH <- ifelse(scenario == "suppl_rnaseq", TRUE, FALSE)

if(scenario == "suppl_rnaseq") {
  label <- read.table(snakemake@params$label, h = T)
  groups <- label[,ncol(label)]
  pdf(paste0(snakemake@output[["rds"]], ".pdf"), w = 12, h = 6)
  colorAnd("large")
  par(mfrow = c(1, 2))
  plot(pcs[,1:2], col=factor(groups), xlab = "PC1", ylab = "PC2")
  plot(eigs, ylab = "Eigen Values", xlab = "PC")
  dev.off()
}
