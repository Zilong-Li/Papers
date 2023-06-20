
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

if(scenario == "suppl_rnaseq") {
  label <- read.table(snakemake@params$label, h = T)
  groups <- label[,ncol(label)]
  pdf(paste0(snakemake@output[["rds"]], ".pdf"), w = 12, h = 6)
  colorAnd("large")
  par(mfrow = c(1, 2))
  plot(pcs[,1:2], col=factor(groups), xlab = "PC1", ylab = "PC2", main = "Samples colored by tissues",cex = 1.3, cex.lab = 1.3, cex.main = 1.5)
  plot(eigs, ylab = "Eigenvalues", xlab = "PC", main = "Distribution of eigenvalues", cex = 1.3, cex.lab = 1.3, cex.main = 1.5)
  dev.off()
}

if(scenario == "suppl_meta") {
  h <- read.table("/maps/projects/greenland-AUDIT/people/rlk420/datahub/metagenome/metahit.depth.txt.gz",head=T,nr=1000)
  r <- read.table("/maps/projects/greenland-AUDIT/people/rlk420/datahub/metagenome/metahit.depth.txt.gz",head=T,colC=sapply(h,class))
  #UID
  speciesID <-sapply( strsplit(r[,1],"|",fixed=T),function(x) gsub("\\].*","",gsub(".*\\[","",x[5])) )
  speciesDepth <- tapply(r$totalAvgDepth,speciesID,sum)
  depthI <- colSums(r[,5:ncol(r)])
  speciesRelAbundanceI <- tapply(1:nrow(r),speciesID,function(x) colSums(r[x,5:ncol(r)])/depthI)
  speciesRelAbundanceI <- do.call(rbind,speciesRelAbundanceI)
  o <- order(speciesDepth,decreasing=T)
  coll <- function(x)
    cut(x,quantile(x,0:5/5),include.lowest=T)
  palette(c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac"))
  i <- 1
  w <- o[i]
  spe <- r[which(speciesID==names(speciesDepth[w]))[1],1]
  pdf(paste0(snakemake@output[["rds"]], ".pdf"), w = 12, h = 6)
  par(mfrow=c(1,2))
  plot(pcs[,c(1, 2)],col=coll(speciesRelAbundanceI[w,]),lwd=2,main=names(speciesDepth[w]), xlab="PC1", ylab="PC2", cex=1.3, cex.lab=1.3, cex.main=1.5)
  ## plot(pcs[,c(2, 3)],col=coll(speciesRelAbundanceI[w,]),lwd=2,main=names(speciesDepth[w]), xlab="PC2", ylab="PC3", cex=1.3, cex.lab=1.3, cex.main=1.5)
  ## plot(pcs[,c(2, 4)],col=coll(speciesRelAbundanceI[w,]),lwd=2,main=names(speciesDepth[w]), xlab="PC2", ylab="PC4", cex=1.3, cex.lab=1.3, cex.main=1.5)
  ## plot(pcs[,c(2, 5)],col=coll(speciesRelAbundanceI[w,]),lwd=2,main=names(speciesDepth[w]), xlab="PC2", ylab="PC5", cex=1.3, cex.lab=1.3, cex.main=1.5)
  ## plot(pcs[,c(2, 6)],col=coll(speciesRelAbundanceI[w,]),lwd=2,main=names(speciesDepth[w]), xlab="PC2", ylab="PC6", cex=1.3, cex.lab=1.3, cex.main=1.5)
  plot(eigs, ylab = "Eigenvalues", xlab = "PC", main = "Distribution of eigenvalues", cex=1.3, cex.lab=1.3, cex.main=1.5)
  dev.off()

}

if(1 == 2) {

  K <- 40
  eigvals.asian <- 2 * read.table("/maps/projects/alab/people/rlk420/pcaone/results/asian2/pcaone.full.eigvals")[1:K,1]
  eigvals.1KGP <- 2 * read.table("/maps/projects/alab/people/rlk420/pcaone/results/1000G/pcaone.full.eigvals")[1:K,1]
  eigvals.hgdp <- 2 * read.table("/maps/projects/alab/people/rlk420/pcaone/results/HGDP/pcaone.full.eigvals")[1:K,1]
  eigvals.ukb <- 2 * read.table("/maps/projects/alab/scratch/zilong/emc/pcaone/paper/ukb/pcaone.a.k40.eigvals")[1:K,1]
  eigvals.brain <- read.table("/maps/projects/alab/scratch/zilong/emc/pcaone/paper/scrnas/brain_v1/pcaone.a.k40.eigvals")[1:K,1]
  eigvals.gtex <- read.table("/maps/projects/alab/people/rlk420/pcaone/results/gtex_rnaseq/pcaone.a.k40.binary.eigvals")[1:K,1]

  pdf("results/supplementary/eigenvalues.pdf", w = 12, h = 9)
  par(mfrow = c(2,3))
  barplot(eigvals.asian, names.arg = 1:K, cex.lab = 1.5,ylab = "Eigenvalues",xlab = "PC", main = "1000 Genomes East Asian Individuals")
  barplot(eigvals.1KGP, names.arg = 1:K, cex.lab = 1.5,ylab = "Eigenvalues",xlab = "PC",main = "1000 Genomes All Individuals")
  barplot(eigvals.hgdp, names.arg = 1:K, cex.lab = 1.5,ylab = "Eigenvalues",xlab = "PC",main = "HGDP All Individuals")
  barplot(eigvals.ukb, names.arg = 1:K, cex.lab = 1.5,ylab = "Eigenvalues",xlab = "PC",main = "UK Biobank All Individuals")
  barplot(eigvals.brain, names.arg = 1:K, cex.lab = 1.5,ylab = "Eigenvalues",xlab = "PC",main = "1.3 Million Mice Brain Single Cells")
  barplot(eigvals.gtex, names.arg = 1:K, cex.lab = 1.5,ylab = "Eigenvalues",xlab = "PC",main = "Gene Expression from GTEx Bulk RNA-seq")
  dev.off()

}
