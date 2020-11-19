## install packages
install.packages('Seurat')
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("cowplot")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")

## load packages
library(Seurat)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(dplyr)
library(ggplot2)
library(tidyr)


## Pre-setting & Sub clustering analysis (NK & T cells)
MPcov <- readRDS("/home/users/kjyi/Projects/covid19/ferret/singlecell/data/Seurat_object_total_cells.Rds")

Idents(MPcov) <- "Annotation"

cov.NK <- subset(MPcov, idents = c("NK cell"))
cov.NK <- RunPCA(cov.NK)
ElbowPlot(cov.NK)
cov.NK <- FindNeighbors(cov.NK, dims = 1:5, force.recalc = T)
cov.NK <- FindClusters(cov.NK, resolution = 0.3)
cov.NK <- RunUMAP(object = cov.NK, dims = 1:5)


cov.CD8 <- subset(MPcov, idents = c("CD8 T cell"))
cov.CD8 <- RunPCA(cov.CD8)
ElbowPlot(cov.CD8)
cov.CD8 <- FindNeighbors(cov.CD8, dims = 1:5, force.recalc = T)
cov.CD8 <- FindClusters(cov.CD8, resolution = 0.2)
cov.CD8 <- RunUMAP(object = cov.CD8, dims = 1:5)


## Fig2a whole umap ; NK

pdf("Fig2a.NK_umap_seurat_clusters.pdf",900/72*0.8, 688/72*0.8) 
DimPlot(cov.NK,group.by = "seurat_clusters", label=T)
dev.off()


## Fig2b bar plot ; NK

pdf("Fig2b.NK_stacked_bar_composition_celltype.pdf",3.8,5)
x = table(cov.NK$seurat_clusters,cov.NK$Experimental_group)
for(i in 1:ncol(x)){
  x[,i] = x[,i] /sum(x[,i])
}
colSums(x)
x = x[nrow(x):1,c("Ct","C2","C5")]
colnames(x) = c("dpi 0","dpi 2","dpi 5")

x <- x[nrow(x):1,]

par(mar=c(11,4,1,1))
barplot(x,las=2)
dev.off()


## Fig2c Vln plot ; NK

pdf("Fig2c.vln_plot_NK.pdf",3.8,5)
VlnPlot(cov.NK, features = "STAT1", pt.size = 0)
VlnPlot(cov.NK, features = "OAS1", pt.size = 0)
VlnPlot(cov.NK, features = "ISG15", pt.size = 0)
VlnPlot(cov.NK, features = "GZMB", pt.size = 0)
VlnPlot(cov.NK, features = "GZMK", pt.size = 0)
VlnPlot(cov.NK, features = "PRF1", pt.size = 0)
def.off()


## Fig2d whole umap ; CD8

pdf("Fig2d.CD8_umap_seurat_clusters.pdf",900/72*0.8, 688/72*0.8)
DimPlot(cov.CD8,label=T,group.by="seurat_clusters", pt.size = 2) + coord_fixed()
dev.off()


## Fig2e Vln plot resident; CD8
pdf("Fig2e.vln_plot_resident_CD8.pdf",3.8,5)
VlnPlot(cov.CD8, features = "CD69", pt.size = 0)
VlnPlot(cov.CD8, features = "S1PR1", pt.size = 0)
VlnPlot(cov.CD8, features = "ITGAE", pt.size = 0)
dev.off()


## Fig2f Vln plot activation; CD8
pdf("Fig2e.vln_plot_activation_CD8.pdf",3.8,5)
VlnPlot(cov.CD8, features = "OAS1", pt.size = 0)
VlnPlot(cov.CD8, features = "ISG15", pt.size = 0)
VlnPlot(cov.CD8, features = "IFNG", pt.size = 0)
VlnPlot(cov.CD8, features = "GZMB", pt.size = 0)
VlnPlot(cov.CD8, features = "PRF1", pt.size = 0)
dev.off()


## Fig2g batch_umap ; CD8

pdf("Fig2g.condition_umap_CD8.pdf",9,15)
par(mfrow=c(1,3), mar = c(1,1,3,1),pty="s", oma=c(0,0,0,0))
for(gr in c("Ct","C2","C5")){
  plot(Embeddings(cov.CD8,reduction="umap")[cov.CD8$Experimental_group!=gr,],col="#00000020",
       pch=16,cex=0.8,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
  points(Embeddings(cov.CD8,reduction="umap")[cov.CD8$Experimental_group==gr,],col="#168C4060",pch=16,cex=1)
  title(c("Ct"="dpi 0","C2"="dpi 2", "C5"="dpi 5")[gr])
}
dev.off()


## Fig2h feature_plot ; cd8

pdf("Fig2h.featureplot.pdf",12,8*0.8)
goi <- c("OAS1", "ISG15")#,
p = FeaturePlot(cov.CD8, features = goi,reduction="umap",ncol=2,combine = FALSE,cols = c("grey80", "#0A64A4"))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p,ncol = 2)
dev.off()