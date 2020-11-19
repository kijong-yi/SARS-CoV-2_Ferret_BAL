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

## load seurat object
MPcov <- readRDS("/home/users/kjyi/Projects/covid19/ferret/singlecell/data/Seurat_object_total_cells.Rds")


## Fig1a experiment scheme

## Fig1b histologic finding

## Fig1c marker gene, dot plot

goi <- c("TIFAB","CD14", "CSF3R", "NRP1", 
         "CEBPE","CPA3","KLRK1","CD3E","TRDC","CD8A","CD4",
         "CD79B","IGHG4","MKI67","WFDC2", "HBA2","EPCAM")

Idents(MPcov) <- "Annotation"

p <- DotPlot(
  MPcov,
  features = rev(goi),
  cols = c("lightgrey", "black"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "size",
  scale.min = NA,
  scale.max = NA
)

p$data$id = factor(p$data$id, levels(MPcov$Annotation) %>% rev)

pdf("Fig1c.dotplot.pdf",7.5*1.2,5*1.2)
p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


## Fig1d whole umap

pdf("Fig1d.whole_umap_seurat_clusters.pdf",900/72*0.8, 688/72*0.8) 
DimPlot(MPcov,group.by = "Annotation", label=T)
dev.off()


## Fig1e whole proportion

colorder = c("Ctrl-1","Ctrl-2","Ctrl-3",
             "C2-1","C2-2","C2-3",
             "C5-1","C5-2", "C5-3","C5-4")

x <- table(MPcov$Annotation,MPcov$Sample_name)
x <- x[, colorder]
x3= t(t(x)/rowSums(t(x)))

x4 = as.data.frame(as.table(t(x3)))
colnames(x4) = c("sample","celltype","Freq")
x4$group = x4$sample %>% str_replace("-.*","")
x4$group = factor(x4$group, levels = c("Ctrl","C2","C5"))


write.csv(x4, "Fig1e.proportion_each_cluster.csv")
