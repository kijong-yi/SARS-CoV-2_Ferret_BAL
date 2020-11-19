## install packages
install.packages('Seurat')
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("cowplot")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("piano")

install.packages("msigdbr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("monocle")

## load packages
library(Seurat)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(dplyr)
library(ggplot2)
library(tidyr)
library(piano)
library(msigdbr)
library(monocle)


## load seurat object
MP <- readRDS("/home/users/kjyi/Projects/covid19/ferret/singlecell/data/Seurat_object_monocyte_macrophage_DC.Rds")


## define cell cluster name 
MP$Annotation <- as.character(MP$Annotation)
MP$Annotation_2 <- as.character(MP$Annotation)

MP$Annotation[MP$Annotation_2 == "APOE pos FABP4 high tissue M2"] <- "APOE+ tissue macrophage"
MP$Annotation[MP$Annotation_2 == "SPP1 high fibrogenic M2"] <- "SPP1hi CHIT1int profibrogenic M2"
MP$Annotation[MP$Annotation_2 == "Transitional M1"] <- "Weakly activated M1"
MP$Annotation[MP$Annotation_2 == "Interferon stimulated M1"] <- "Highly activated M1"
MP$Annotation[MP$Annotation_2 == "Infiltrating macrophage"] <- "Monocyte-derived infiltrating macrophage"
MP$Annotation[MP$Annotation_2 == "APOE pos FABP4 high tissue M2"] <- "APOE+ tissue macrophage"


## Fig5 b,c

### calculate pseudotime

SO.traj1 <- MP[,MP$Annotation %in% c("Monocyte-derived infiltrating macrophage", "Weakly activated M1", "Highly activated M1")]

### 
# ##### make input data 
pd <- new('AnnotatedDataFrame', data = SO.traj1@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(SO.traj1), row.names = row.names(SO.traj1)))

cds <- newCellDataSet(as(SO.traj1@assays$SCT@data, "matrix"), phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

cds$clusters <- SO.traj1$Annotation

### 
# ##### Estimate size factors and dispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# ##### Filtering low-quality cells 
set.seed(123)
cds <- detectGenes(cds, min_expr = 0.1)

# ##### Select genes use Seurat::FindAllMarkers function
markers <- FindAllMarkers(object = SO.traj1, min.pct = 0.25, thresh.use = 0.25)
markers <- subset(markers, p_val_adj < 0.05)

order.genes <- unique(as.character(markers$gene))

# ##### Constructing Single Cell Trajectories
cds <- setOrderingFilter(cds, order.genes)
cds <- reduceDimension(cds = cds, max_components = 3,method = 'DDRTree')
cds <- orderCells(cds)


## Fig 5b Pseudotime Trajectory

pdf("Fig5d.MP_M1_pseudotime_trajectory.pdf", 400/72,200/72)
plot_cell_trajectory(cds, color_by = "clusters", show_branch_points = F, cell_size = 0.5) + scale_x_reverse()
dev.off()


## Fig 5c Pseudotime Heatmap
order.genes <- order.genes[!grepl("ENSMPUG", order.genes)]
pdf("Fig5e.MP_M1_pseudotime_heatmap.pdf", 5/2.54*2,9/2.54*1.5)
plot_pseudotime_heatmap(cds[order.genes,], num_clusters = 4, cores = 1, show_rownames = T, return_heatmap = T)
dev.off()





## Fig5 d,e

### calculate pseudotime

SO.traj2 <- MP[,MP$Annotation %in% c("Monocyte-derived infiltrating macrophage", "SPP1hi CHIT1int profibrogenic M2")]

### 
# ##### make input data 
pd <- new('AnnotatedDataFrame', data = SO.traj2@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(SO.traj2), row.names = row.names(SO.traj2)))

cds <- newCellDataSet(as(SO.traj2@assays$SCT@data, "matrix"), phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

cds$clusters <- SO.traj2$Annotation

### 
# ##### Estimate size factors and dispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# ##### Filtering low-quality cells 
set.seed(123)
cds <- detectGenes(cds, min_expr = 0.1)

# ##### Select genes use Seurat::FindAllMarkers function
markers <- FindAllMarkers(object = SO.traj2, min.pct = 0.25, thresh.use = 0.25)
markers <- subset(markers, p_val_adj < 0.05)

order.genes <- unique(as.character(markers$gene))

# ##### Constructing Single Cell Trajectories
cds <- setOrderingFilter(cds, order.genes)
cds <- reduceDimension(cds = cds, max_components = 3,method = 'DDRTree')
cds <- orderCells(cds)

## Fig 5d Pseudotime Trajectory

pdf("Fig5d.MP_M2_pseudotime_trajectory.pdf", 400/72,200/72)
plot_cell_trajectory(cds, color_by = "clusters", show_branch_points = F, cell_size = 0.5) + scale_x_reverse()
dev.off()

## Fig 5e Pseudotime Heatmap
order.genes <- order.genes[!grepl("ENSMPUG", order.genes)]
pdf("Fig5e.MP_M2_pseudotime_heatmap.pdf", 5/2.54*2,9/2.54*1.5)
plot_pseudotime_heatmap(cds[order.genes,], num_clusters = 4, cores = 1, show_rownames = T, return_heatmap = T)
dev.off()


