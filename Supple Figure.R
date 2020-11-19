## install packages
install.packages('Seurat')

install.packages("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("piano")

install.packages("msigdbr")

install.packages("ggpubr")

install.packages("cowplot")

## load packages
library(Seurat)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(piano)
library(msigdbr)


## load seurat object
MPcov <- readRDS("/home/users/kjyi/Projects/covid19/ferret/singlecell/data/Seurat_object_total_cells.Rds")
MP <- readRDS("/home/users/kjyi/Projects/covid19/ferret/singlecell/data/Seurat_object_monocyte_macrophage_DC.Rds")

## FigS1a whole_umap

pdf("FigS1a.whole_umap_seurat_clusters.pdf",900/72*0.8, 688/72*0.8) 
DimPlot(MPcov,group.by = "seurat_clusters", label=T)
dev.off()


## FigS1b whole_umap

pdf("FigS1b.whole_stacked_bar_composition_celltype.pdf",3.8,5)
x = table(MPcov$Annotation,MPcov$Experimentwlrmal_group)
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


## FigS1c_Feature plot

pdf("FigS1c.featureplot.pdf",12,8*0.8)
goi <- c("EPCAM", "HBA2", "CD3E", "CD8A", "CD4",
         "TRDC", "CD19", "CD79A", "CD14", "FCGR3B")
p = FeaturePlot(MPcov, features = goi,reduction="umap",ncol=2,combine = FALSE,cols = c("grey80", "#0A64A4"))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p,ncol = 2)
dev.off()


## FigS1d_batch ; whole

pdf("FigS1d.condition_umap_whole.pdf",9,15)
par(mfrow=c(1,3), mar = c(1,1,3,1),pty="s", oma=c(0,0,0,0))
for(gr in c("Ct","C2","C5")){
  plot(Embeddings(MPcov,reduction="umap")[MPcov$Experimental_group!=gr,],col="#00000020",
       pch=16,cex=0.8,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
  points(Embeddings(MPcov,reduction="umap")[MPcov$Experimental_group==gr,],col="#168C4060",pch=16,cex=1)
  title(c("Ct"="dpi 0","C2"="dpi 2", "C5"="dpi 5")[gr])
}
dev.off()


## FigS1E_virus containing cells

MPcov$viral_count = MPcov@assays$RNA@counts["SARS-CoV-2",]
pdf("FigS1e.virus_containing_whole.pdf",9,15)
plot(Embeddings(MPcov, "umap"), pch=16, col="grey", bty="n", cex=1, xaxt="n", yaxt="n")
points(Embeddings(MPcov,"umap")[MPcov$viral_count>0,],col="red",pch=2)
dev.off()


## FigS2a_batch ; NK

pdf("FigS2a.condition_umap_NK.pdf",9,15)
par(mfrow=c(1,3), mar = c(1,1,3,1),pty="s", oma=c(0,0,0,0))
for(gr in c("Ct","C2","C5")){
  plot(Embeddings(cov.NK,reduction="umap")[cov.NK$Experimental_group!=gr,],col="#00000020",
       pch=16,cex=0.8,bty="n",xaxt="n",yaxt="n",cov.NKab="",ylab="")
  points(Embeddings(cov.NK,reduction="umap")[cov.NK$Experimental_group==gr,],col="#168C4060",pch=16,cex=1)
  title(c("Ct"="dpi 0","C2"="dpi 2", "C5"="dpi 5")[gr])
}
dev.off()


## FigS2b_IFN response

markers.NK <- FindAllMarkers(cov.NK, only.pos = T)
write.csv(markers.NK, "markers.NK.csv")


## FigS2c_bar plot ; CD8

pdf("FigS2c.CD8_stacked_bar_composition_celltype.pdf",3.8,5)
x = table(cov.CD8$seurat_clusters,cov.CD8$Experimental_group)
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


## FigS3a_whole_umap ; MP

pdf("FigS3a.whole_umap_seurat_clusters_MP.pdf",900/72*0.8, 688/72*0.8) 
DimPlot(MP, group.by = "seurat_clusters", label = T)
dev.off()


## FigS3b_bar plot ; MP

pdf("FigS3b.MP_stacked_bar_composition_celltype.pdf",3.8,5)
x = table(MP$Annotation,MP$Experimental_group)
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


## FigS3c_Feature plot

pdf("FigS3c.featureplot_MP.pdf",12,8*0.8)
goi <- c("CDKN2A", "APOE", "SPP1", "DDX60", "FCGR3B",
         "MARCO", "S100A8", "APOBEC3G", "CSF3R", "CD163", 
         "PPARG", "RGS2", "TOP2A", "STAB1", "WFDC2")
p = FeaturePlot(MPcov, features = goi,reduction="umap",ncol=2,combine = FALSE,cols = c("grey80", "#0A64A4"))
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p,ncol = 2)
dev.off()


## FigS3d_batch ; MP

pdf("FigS3d.condition_umap_MP.pdf",9,15)
par(mfrow=c(1,3), mar = c(1,1,3,1),pty="s", oma=c(0,0,0,0))
for(gr in c("Ct","C2","C5")){
  plot(Embeddings(MP,reduction="umap")[MP$Experimental_group!=gr,],col="#00000020",
       pch=16,cex=0.8,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
  points(Embeddings(MP,reduction="umap")[MP$Experimental_group==gr,],col="#168C4060",pch=16,cex=1)
  title(c("Ct"="dpi 0","C2"="dpi 2", "C5"="dpi 5")[gr])
}
dev.off()


## FigS3E_virus containing cells

MP$viral_count = MP@assays$RNA@counts["SARS-CoV-2",]

pdf("FigS3e.virus_containing_MP.pdf",9,15)
plot(Embeddings(MP, "umap"), pch=16, col="grey", bty="n", cex=1, xaxt="n", yaxt="n")
points(Embeddings(MP,"umap")[MP$viral_count>0,],col="red",pch=2)
dev.off()


## FigS3f GSEA_MP

### Define genes to show
markers_to_show3 <- markers_to_show %>%
  arrange(desc(avg_logFC)) %>% 
  {.[!duplicated(.$gene),]} %>%
  dplyr::filter(!gene %in% MP@misc$Dump$gene,
                !grepl("ENSMPUG",gene),
                !cluster == "Unclassified",
                pct.2<0.8) %>%
  mutate(cluster = factor(cluster, levels=levels(MP$Annotation))) %>%
  group_by(cluster) %>% dplyr::slice(1:50)


### load gene set
m_go.bp <- msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'BP') # Gene Ontology: biologic process
pat2 = "T_HELPER|T_CELL|EOSINOPHIL|POSITIVE|NEGATIVE"
m_go.bp <- m_go.bp[!grepl(pat2, m_go.bp$gs_name),]
m_go.bp <- piano::loadGSC(m_go.bp[,c("human_gene_symbol","gs_name")])


### depict GSEA bar plot

pdf("FigS3f.MP.goterm.enrichment_tight_set_MP.pdf",13,10)
par(mfrow=c(5,2),mar=c(3,35,2,1))
for(gr in unique(markers_to_show3$cluster)[c("Weakly activated M1", "Highly activated M1", "Proliferating macrophage", "Engulfing macrophage")]){
  goi <- markers_to_show3$gene[markers_to_show3$cluster == gr]
  universe = intersect(unique(unlist(m_go.bp$gsc)), rownames(MP))
  x <- runGSAhyper(genes = goi, gsc = m_go.bp, universe = universe, gsSizeLim = c(1,Inf), adjMethod = "BH")
  x$resTab[order(x$resTab[,"p-value"]),][10:1,1] %>% {-log10(.)} %>% 
  {names(.) = str_replace_all(names(.),"GO_","");
  names(.) = str_replace_all(names(.),"_"," ") %>% tolower() %>% Hmisc::capitalize(); .} %>%
    barplot(horiz=T,las=1)
  mtext(gr)
}
dev.off()