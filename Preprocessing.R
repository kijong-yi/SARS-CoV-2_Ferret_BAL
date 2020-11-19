

library(tidyverse)
library(Seurat)
library(SoupX)

# ------------------------------------------------------------------------------
#  1. SoupX
# ------------------------------------------------------------------------------

Ctrl_1 <- load10X(dataDir = "cellranger/Ctrl-1/outs")
Ctrl_2 <- load10X(dataDir = "cellranger/Ctrl-2/outs")
Ctrl_3 <- load10X(dataDir = "cellranger/Ctrl-3/outs")
C2_1 <- load10X(dataDir = "cellranger/C2-1/outs")
C2_2 <- load10X(dataDir = "cellranger/C2-2/outs")
C2_3 <- load10X(dataDir = "cellranger/C2-3/outs")
C5_1 <- load10X(dataDir = "cellranger/C5-1/outs")
C5_2 <- load10X(dataDir = "cellranger/C5-2/outs")
C5_3 <- load10X(dataDir = "cellranger/C5-3/outs")
C5_4 <- load10X(dataDir = "cellranger/C5-4/outs")


Ctrl_1.soupx = autoEstCont(Ctrl_1) # global rho = 0.01
Ctrl_2.soupx = autoEstCont(Ctrl_2) # global rho = 0.01
Ctrl_3.soupx = autoEstCont(Ctrl_3) # global rho = 0.03
C2_1.soupx = autoEstCont(C2_1) # global rho = 0.01
C2_2.soupx = autoEstCont(C2_2) # global rho = 0.02
C2_3.soupx = autoEstCont(C2_3) # global rho = 0.03
C5_1.soupx = autoEstCont(C5_1) # global rho = 0.01
C5_2.soupx = autoEstCont(C5_2) # global rho = 0.02
C5_3.soupx = autoEstCont(C5_3) # global rho = 0.02
C5_4.soupx = autoEstCont(C5_4) # global rho = 0.01

Ctrl_1.soupx %>% write_rds("data/soupX.autoEstCont.Ctrl_1.Rds")
Ctrl_2.soupx %>% write_rds("data/soupX.autoEstCont.Ctrl_2.Rds")
Ctrl_3.soupx %>% write_rds("data/soupX.autoEstCont.Ctrl_3.Rds")
C2_1.soupx %>% write_rds("data/soupX.autoEstCont.C2_1.Rds")
C2_2.soupx %>% write_rds("data/soupX.autoEstCont.C2_2.Rds")
C2_3.soupx %>% write_rds("data/soupX.autoEstCont.C2_3.Rds")
C5_1.soupx %>% write_rds("data/soupX.autoEstCont.C5_1.Rds")
C5_2.soupx %>% write_rds("data/soupX.autoEstCont.C5_2.Rds")
C5_3.soupx %>% write_rds("data/soupX.autoEstCont.C5_3.Rds")
C5_4.soupx %>% write_rds("data/soupX.autoEstCont.C5_4.Rds")

Ctrl_1 <- adjustCounts(Ctrl_1)
Ctrl_2 <- adjustCounts(Ctrl_2)
Ctrl_3 <- adjustCounts(Ctrl_3)
C2_1 <- adjustCounts(C2_1)
C2_2 <- adjustCounts(C2_2)
C2_3 <- adjustCounts(C2_3)
C5_1 <- adjustCounts(C5_1)
C5_2 <- adjustCounts(C5_2)
C5_3 <- adjustCounts(C5_3)
C5_4 <- adjustCounts(C5_4)



Ctrl_1 %>% write_rds("data/SoupX.adjCount.out.Ctrl_1.out.Rds")
Ctrl_2 %>% write_rds("data/SoupX.adjCount.out.Ctrl_2.out.Rds")
Ctrl_3 %>% write_rds("data/SoupX.adjCount.out.Ctrl_3.out.Rds")
C2_1 %>% write_rds("data/SoupX.adjCount.out.C2_1.out.Rds")
C2_2 %>% write_rds("data/SoupX.adjCount.out.C2_2.out.Rds")
C2_3 %>% write_rds("data/SoupX.adjCount.out.C2_3.out.Rds")
C5_1 %>% write_rds("data/SoupX.adjCount.out.C5_1.out.Rds")
C5_2 %>% write_rds("data/SoupX.adjCount.out.C5_2.out.Rds")
C5_3 %>% write_rds("data/SoupX.adjCount.out.C5_3.out.Rds")
C5_4 %>% write_rds("data/SoupX.adjCount.out.C5_4.out.Rds")


# ------------------------------------------------------------------------------
#  2. Gene name annotation
# ------------------------------------------------------------------------------


library(readr)
library(dplyr)
library(Seurat)

load_biomart <- function(xml="",col_names=F) {
  xml <- stringr::str_replace_all(xml,"\n","")
  xml <- stringr::str_replace_all(xml,"\t","")
  cmd <- paste0("wget -O - 'http://www.ensembl.org/biomart/martservice?query=",xml,"'")
  readr::read_tsv(system(cmd,intern=T),col_names=col_names)
}

biomart <- '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

<Dataset name = "mpfuro_gene_ensembl" interface = "default" >
<Attribute name = "ensembl_gene_id" />
<Attribute name = "hsapiens_homolog_ensembl_gene" />
<Attribute name = "hsapiens_homolog_associated_gene_name" />
</Dataset>
</Query>' %>% load_biomart()



gene_list = read_tsv("cellranger/C5-3/outs/filtered_feature_bc_matrix/features.tsv.gz",
										 col_names = c("rawid","rawname","type"))
gene_list$id = str_replace(gene_list$rawid,"MusPutFur1.0_","")
gene_list$name = str_replace(gene_list$rawname,"MusPutFur1.0_","")
gene_list

biomart
biomart1 = structure(biomart$X3,names=biomart$X1)

gene_list$ortholog = biomart1[str_replace(gene_list$id,"MusPutFur1.0_","")]

gene_list$ortholog[1:4] = c("SARS-CoV-2-5UTR","SARS-CoV-2","SARS-CoV-5UTR","SARS-CoV")
gene_list$id[1:4] = c("SARS-CoV-2-5UTR","SARS-CoV-2","SARS-CoV-5UTR","SARS-CoV")
gene_list$name[1:4] = c("SARS-CoV-2-5UTR","SARS-CoV-2","SARS-CoV-5UTR","SARS-CoV")

table(!grepl("ENSMPUG",gene_list$name)) # annotated by ensembl 16805/32059
table(!is.na(gene_list$ortholog)) # annotated by ensembl 16850/32059

table(!grepl("ENSMPUG",gene_list$name),!is.na(gene_list$ortholog))
table(duplicated(na.omit(gene_list$ortholog)))

gene_list$final = gene_list$id

gene_list$final[!grepl("ENSMPUG",gene_list$name)] = gene_list$name[!grepl("ENSMPUG",gene_list$name)]

gene_list$final[!grepl("ENSMPUG",gene_list$final)] = gene_list$name[!grepl("ENSMPUG",gene_list$name)]

for(i in 1:nrow(gene_list)){
	if(!grepl("ENSMPUG",gene_list$final[i])){next}
	if(!is.na(gene_list$ortholog[i])){
		if(!gene_list$ortholog[i] %in% gene_list$name & 
			 !gene_list$ortholog[i] %in% gene_list$ortholog[-i]){
			gene_list$final[i] = gene_list$ortholog[i]
		}
	}
}

table(!grepl("ENSMPUG",gene_list$name))
table(!grepl("ENSMPUG",gene_list$final))
grepl("ENSMPUG",gene_list$final)
# 1322 genes were additionally annotated

gene_list$final2 <- make.unique(gene_list$final)
gene_list$name_ = str_replace_all(gene_list$rawname, "MusPutFur1.0_", "MusPutFur1.0-")
gene_list
write_rds(gene_list,"data/gene_annotation.Rds")



# ------------------------------------------------------------------------------
# 3. Seurat processing and annotation 
# ------------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
# library(future)
# plan("multiprocess", workers = 6)
# options(future.globals.maxSize = 200 * 1000 * 1024^2)


Ctrl_1 <- read_rds('data/SoupX.adjCount.out.Ctrl_1.out.Rds')
Ctrl_2 <- read_rds('data/SoupX.adjCount.out.Ctrl_2.out.Rds')
Ctrl_3 <- read_rds('data/SoupX.adjCount.out.Ctrl_3.out.Rds')
C2_1 <- read_rds('data/SoupX.adjCount.out.C2_1.out.Rds')
C2_2 <- read_rds('data/SoupX.adjCount.out.C2_2.out.Rds')
C2_3 <- read_rds('data/SoupX.adjCount.out.C2_3.out.Rds')
C5_1 <- read_rds('data/SoupX.adjCount.out.C5_1.out.Rds')
C5_2 <- read_rds('data/SoupX.adjCount.out.C5_2.out.Rds')
C5_3 <- read_rds('data/SoupX.adjCount.out.C5_3.out.Rds')
C5_4 <- read_rds('data/SoupX.adjCount.out.C5_4.out.Rds')


gene_list <- read_rds("data/gene_annotation.Rds")

# change virus gene name
gene_list$final2[1:(nrow(gene_list))] %>% head
gene_list$final2[1:(nrow(gene_list))] %>% tail
gene_list$final2[1:4] <- c("SARS-CoV-2_5","SARS-CoV-2","SARS-CoV_5","SARS-CoV")


rownames(Ctrl_1) <- gene_list$final2
rownames(Ctrl_2) <- gene_list$final2
rownames(Ctrl_3) <- gene_list$final2
rownames(C2_1) <- gene_list$final2
rownames(C2_2) <- gene_list$final2
rownames(C2_3) <- gene_list$final2
rownames(C5_1) <- gene_list$final2
rownames(C5_2) <- gene_list$final2
rownames(C5_3) <- gene_list$final2
rownames(C5_4) <- gene_list$final2

colnames(Ctrl_1)[1:3]
colnames(Ctrl_1) = paste0("Ctrl-1:",colnames(Ctrl_1),"x")
colnames(Ctrl_2) = paste0("Ctrl-2:",colnames(Ctrl_2),"x")
colnames(Ctrl_3) = paste0("Ctrl-3:",colnames(Ctrl_3),"x")
colnames(C2_1) = paste0("C2-1:",colnames(C2_1),"x")
colnames(C2_2) = paste0("C2-2:",colnames(C2_2),"x")
colnames(C2_3) = paste0("C2-3:",colnames(C2_3),"x")
colnames(C5_1) = paste0("C5-1:",colnames(C5_1),"x")
colnames(C5_2) = paste0("C5-2:",colnames(C5_2),"x")
colnames(C5_3) = paste0("C5-3:",colnames(C5_3),"x")
colnames(C5_4) = paste0("C5-4:",colnames(C5_4),"x")


merged_mat = cbind(C0_1,C0_2,C0_3,C2_1,C2_2,C2_3,C5_1,C5_2,C5_3,C5_4)

MPcov <- CreateSeuratObject(counts = merged_mat)

MPcov$ID_1  = colnames(MPcov) %>% substr(1,4)
MPcov$expgr = colnames(MPcov) %>% substr(1,2)
MPcov$disease = c("C0"="0","C2"="C","C5"="C","S2"="S","S5"="S")[substr(colnames(MPcov),1,2)]

MPcov$ID_1  = colnames(MPcov) %>% str_replace(":.*","")

findallmarkers_mc <- function(seurat_obj, ...){
  n_clust <- unique(Idents(seurat_obj))
  mcFindMarkers <- function(i){
    ident1 <- i
    ident2 <- n_clust[n_clust != i]
    table <- FindMarkers(seurat_obj,
                         ident.1 = ident1, ident.2 = ident2, ...)
    table$Gene.name.uniq <- rownames(table)
    table$cluster <- rep(i, nrow(table))
    return(table)
  }
  marker_results <- list()[n_clust]
  marker_results <- parallel::mclapply(n_clust, mcFindMarkers, mc.cores = 1)
  marker_results
}


MPcov <- SCTransform(MPcov)
MPcov <- RunPCA(MPcov)
p1 <- ElbowPlot(MPcov)
MPcov <- JackStraw(object = MPcov, num.replicate = 100)
MPcov <- ScoreJackStraw(object = MPcov, dims = 1:20)
p2 <- JackStrawPlot(object = MPcov, dims = 1:20)
CombinePlots(plot = list(p1, p2))
MPcov <- FindNeighbors(MPcov, dims = 1:20)
MPcov <- FindClusters(MPcov,resolution = 0.5)
MPcov <- FindClusters(MPcov,resolution = 0.3)
MPcov <- FindClusters(MPcov,resolution = 0.2)
MPcov <- FindClusters(MPcov,resolution = 0.1)
MPcov <- FindClusters(MPcov,resolution = 0.4)
MPcov <- RunUMAP(MPcov, dims = 1:20)
MPcov <- RunTSNE(MPcov, dims = 1:20)

Idents(MPcov) <- "SCT_snn_res.0.3"
MPcov@misc$MPcov_res0.3_logFC0.2_minpct0.2 <- findallmarkers_mc(MPcov, min.pct = 0.2, logfc.threshold = 0.2, only.pos = TRUE)



SO.ferret %>% write_rds("data/SO.ferret.Rds")



# subset: monocyte-macrophage-DC, Lymphocytes

MP <- MPcov[,MPcov$annot %in% c("tissue macrophage",
                                "Siglec F Macrophage",
                                "infiltrating macrophage",
                                "Mast cell",
                                "Dendritic cell","pacman")]
L <- MPcov[,MPcov$annot %in% c("CD4 T cell",
                               "CD8 T cell",
                               "Dividing T cells",
                               "rd T cell",
                               "NK cell")]


MP <- RunPCA(MP)
ElbowPlot(MP)
MP <- FindNeighbors(MP, dims = 1:12)
MP <- FindClusters(MP,resolution = 0.5)
MP <- FindClusters(MP,resolution = 0.4)
MP <- FindClusters(MP,resolution = 0.2)
MP <- FindClusters(MP,resolution = 0.1)
MP <- FindClusters(MP,resolution = 0.3)
MP <- RunUMAP(MP, dims = 1:12)
MP <- RunUMAP(MP,dims = 1:7,reduction.name = "umap_7")
MP <- RunUMAP(MP,dims = 1:10,reduction.name = "umap_10")
MP <- RunUMAP(MP,dims = 1:12,reduction.name = "umap_12")
MP <- RunUMAP(MP,dims = 1:15,reduction.name = "umap_15")
MP <- RunUMAP(MP,dims = 1:20,reduction.name = "umap_20")



L <- RunPCA(L)
ElbowPlot(L)
L <- FindNeighbors(L, dims = 1:7)
L <- FindClusters(L,resolution = 0.5)
L <- FindClusters(L,resolution = 0.4)
L <- FindClusters(L,resolution = 0.2)
L <- FindClusters(L,resolution = 0.1)
L <- FindClusters(L,resolution = 0.3)
L <- RunUMAP(L, dims = 1:7)

write_rds(MP,"data/SO.ferret.MP.rds")
write_rds(L,"data/SO.ferret.L.rds")



MPcov$ID_1 %>% table
Idents(MPcov) = "ID_1"
Idents(MP) = "ID_1"
case_specific_markers <- findallmarkers_mc(MPcov, test.use="roc", min.pct = "0.1", only.pos=T)
case_specific_markers.MP <- findallmarkers_mc(MP, test.use="roc", min.pct = "0.1", only.pos=T)
# write_csv(rownames_to_column(do.call(rbind,case_specific_markers),"gene"),"case_specific_markers_whole.csv")

genes_to_exclude = c("HLA-DQA1", "ENSMPUG00000007244")

MPcov = MPcov[rownames(MPcov)[!rownames(MPcov) %in% c("HLA-DQA1", "ENSMPUG00000007244")], ]
MP = MP[rownames(MP)[!rownames(MP) %in% c("HLA-DQA1", "ENSMPUG00000007244")], ]
L = L[rownames(L)[!rownames(L) %in% c("HLA-DQA1", "ENSMPUG00000007244")], ]

MPcov <- RunPCA(MPcov)
MPcov <- FindNeighbors(MPcov, dims = 1:20)
MPcov <- FindClusters(MPcov, resolution = 0.5)
MPcov <- RunUMAP(MPcov, dims = 1:20)

MPcov$finalannot_nonum = "NA"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "0"]  = "M2 macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "1"]  = "M2 macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "2"]  = "CD4 T cell"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "3"]  = "M2 macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "4"]  = "M1 macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "5"]  = "M2 macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "6"]  = "Proliferating macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "7"]  = "M1 macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "8"]  = "M2 macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "9"]  = "CD8 T cell"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "10"] = "Gamma delta T cell"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "11"] = "Plasma cell"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "12"] = "B cell"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "13"] = "M2 macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "14"] = "Engulfing macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "15"] = "Engulfing macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "16"] = "Mast cell"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "17"] = "Doublet"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "18"] = "Doublet"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "19"] = "Epithelial cell"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "20"] = "Dendritic cell"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "21"] = "NK cell"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "22"] = "Proliferating T cell"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "23"] = "Granulocyte"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "24"] = "M2 macrophage"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "25"] = "Doublet"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "26"] = "RBC"
MPcov$finalannot_nonum[MPcov$SCT_snn_res.0.4 == "27"] = "Epithelial cell"
# MPcov %>% write_rds("data/SO.ferret.Rds")



# write_rds(MPcov, "data/SO.ferret.Rds")

MP <- RunPCA(MP)
MP <- FindNeighbors(MP, dims = 1:13)
MP <- FindClusters(MP, resolution = 0.5)
# MP <- RunUMAP(MP, dims = 1:12, min.dist = 0.2, n.neighbors = 30)
MP <- RunUMAP(MP, dims = 1:10)

# write_rds(SO.ferret.MP, "data/SO.ferret.MP2.Rds")

L <- RunPCA(L)
L <- FindNeighbors(L, dims = 1:7)
L <- FindClusters(L,resolution = 0.3)
L <- RunUMAP(L, dims = 1:7)

MP$Annotation = "NA"
MP$Annotation[MP$SCT_snn_res.0.6 == "0"]  = "SPP1 low APOE neg FABP4 high tissue M2"
MP$Annotation[MP$SCT_snn_res.0.6 == "1"]  = "Interferon stimulated M2"
MP$Annotation[MP$SCT_snn_res.0.6 == "2"]  = "SPP1 low APOE neg FABP4 high tissue M2"
MP$Annotation[MP$SCT_snn_res.0.6 == "3"]  = "SPP1 high fibrogenic M2"
MP$Annotation[MP$SCT_snn_res.0.6 == "4"]  = "APOE pos FABP4 high tissue M2"
MP$Annotation[MP$SCT_snn_res.0.6 == "5"]  = "Interferon stimulated M1"
MP$Annotation[MP$SCT_snn_res.0.6 == "6"]  = "APOE pos FABP4 high tissue M2"
MP$Annotation[MP$SCT_snn_res.0.6 == "7"]  = "Interferon stimulated M2"
MP$Annotation[MP$SCT_snn_res.0.6 == "8"]  = "Transitional M1"
MP$Annotation[MP$SCT_snn_res.0.6 == "9"]  = "Infiltrating macrophage"
MP$Annotation[MP$SCT_snn_res.0.6 == "10"] = "Engulfing macrophage"
MP$Annotation[MP$SCT_snn_res.0.6 == "11"] = "Proliferating macrophage"
MP$Annotation[MP$SCT_snn_res.0.6 == "12"] = "Proliferating macrophage"
MP$Annotation[MP$SCT_snn_res.0.6 == "13"] = "SPP1 high fibrogenic M2"
MP$Annotation[MP$SCT_snn_res.0.6 == "14"] = "Engulfing macrophage"
MP$Annotation[MP$SCT_snn_res.0.6 == "15"] = "Interferon stimulated M2"
MP$Annotation[MP$SCT_snn_res.0.6 == "16"] = "Unclassified"
MP$Annotation = factor(MP$Annotation,levels=c("SPP1 low APOE neg FABP4 high tissue M2",
                                              "APOE pos FABP4 high tissue M2",
                                              "SPP1 high fibrogenic M2",
                                              "Interferon stimulated M2",
                                              "Transitional M1",
                                              "Interferon stimulated M1",
                                              "Infiltrating macrophage",
                                              "Engulfing macrophage",
                                              "Proliferating macrophage",
                                              "Unclassified"))

# write_rds(MP, "data/SO.ferret.MP.rds")

L <- RunPCA(L)
ElbowPlot(L)
L <- FindNeighbors(L, dims = 1:7)
L <- FindClusters(L,resolution = 0.5)
L <- FindClusters(L,resolution = 0.4)
L <- FindClusters(L,resolution = 0.2)
L <- FindClusters(L,resolution = 0.1)
L <- FindClusters(L,resolution = 0.3)
L <- RunUMAP(L, dims = 1:7)


# End of script

