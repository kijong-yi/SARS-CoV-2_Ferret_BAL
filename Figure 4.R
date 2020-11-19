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


## load seurat object
MP <- readRDS("/home/users/kjyi/Projects/covid19/ferret/singlecell/data/Seurat_object_monocyte_macrophage_DC.Rds")


## define cell cluster name 
MP$Annotation <- as.character(MP$Annotation)
MP$Annotation <- as.character(MP$Annotation_2)

MP$Annotation[MP$Annotation_2 == "APOE pos FABP4 high tissue M2"] <- "APOE+ tissue macrophage"
MP$Annotation[MP$Annotation_2 == "SPP1 high fibrogenic M2"] <- "SPP1hi CHIT1int profibrogenic M2"
MP$Annotation[MP$Annotation_2 == "Transitional M1"] <- "Weakly activated M1"
MP$Annotation[MP$Annotation_2 == "Interferon stimulated M1"] <- "Highly activated M1"
MP$Annotation[MP$Annotation_2 == "Infiltrating macrophage"] <- "Monocyte-derived infiltrating macrophage"
MP$Annotation[MP$Annotation_2 == "APOE pos FABP4 high tissue M2"] <- "APOE+ tissue macrophage"


## Fig4a,b volcano plot

### Find DEG betwween C2 vs C5

namelist <- c("Resting tissue macrophage",
              "APOE+ tissue macrophage",
              "SPP1hi CHIT1int profibrogenic M2",
              "Activated tissue macrophage",
              "Weakly activated M1",
              "Highly activated M1",
              "Monocyte-derived infiltrating macrophage",
              "Engulfing macrophage",
              "Proliferating macrophage",
              "Unclassified")

Idents(MP) <- "Experimental_group"

x.list <- foreach(i = 1:10) %dopar% {
  FindMarkers(MP[, MP$Annotation == c(namelist[[i]])],
              only.pos = F,
              ident.1 = "C2", 
              ident.2 = "C5",
              min.pct = 0, logfc.threshold = 0)
}


### Depict volano plot

volcano.plot.d2d5 <- list()


for(j in 1:10){
  gene_list <- x.list[[j]]
  gene_list <- rownames_to_column(gene_list, var = "gene")
  rownames(gene_list) <- gene_list$gene
  gene_list <- gene_list %>% mutate(., rank = -log(p_val)*avg_logFC)
  
  gene_list$rank2 <- gene_list$rank
  table(grepl("ENSMPUG",gene_list$gene))
  table(gene_list$avg_logFC < 0.4 & gene_list$avg_logFC > -0.4)
  gene_list$rank2[grepl("ENSMPUG",gene_list$gene)] <- NA
  gene_list$rank2[gene_list$avg_logFC < 0.4 & gene_list$avg_logFC > -0.4] <- NA
  
  upgene <- arrange(gene_list, -rank2) %>% filter(avg_logFC > 0) %>% head(8) %>% na.omit() %>% .$gene
  dwgene <- arrange(gene_list, rank2) %>% filter(avg_logFC < 0) %>% head(8) %>%  na.omit() %>% .$gene
  genestoshow <- c(upgene, dwgene)
  print(genestoshow)

  gene_list_ordered   <- gene_list
  gene_list_ordered$gene <- ifelse(gene_list_ordered$gene %in% genestoshow, 
                                   gene_list_ordered$gene, "")
  
  gene_list_ordered$threshold = as.factor(ifelse(gene_list_ordered$p_val > 0.05, 'Not', 
                                                 ifelse(gene_list_ordered$p_val < 0.05 & gene_list_ordered$avg_logFC <= -0.4, 'Down',
                                                        ifelse(gene_list_ordered$p_val < 0.05 & gene_list_ordered$avg_logFC >= 0.4, 'Up', 'Not'))))
  
  
  volcano.plot.d2d5[[j]] <-  ggplot(gene_list_ordered, aes(x=-avg_logFC, y=-log10(p_val))) +
    ggrastr::geom_point_rast(aes(color = threshold), shape = 20, alpha = 1, size = 1.75) +
    scale_x_continuous(limits = c(-2, 2)) +
    scale_color_manual(values = c('blue', 'grey', 'red')) +
    ggtitle(paste(namelist[j], "d2 vs d5")) +
    xlab("log2 fold change") + 
    ylab("-log10 p-value") +
    geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.4, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.4, colour="#990000", linetype="dashed") +
    theme_classic() + 
    theme(legend.position = "none") +
    geom_text_repel(data=filter(gene_list_ordered, p_val<0.05), aes(label = gene, color = threshold), size = 4)
}

print(volcano.plot.d2d5[[1]])


pdf("Fig4a.volcano_plot_MP.pdf",12,8)
cowplot::plot_grid(plotlist = volcano.plot.d2d5[c(1,4,7,5,6,3)], ncol = 3)
dev.off()