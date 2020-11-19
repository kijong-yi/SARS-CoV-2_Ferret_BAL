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


## Fig3a MP umap

pdf("Fig3a.MP_umap_seurat_clusters.pdf",900/72*0.8, 688/72*0.8) 
DimPlot(MP,group.by = "Annotation", label=T)
dev.off()


## Fig3b marker gene, dot plot

goi <- c("FABP4", "PPARG", "APOE", "APOC1", "DDX60", 
         "SPP1", "CSF3R", "IL1B", "RGS2", "ISG15", 
         "CHIT1", "STAB1", "TOP2A", "WFDC2")

Idents(MP) <- "Annotation"

p <- DotPlot(MP,
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

p$data$id = factor(p$data$id, levels(MP$Annotation) %>% rev)

pdf("Fig3b.dotplot.pdf",7.5*1.2,5*1.2)
p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


## Fig3c MP proportion

colorder = c("Ctrl-1","Ctrl-2","Ctrl-3",
             "C2-1","C2-2","C2-3",
             "C5-1","C5-2", "C5-3","C5-4")

x <- table(MP$Annotation,MP$Sample_name)
x <- x[, colorder]
x3= t(t(x)/rowSums(t(x)))

x4 = as.data.frame(as.table(t(x3)))
colnames(x4) = c("sample","celltype","Freq")
x4$group = x4$sample %>% str_replace("-.*","")
x4$group = factor(x4$group, levels = c("Ctrl","C2","C5"))

write.csv(x4, "Fig3c.proportion_each_cluster.csv")


## Fig3d Heatmap; MP

### define dump genes

MPcov <- readRDS("/home/users/kjyi/Projects/covid19/ferret/singlecell/data/Seurat_object_total_cells.Rds")

Idents(MPcov) <- "Annotation"
Dump <- c("CD8 T cell","CD4 T cell","Epithelial cell","RBC","Plasma cell","B cell") %>%
	lapply(function(x){FindMarkers(MPcov, ident.1 = x,only.pos = T, min.pct = 0.3)}) %>%
	do.call(rbind,.)

rm(MPcov)

Dump$gene = rownames(Dump)

MP@misc$Dump = Dump


### define genes to show

markers_to_show <- FindAllMarkers(MP, only.pos = T)

markers_to_show_2 <- markers_to_show %>% arrange(desc(avg_logFC)) %>% 
  {.[!duplicated(.$gene),]} %>%
  dplyr::filter(!gene %in% MP@misc$Dump$gene,
                !grepl("ENSMPUG",gene),
                !cluster == "Unclassified",
                pct.2<0.8) %>%
                group_by(cluster) %>% dplyr::slice(1:20)


### define DoMultiBarHeatmap function

suppressPackageStartupMessages({
  library(rlang)
})

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

pdf("Fig3d.heatmap_MP.pdf",13,9)
DoMultiBarHeatmap(MP,
                  features = markers_to_show_2$gene, 
                  group.by = 'Annotation',
                  disp.min = -2.5,disp.max = 2.5,
                  additional.group.by = 'Experimental_group',
                  size = 3,
                  label = F) + 
  scale_fill_gradient2(low = "magenta", 
                       mid = "black", 
                       high = "yellow", 
                       midpoint = 0, guide = "colourbar", aesthetics = "fill")+
  theme(axis.text.y = element_text(size = 7))
dev.off()


## Fig3e GSEA_MP

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

pdf("Fig3e.MP.goterm.enrichment_tight_set_MP.pdf",13,10)
par(mfrow=c(5,2),mar=c(3,35,2,1))
for(gr in unique(markers_to_show3$cluster)[c("Resting tissue macrophage", "APOE+ tissue macrophage", "SPP1hi CHIT1int profibrogenic M2", "Activated tissue macrophage", "Monocyte-derived infiltrating macrophage")]){
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
