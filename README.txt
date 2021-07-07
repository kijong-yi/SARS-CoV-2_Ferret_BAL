[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5078140.svg)](https://doi.org/10.5281/zenodo.5078140)


Processed data files (Rds) can be downloaded from links below
 https://drive.google.com/file/d/1NG9l_utCUBgJRdOCdEBHm1FIkh7LOjGV/view?usp=sharing
 https://drive.google.com/file/d/1D-T7fTjwc34hmJ68SBqaOMrTQwZm3akp/view?usp=sharing

Raw sequencing data will be available later in GEO database.



My environment: (not all are neccessary)

```
R version 3.6.0 (2019-04-26)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] splines   stats4    parallel  grid      stats     graphics  grDevices utils     datasets 
[10] methods   base     

other attached packages:
 [1] piano_2.2.0                 doMC_1.3.5                  iterators_1.0.12           
 [4] foreach_1.5.0               cowplot_1.0.0               GSVA_1.34.0                
 [7] kjyi_0.1.0                  GOplot_1.0.2                gridExtra_2.3              
[10] ggdendro_0.1.22             pathview_1.26.0             org.Hs.eg.db_3.10.0        
[13] AnnotationDbi_1.48.0        msigdbr_7.2.1               openxlsx_4.1.4             
[16] readxl_1.3.1                viridis_0.5.1               viridisLite_0.3.0          
[19] RColorBrewer_1.1-2          reshape_0.8.8               ggrepel_0.8.2              
[22] DESeq2_1.26.0               SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
[25] BiocParallel_1.20.1         matrixStats_0.56.0          GenomicRanges_1.38.0       
[28] GenomeInfoDb_1.22.1         IRanges_2.20.2              S4Vectors_0.24.4           
[31] pheatmap_1.0.12             sctransform_0.2.1           SeuratWrappers_0.1.0       
[34] ggpubr_0.2.5                magrittr_1.5                monocle_2.14.0             
[37] DDRTree_0.1.5               irlba_2.3.3                 VGAM_1.1-2                 
[40] Biobase_2.46.0              BiocGenerics_0.32.0         Matrix_1.2-17              
[43] circlize_0.4.8              ComplexHeatmap_2.2.0        immunarch_0.6.6            
[46] patchwork_1.0.0             data.table_1.12.8           dtplyr_1.0.1               
[49] forcats_0.5.0               stringr_1.4.0               dplyr_0.8.5                
[52] purrr_0.3.4                 readr_1.3.1                 tidyr_1.0.2                
[55] tibble_3.0.0                ggplot2_3.3.0               tidyverse_1.3.0            
[58] Seurat_3.1.4               

loaded via a namespace (and not attached):
  [1] rsvd_1.0.3             Hmisc_4.4-1            ica_1.0-2              class_7.3-15          
  [5] lmtest_0.9-37          crayon_1.3.4           MASS_7.3-51.4          nlme_3.1-139          
  [9] backports_1.1.6        qlcMatrix_0.9.7        reprex_0.3.0           rlang_0.4.5           
 [13] XVector_0.26.0         ROCR_1.0-7             limma_3.42.2           sets_1.0-18           
 [17] rjson_0.2.20           bit64_4.0.5            glue_1.4.0             UpSetR_1.4.0          
 [21] shinydashboard_0.7.1   haven_2.2.0            tidyselect_1.0.0       fitdistrplus_1.0-14   
 [25] XML_3.99-0.3           zoo_1.8-7              xtable_1.8-4           evaluate_0.14         
 [29] bibtex_0.4.2.2         Rdpack_0.11-1          cli_2.0.2              zlibbioc_1.32.0       
 [33] sn_1.6-1               rstudioapi_0.11        rpart_4.1-15           fastmatch_1.1-0       
 [37] shiny_1.4.0.2          xfun_0.13              clue_0.3-57            multtest_2.42.0       
 [41] cluster_2.0.8          caTools_1.18.0         KEGGREST_1.26.1        ape_5.3               
 [45] listenv_0.8.0          Biostrings_2.54.0      png_0.1-7              future_1.17.0         
 [49] withr_2.1.2            bitops_1.0-6           slam_0.1-47            plyr_1.8.6            
 [53] cellranger_1.1.0       GSEABase_1.48.0        sparsesvd_0.2          pillar_1.4.3          
 [57] gplots_3.0.3           GlobalOptions_0.1.1    multcomp_1.4-13        fs_1.4.1              
 [61] flexmix_2.3-15         kernlab_0.9-29         GetoptLong_0.1.8       vctrs_0.2.4           
 [65] ellipsis_0.3.0         generics_0.0.2         metap_1.3              tools_3.6.0           
 [69] foreign_0.8-71         munsell_0.5.0          fgsea_1.12.0           fastmap_1.0.1         
 [73] HSMMSingleCell_1.6.0   compiler_3.6.0         httpuv_1.5.2           plotly_4.9.2.1        
 [77] GenomeInfoDbData_1.2.2 lattice_0.20-38        visNetwork_2.0.9       relations_0.6-9       
 [81] mutoss_0.1-12          utf8_1.1.4             later_1.0.0            jsonlite_1.6.1        
 [85] scales_1.1.0           docopt_0.6.1           graph_1.64.0           pbapply_1.4-2         
 [89] genefilter_1.68.0      lazyeval_0.2.2         promises_1.1.0         latticeExtra_0.6-29   
 [93] reticulate_1.15        checkmate_2.0.0        rmarkdown_2.1          sandwich_2.5-1        
 [97] Rtsne_0.15             uwot_0.1.8             igraph_1.2.5           survival_3.1-12       
[101] numDeriv_2016.8-1.1    yaml_2.2.1             plotrix_3.7-7          prabclus_2.3-2        
[105] htmltools_0.4.0        memoise_1.1.0          modeltools_0.2-23      locfit_1.5-9.4        
[109] digest_0.6.25          assertthat_0.2.1       mime_0.9               rappdirs_0.3.1        
[113] densityClust_0.3       npsurv_0.4-0           RSQLite_2.2.1          future.apply_1.4.0    
[117] lsei_1.2-0             remotes_2.2.0          blob_1.2.1             fastICA_1.2-2         
[121] shinythemes_1.1.2      Formula_1.2-3          labeling_0.3           fpc_2.2-5             
[125] RCurl_1.98-1.1         broom_0.7.0            hms_0.5.3              modelr_0.1.6          
[129] colorspace_1.4-1       base64enc_0.1-3        BiocManager_1.30.10    mnormt_1.5-6          
[133] shape_1.4.4            nnet_7.3-12            Rcpp_1.0.4.6           mclust_5.4.6          
[137] RANN_2.6.1             mvtnorm_1.1-0          ggseqlogo_0.1          fansi_0.4.1           
[141] R6_2.4.1               factoextra_1.0.7       ggridges_0.5.2         lifecycle_0.2.0       
[145] zip_2.0.4              TFisher_0.2.0          ggsignif_0.6.0         gdata_2.18.0          
[149] leiden_0.3.3           robustbase_0.93-6      RcppAnnoy_0.0.16       TH.data_1.0-10        
[153] htmlwidgets_1.5.1      marray_1.64.0          rvest_0.3.5            globals_0.12.5        
[157] htmlTable_2.1.0        KEGGgraph_1.46.0       codetools_0.2-16       lubridate_1.7.8       
[161] FNN_1.1.3              gtools_3.8.2           dbplyr_1.4.2           RSpectra_0.16-0       
[165] gtable_0.3.0           tsne_0.1-3             DBI_1.1.0              ggalluvial_0.11.1     
[169] httr_1.4.1             KernSmooth_2.23-15     stringi_1.4.6          reshape2_1.4.4        
[173] farver_2.0.3           diptest_0.75-7         annotate_1.64.0        Rgraphviz_2.30.0      
[177] DT_0.13                xml2_1.3.1             combinat_0.0-8         shinyjs_1.1           
[181] geneplotter_1.64.0     DEoptimR_1.0-8         bit_4.0.4              jpeg_0.1-8.1          
[185] pkgconfig_2.0.3        gbRd_0.4-11            knitr_1.28        
```
