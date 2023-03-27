## Export RNA counts and data slots with gene expression and cluster identity for Cre3 and Lox2 samples 
# Run by: António Sousa - UBI-IGC (agsousa@igc.gulbenkian.pt)
# Date: 11/06/2021

# Set seed 
set.seed(seed = 1024) # to keep reproducibility

## Import packages 
library("dplyr")
library("Seurat")
library("Matrix")

# Import Seurat object with the Cre3 & Lox2 integrated samples
seu <- readRDS(file = "../results/int_28_05_21/R_objects/seu.rds")

# Split obj by sample
sobj <- SplitObject(object = seu, split.by = "orig.ident")

# Check if RNA is the active default assay
stopifnot( all( c(DefaultAssay(sobj$Cre3), DefaultAssay(sobj$Lox2) ) == "RNA" ) )

# Get and export matrices
gene_exp[["Cre3"]] <- gene_exp[["Lox2"]] <- gene_exp <- list()
gene_exp[["Cre3"]][["counts"]] <- as.matrix( GetAssayData(sobj$Cre3[["RNA"]], slot = "counts") ) 
gene_exp[["Cre3"]][["data"]] <- as.matrix( GetAssayData(sobj$Cre3[["RNA"]], slot = "data") ) 
gene_exp[["Lox2"]][["counts"]] <- as.matrix( GetAssayData(sobj$Lox2[["RNA"]], slot = "counts") ) 
gene_exp[["Lox2"]][["data"]] <- as.matrix( GetAssayData(sobj$Lox2[["RNA"]], slot = "data") ) 

# Add cluster identity to cell barcode ids
samples <- names(gene_exp)
type_counts <- names(gene_exp$Lox2)
for ( i in samples ) {
  cell_ids <- row.names(sobj[[i]]@meta.data)
  seurat_clusters <- as.character(sobj[[i]]@meta.data$seurat_clusters)
  cluster_cell <- paste0(paste0("Clt_", seurat_clusters, "_"), cell_ids)
  for ( x in type_counts ) {
    cell_barcodes <- colnames(gene_exp[[i]][[x]])
    stopifnot( all( cell_barcodes == cell_ids ) )
    colnames(gene_exp[[i]][[x]]) <- cluster_cell
  }
}

# Export matrices RNA counts and data (normalized)
results_folder <- "../results/export_int_mtx_11_06_21"
if ( ! dir.exists(results_folder) ) dir.create(results_folder)
for ( i in samples ) {
  for ( x in type_counts ) {
    write.table(x = cbind("Geneid" = row.names(gene_exp[[i]][[x]]), gene_exp[[i]][[x]]), 
                file = paste(results_folder, paste0(i, "_", x, "_gene_expression_matrix.tsv"), sep = "/"), 
                quote = FALSE, sep = "\t", row.names = FALSE)
  }
}

# Print R session
sessionInfo()
# > sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.5 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Matrix_1.3-2       SeuratObject_4.0.0 Seurat_4.0.0       dplyr_0.8.5       
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-152         matrixStats_0.57.0   RcppAnnoy_0.0.18     RColorBrewer_1.1-2  
# [5] httr_1.4.2           sctransform_0.3.2    tools_3.6.3          R6_2.5.0            
# [9] irlba_2.3.3          rpart_4.1-15         KernSmooth_2.23-20   uwot_0.1.10         
# [13] mgcv_1.8-35          lazyeval_0.2.2       colorspace_1.4-1     tidyselect_1.1.0    
# [17] gridExtra_2.3        compiler_3.6.3       plotly_4.9.3         scales_1.1.1        
# [21] lmtest_0.9-38        spatstat.data_1.7-0  ggridges_0.5.3       pbapply_1.4-3       
# [25] spatstat_1.64-1      goftest_1.2-2        stringr_1.4.0        digest_0.6.27       
# [29] spatstat.utils_2.0-0 pkgconfig_2.0.3      htmltools_0.5.1.1    parallelly_1.22.0   
# [33] fastmap_1.1.0        htmlwidgets_1.5.3    rlang_0.4.10         rstudioapi_0.11     
# [37] shiny_1.6.0          zoo_1.8-8            jsonlite_1.7.2       ica_1.0-2           
# [41] magrittr_2.0.1       patchwork_1.1.1      Rcpp_1.0.6           munsell_0.5.0       
# [45] abind_1.4-5          reticulate_1.18      lifecycle_0.2.0      stringi_1.5.3       
# [49] MASS_7.3-54          Rtsne_0.15           plyr_1.8.6           grid_3.6.3          
# [53] parallel_3.6.3       listenv_0.8.0        promises_1.1.1       ggrepel_0.9.1       
# [57] crayon_1.3.4         miniUI_0.1.1.1       deldir_0.2-9         lattice_0.20-44     
# [61] cowplot_1.1.1        splines_3.6.3        tensor_1.5           pillar_1.4.7        
# [65] igraph_1.2.6         future.apply_1.7.0   reshape2_1.4.4       codetools_0.2-18    
# [69] leiden_0.3.7         glue_1.4.2           data.table_1.13.2    vctrs_0.3.6         
# [73] png_0.1-7            httpuv_1.5.5         gtable_0.3.0         RANN_2.6.1          
# [77] purrr_0.3.4          polyclip_1.10-0      tidyr_1.1.2          scattermore_0.7     
# [81] future_1.21.0        assertthat_0.2.1     ggplot2_3.3.2        mime_0.9            
# [85] xtable_1.8-4         later_1.1.0.1        survival_3.2-11      viridisLite_0.3.0   
# [89] tibble_3.0.6         cluster_2.1.2        globals_0.14.0       fitdistrplus_1.1-3  
# [93] ellipsis_0.3.1       ROCR_1.0-11   
