#-------------------------------------------------------------------------------
# Description: Several helper functions 

# Author: António Sousa - UBI-IGC (e-mail: agsousa@igc.gulbenkian.pt)

# Date: 22/04/2021
# Last update: 23/04/2021
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
## Function to plot a volcano plot based on the table of DGE obtained with 
#'FindMarkers()' 
#
volcano_plot <- function(dge_tbl, log2FC_cutoff = 1, padj_cutoff = 0.05, n_top = NULL) {
  
  # 'volcano_plot()': plot volcano plot based on 
  #a table of differentially expressed genes 
  #obtained with the function 'FindMarkers()' 
  #from Seurat across two conditions/samples. 
  # 'dge_tbl' (mandatory): differentially 
  #expressed gene table obtained with the 
  #function 'FindMarkers()'. 
  # 'log2FC_cutoff' (mandatory): absolute 
  #log2 fold change cutoff to apply to the 
  #data and highlight genes as differentially 
  #expressed in the volcano plot. Default is 
  # 1 meaning that only genes with a log2 
  #fold change higher or lower than 1 are 
  #considered differentially expressed. Numeric 
  #of length one. 
  # 'padj_cutoff' (mandatory): adjusted p-value
  #cutoff to apply to the data and highlight genes 
  #as differentially expressed in the volcano plot. 
  #Default is 0.05 meaning that only genes with a 
  #adjusted p-value lower than 0.05 are considered 
  #differentially expressed. Numeric of length one. 
  # 'n_top' (optional): by default is 'NULL' 
  #meaning that any gene considered differentially 
  #expressed will be labeled by 'Geneid' - rownames. 
  #If given an integer of length one, it will be 
  #highlighted only these genes ranked from the 
  #highest and lowest log2 fold change for up- and 
  #downregulated genes. If given an integer of length
  #2, the first element will be applied for up- and 
  #the second element for downregulated genes. 
  
  ## required packages
  #require("dplyr")
  #require("ggplot2")
  #require("cowplot")
  #require("ggrepel")
  
  ## Check
  stopifnot(is.data.frame(dge_tbl))
  stopifnot( all( c("avg_log2FC", "p_val_adj") %in% colnames(dge_tbl) ) )
  stopifnot( (is.numeric(log2FC_cutoff) & length(log2FC_cutoff)==1) & (is.numeric(padj_cutoff) & length(padj_cutoff)==1) )
  stopifnot( is.null(n_top) | (is.integer(as.integer(n_top)) & length(n_top)<=2) )
  
  ## Parse data
  dge_parse <- dge_tbl 
  dge_parse[,"Geneid"] <- rownames(dge_tbl)
  dge_parse <- dge_parse %>% 
    #filter(p_val_adj < 0.05) %>% 
    arrange(desc(avg_log2FC), p_val_adj) %>% 
    mutate("Expression" = ifelse((avg_log2FC > log2FC_cutoff) & (p_val_adj < padj_cutoff), "Upregulated",
                                 ifelse((avg_log2FC < (- log2FC_cutoff)) & (p_val_adj < padj_cutoff), "Downregulated", "non_DE"))) %>%
    mutate("Expression" = factor(Expression, levels = c("Upregulated", "Downregulated", "non_DE"))) %>%
    mutate("most_logFC" = ifelse( (avg_log2FC > log2FC_cutoff | avg_log2FC < (- log2FC_cutoff)) & (p_val_adj < padj_cutoff) , Geneid, NA ) )
  
  ## Calculate how many up, down & non-DE
  dge_freq <- table(dge_parse$Expression)
  up <- dge_freq[["Upregulated"]]
  down <- dge_freq[["Downregulated"]]
  non_DE <- dge_freq[["non_DE"]]
  
  ## Select up and down to show label
  dge_parse[,"Label"] <- dge_parse$most_logFC
  n_up_genes <- n_down_genes <- NULL
  if ( ! is.null(n_top) ) {
    if ( length(n_top) == 1 ) {
      n_up <- n_down <- n_top
    } else {
      n_up <- n_top[1]; n_down <- n_top[2];
    }
    if ( up != 0 ) {
      if ( n_up == 0 ) { # do not label
        n_up_genes <- NA
      } else { # label n_up
        n_up_genes <- dge_parse %>% 
          filter(Expression == "Upregulated") %>% 
          arrange(desc(avg_log2FC)) %>% 
          pull(most_logFC) %>% 
          head(n_up)
      } 
    } 
    if ( down != 0 ) {
      if( n_down == 0 ) { # do not label
        n_down_genes = NA
      } else { # label n_down
        n_down_genes <- dge_parse %>% 
          filter(Expression == "Downregulated") %>% 
          arrange(avg_log2FC) %>% 
          pull(most_logFC) %>% 
          head(n_down)
      }
    }  
  }
  if( any( c(!is.null(n_up_genes), !is.null(n_down_genes)) )) {
    genes_2_highlight <- c(n_up_genes, n_down_genes)
    dge_parse <- dge_parse %>% 
      mutate(Label = ifelse(most_logFC %in% genes_2_highlight, as.character(most_logFC), NA))
  }
  
  ## Plot
  colors_plot <- c("#E64B35FF", "#4DBBD5FF", "grey") # colors to use
  volcano <- ggplot(data = dge_parse, 
                    aes(x = avg_log2FC, y = -log10(p_val_adj), label = Label)) + 
    geom_point( 
      aes(fill = Expression), size = 1.5, shape = 21, stroke = 0.25, 
      alpha = ifelse(!is.na(dge_parse$Label), 1, 0.5),
      color = ifelse(!is.na(dge_parse$Label), "black",
                     ifelse(dge_parse$Expression == "Upregulated", "#E64B35FF",
                            ifelse(dge_parse$Expression == "Downregulated", "#4DBBD5FF",
                                   "grey")))
      ) +
    scale_fill_manual(values = colors_plot,
                      name = "", breaks = c("Upregulated", "Downregulated", "non_DE"),
                      labels = c(
                        paste0("Up (>", log2FC_cutoff, " = ", up, ")"), 
                        paste0("Down (<-", log2FC_cutoff, " = ", down, ")"), 
                        paste0("non_DE (>=-", log2FC_cutoff," & <=", log2FC_cutoff," = ", non_DE, ")") 
                        )) +
    cowplot::theme_cowplot() +
    ggrepel::geom_text_repel(size = 3, min.segment.length = 0, seed = 42,
                             box.padding = 0.5, max.overlaps = Inf,
                             segment.curvature = -0.1, nudge_y = 1,
                             nudge_x = 1, segment.ncp = 1, segment.angle = 20,
                             segment.size = 0.2) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype="dotted",
               color = c("#0000ff"), size=.25) +
    geom_vline(xintercept = c(- log2FC_cutoff, log2FC_cutoff), linetype="dotted",
               color = c("dodgerblue4"), size=.25) +
    theme(axis.title = element_text(size = 12, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          text = element_text(size = 9, color = "black"),
          legend.position = c(0.01, 0.9)) +
    ylab("Adjusted p-value (-log10)") +
    xlab("Fold expression (log2)")
    
  return(volcano)
}
#-------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------#
#   
# Plot pseudotime with Slingshot package using the seurat obj as input

slingShot_pseudoTime <- function( seuratObj, starts = NULL, ends = NULL, dim_var = "UMAP" ) {
  
  # seuratObj: seurat obj (S4 class) with 'seurat_clusters'
  # starts: cluster to start with
  # ends: cluster to end with
  # dim_var: one of "UMAP", "t-SNE", "PC" (default: "UMAP") 
  
  #-----------------------------------------------------------------------------------------
  #
  ## code retrieve from: https://bustools.github.io/BUS_notebooks_R/slingshot.html
  require("scales")
  cell_pal <- function(cell_vars, pal_fun,...) {
    if (is.numeric(cell_vars)) {
      pal <- pal_fun(100, ...)
      return(pal[cut(cell_vars, breaks = 100)])
    } else {
      categories <- sort(unique(cell_vars))
      pal <- setNames(pal_fun(length(categories), ...), categories)
      return(pal[cell_vars])
    }
  }
  #
  #---------------------------------------------------------------------------------------
  
  var_dim_1 <- paste0( dim_var, "_1" ); var_dim_2 <- paste0( dim_var, "_2" );
  dim_data <- FetchData( seuratObj, 
                         vars = c( var_dim_1, var_dim_2 ) ) # retrieve dimensional data from Seurat
  colNames <- colnames( seuratObj@meta.data ) # col names of sobj
  clts <- seuratObj@meta.data[ , grep( "res", colNames )] # retrieve cluster names
  slshot_obj <- slingshot( data = dim_data, 
                           clusterLabels = clts, 
                           start.clus = starts, end.clus = ends ) # run slingshot
  liNe <- getLineages( data = dim_data, clusterLabels = clts, 
                       start.clus = starts, end.clus = ends ) # add lineage line
  colorCells <- cell_pal( cell_vars = clts, 
                          pal_fun = hue_pal() ) # add col pal 
  plot( dim_data, col = colorCells, pch = 19, asp = 1, 
        cex = 0.5 ); lines( liNe, lwd = 1, type = 'lineages', col = 'black', 
                            show.constraints = TRUE ); 
  plotOut <- recordPlot()
  
  return( plotOut )
}
#                                                                                                                    #
#--------------------------------------------------------------------------------------------------------------------#
#   
# Plot clusters by Principal Components after performing SlingShot

plot_cells_by_slgshot_pcs <- function( seuratObj, group_var = NULL, 
                                       ord_var = NULL, dim_var = "umap", 
                                       dim_axis = "UMAP_1" ) {
  
  ## Analysis adapted from: https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html
  
  # seuratObj: seurat obj (S4 class) with 'seurat_clusters'
  # ord_var: order the levels of the factor variable highlighted in the y-axis
  # dim_var: one of "umap", "tsne", "pca" (default: "UMAP") 
  # dim_axis: one of the two PCs axis, i.e., if `dim_var = "umap"`, "UMAP_1" or "UMAP_2"
  
  require("gam")
  require("slingshot")
  require("Seurat")
  require("ggbeeswarm")
  require("dplyr")
  
  if ( is.null( group_var ) ) { # give the ident of sobj to the "group_var"
    seuratObj$group_var <- Idents( object = seuratObj )
  } 
  
  label_var <- seuratObj[[ group_var, drop = TRUE ]]
  
  if ( ( !is.null( group_var ) ) & ( !is.null( ord_var ) ) ) { # order labels of "group_var"
    label_var <- factor( label_var, levels = ord_var )
  }
  
  # import dimensions of Sobj into Slingshot obj
  slgShot_obj <- slingshot( Embeddings( object = seuratObj,
                                        reduction = "umap" ), 
                            clusterLabels = label_var )
  
  ## Combine "group_var" with "dim_var"
  #
  dim_df <- slgShot_obj@reducedDim[ , dim_axis] %>% 
    as.data.frame(.) %>%
    mutate( "Cellid" = names( slgShot_obj@reducedDim[, dim_axis] ) ) %>% 
    `colnames<-`( c( dim_axis, "Cell_id" ) )
  dim_df <- dim_df[ , c(2, 1) ] %>% 
    arrange( .[[2]] )
  
  group_df <- label_var %>% 
    as.data.frame(.) %>%
    mutate( "Cellid" = names( label_var ) ) %>% 
    `colnames<-`( c("Cells", "Cell_id") )
  group_df <- group_df[ , c(2, 1) ]
  
  df <- left_join( x = dim_df, y = group_df, 
                   by = "Cell_id" ) # merged tbl
  
  ## Plot cell groups by PCs
  #
  plot_cells_by_pcs <- ggplot( df, aes(x = UMAP_1, y = Cells, colour = Cells)) +
    ggbeeswarm::geom_quasirandom( groupOnX = FALSE ) + 
    theme_classic()
  
  return(plot_cells_by_pcs)
}
#                                                                                                                    #
#--------------------------------------------------------------------------------------------------------------------#
#   
# Plot clusters by Principal Components from Seurat

plot_cells_by_seurat_pcs <- function( seuratObj, group_var = NULL, 
                                      ord_var = NULL, dim_var = "umap", 
                                      dim_axis = "UMAP_1" ) {
  
  ## Analysis adapted from: https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html
  
  # seuratObj: seurat obj (S4 class) with 'seurat_clusters'
  # ord_var: order the levels of the factor variable highlighted in the y-axis
  # dim_var: one of "umap", "tsne", "pca" (default: "UMAP") 
  # dim_axis: one of the two PCs axis, i.e., if `dim_var = "umap"`, "UMAP_1" or "UMAP_2"
  
  require("gam")
  require("slingshot")
  require("Seurat")
  require("ggbeeswarm")
  require("dplyr")
  
  if ( is.null( group_var ) ) { # give the ident of sobj to the "group_var"
    seuratObj$group_var <- Idents( object = seuratObj )
  } 
  
  label_var <- seuratObj[[ group_var, drop = TRUE ]]
  
  if ( ( !is.null( group_var ) ) & ( !is.null( ord_var ) ) ) { # order labels of "group_var"
    label_var <- factor( label_var, levels = ord_var )
  }
  
  ## Combine "group_var" with "dim_var"
  #
  coord <- seuratObj@reductions[[ dim_var ]][[, dim_axis, drop = TRUE]]
  dim_df <- coord %>% 
    as.data.frame(.) %>%
    mutate( "Cellid" = names( coord ) ) %>% 
    `colnames<-`( c( dim_axis, "Cell_id" ) )
  dim_df <- dim_df[ , c(2, 1) ] %>% 
    arrange( .[[2]] )
  
  group_df <- label_var %>% 
    as.data.frame(.) %>%
    mutate( "Cellid" = names( label_var ) ) %>% 
    `colnames<-`( c("Cells", "Cell_id") )
  group_df <- group_df[ , c(2, 1) ]
  
  df <- left_join( x = dim_df, y = group_df, 
                   by = "Cell_id" ) # merged tbl
  colnames( df ) <- c( "Cell_id", "Dim", "Cells" )
  
  ## Plot cell groups by PCs
  #
  plot_cells_by_pcs <- ggplot( df, aes(x = Dim, y = Cells, colour = Cells)) +
    ggbeeswarm::geom_quasirandom( groupOnX = FALSE ) + 
    xlab( dim_axis ) +
    theme_classic()
  
  return(plot_cells_by_pcs)
}
#                                                                                                                    #
#--------------------------------------------------------------------------------------------------------------------#
#   
# Plot diffusion map and clusters by DC's 

sobjs2difusionMap <- function( seuratObj, exp_mtx = "integrated", 
                               label_cells = "seurat_clusters", 
                               cell_ord = NULL, 
                               dim_axis = "DC1",
                               seed = 2020 ) {
  
  # seuratObj: seurat obj (S4 class) with 'seurat_clusters'
  # label_cells: factor variable to group cells
  # cell_ord: order the levels of the factor variable highlighted in the y-axis, i.e., "label_cells"
  # dim_axis: one of the two PCs axis, i.e., if `dim_var = "umap"`, "UMAP_1" or "UMAP_2"
  
  set.seed( seed = seed )
  
  require( "SingleCellExperiment" )
  require( "destiny" )
  require( "Seurat" )
  require( "dplyr" )
  require( "ggbeeswarm" )
  
  sc_intObj <- as.SingleCellExperiment( x = seuratObj, assay = exp_mtx )
  log_mtx <- logcounts( object = sc_intObj )
  label_cells <- seuratObj[[ label_cells, drop = TRUE ]]
  colnames(log_mtx) <- label_cells # change barcode id cell labels to "label_cells"
  
  # run diffusion map 
  dm <- DiffusionMap( t( as.matrix( log_mtx ) ) )
  
  # df to plot the two first DC (Diffusion Components)
  df_1 <- data.frame( "DC1" = eigenvectors( dm )[ , 1],
                      "DC2" = eigenvectors( dm )[ , 2],
                      "Cells" = colnames( log_mtx ) )
  
  # plot
  plot_1 <- 
    ggplot( df_1, aes( x = DC1, y = DC2, colour = Cells) ) +
    geom_point() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic()
  
  #-----------------------------------------------------------------------
  
  ## Plot 
  
  sc_intObj$dm_eigenvalues <- rank( eigenvectors( dm )[ , dim_axis ]) 
  
  df_2 <- colData( sc_intObj ) %>% 
    as.data.frame(.)
  
  if ( !is.null( cell_ord ) ) {
    df_2[ , label_cells, drop = TRUE ] <- factor( df_2[ , label_cells, drop = TRUE ], 
                                                  levels = cell_ord )
  }
  
  # plot 
  
  plot_2 <- 
    ggplot( df_2, aes( x = dm_eigenvalues, y = label_cells, color = label_cells ) ) +
    geom_quasirandom(groupOnX = FALSE) + 
    theme_classic() + 
    xlab( paste0( "Diffusion component ", gsub( "DC", "", dim_axis ), 
                  " ", dim_axis ) ) + 
    ylab("Cells")
  
  plots <- list( "diffusion_map" = plot_1, 
                 "Cell_groups_by_dm" = plot_2 )
  
  return( plots )
}
#--------------------------------------------------------------------------------------------------------------------#
### Function to run pseudotime from Seurat objs and target gene list
# ### Perform pseudotime analysis

#This pseudotime analysis is based on the code and solution presented at: 
#https://github.com/satijalab/seurat/issues/1658

pseudoTimeIntegrated <- function(seuratObj, targetGeneList = NULL) {
  
  # seuratObj: seurat obj (S4 class) with 'seurat_clusters'
  # targetGeneList: vector with a gene list to perform  
  #pseudotime
  
  ## Import packages
  require("monocle")
  require("Seurat")
  require("dplyr")
  
  # check statements
  stopifnot(is(seuratObj, "Seurat"))
  stopifnot("seurat_clusters" %in% colnames(seuratObj@meta.data))
  
  # gene mtx, phenotype and feature data from the SeuratObject
  geneMtx <- as(as.matrix(seuratObj@assays$integrated@data), 
                'sparseMatrix') # gene mtx
  phenoDat <- new('AnnotatedDataFrame', 
                  data = seuratObj@meta.data) # phenotype data
  featureDat <- data.frame(gene_short_name = row.names(geneMtx), 
                           row.names = row.names(geneMtx))
  featureDat <- new('AnnotatedDataFrame', data = featureDat) # feature data
  
  # import data into monocle
  monocleObj <- newCellDataSet(cellData = geneMtx, 
                               phenoData = phenoDat, 
                               featureData = featureDat, 
                               expressionFamily = uninormal() # 'uninormal()': 
                               #because gene expression 
                               #is already normalized (from Seurat)
  )
  
  ## Define target genes to perform pseudotime analysis
  #
  if ( is.null(targetGeneList) ) {
    targetGenes <- seuratObj@assays$integrated@var.features 
  } else {
    targetGenes <- targetGeneList
  }
  
  # add/order target genes info to monocle obj 
  monocleObj <- setOrderingFilter(monocleObj, targetGenes)  
  
  ## Reduce dims
  monocleObj <- reduceDimension(monocleObj, 
                                norm_method = "none", 
                                reduction_method = "DDRTree",
                                max_components = 2,
                                scaling = TRUE,
                                verbose = FALSE,
                                pseudo_expr = 0)
  
  ## Order cells
  monocleObj <- orderCells(monocleObj)
  
  ## Plot trajectories
  plotPseudotime <- plot_cell_trajectory(monocleObj, 
                                         color_by = "seurat_clusters", 
                                         cell_size = .5) +
    theme(legend.position = "right")
  
  rm(list = c("monocleObj", "geneMtx", 
              "targetGenes", "phenoDat", 
              "featureDat")) # remove data from env space
  return(plotPseudotime)
}
#--------------------------------------------------------------------------------------------------------------------#