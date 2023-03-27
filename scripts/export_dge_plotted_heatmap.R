## Export DGE tables used in the heatmaps 
# Run by: António Sousa - UBI-IGC (agsousa@igc.gulbenkian.pt)
# Date: 22/07/2021

## Set seed 
#
set.seed(seed = 1024) # to keep reproducibility

## Import packages
#
library("dplyr")

## Import DGE tables
#
dge_tables <- "../results/int_28_05_21/tables/dge_tables"
dge2import <- list.files(path = dge_tables, pattern = "clt_", full.names = TRUE)
names(dge2import) <- lapply(basename(dge2import), function(x) strsplit(x = x, split = "_")[[1]][1:2] %>% 
                              unlist(.) %>% paste(., collapse = "_")) %>% unlist(.)
dge <- lapply(dge2import, function(x) read.table(file = x, header = TRUE, sep = "\t", 
                                                 stringsAsFactors = FALSE))

## Save tables 
#
table2save <- "../results/export_dge_tbls_used_2_heatmap_22_07_21"
if ( ! dir.exists(table2save) ) dir.create(table2save)

## loop over list of DGE tbls by cluster, parse them and export them
clts <- paste0("clt_", 0:7)
for ( clt in clts ) {
  sub_df <- dge[[clt]] %>% 
    filter(p_val_adj<0.05 & abs(avg_log2FC)>0) %>% 
    arrange(desc(avg_log2FC))
  write.table(x = sub_df, file = paste(table2save, paste0(clt, "_dge_table_p_val_adj_005.tsv"), 
                                       sep = "/"), row.names = FALSE, quote = FALSE, sep = "\t")
} 

















