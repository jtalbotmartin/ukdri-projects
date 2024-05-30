file <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/selective_enrichment/enrichment_res_gazestani/Micro_res_dream_aBeta.tsv"

dir_path <-"/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/selective_enrichment/enrichment_res_gazestani"

# List all files in the directory
all_files <- list.files(dir_path, full.names = TRUE)

# Select only files ending with "aBeta.tsv"
aBeta_files <- all_files[grepl("aBeta\\.tsv$", all_files)]

celltype_output_dfs <- list()

for (file in aBeta_files){
  
  celltype <- strsplit(basename(file), "_")[[1]][[1]]
  dream_data <- read.table(file, header=TRUE, sep = "\t")
  dream_data <- dream_data |>
    mutate("celltype" = celltype)
  celltype_output_dfs <- c(celltype_output_dfs, list(dream_data))
}

combined_celltype_df <- do.call(rbind, celltype_output_dfs)

heatmap_data <- combined_celltype_df |>
  select(c(logFC, celltype)) |>
  as.matrix()

heatmap_data <- reshape2::dcast(heatmap_data, celltype ~ gene, value.var = "logFC")
rownames(heatmap_data) <- heatmap_data$celltype
heatmap_data <- heatmap_data[, -1]  # Remove the 'celltype' column for the heatmap

heatmap(t(heatmap_data))
pheatmap::pheatmap(heatmap_data, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, 
                   color = colorRampPalette(c("blue", "white", "red"))(50))
