dir_path <-"/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/selective_enrichment/enrichment_res"
file <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/selective_enrichment/enrichment_res_gazestani/Micro_res_dream_aBeta.tsv"

# List all files in the directory

all_files <- list.files(dir_path, full.names = TRUE, recursive = TRUE)

# # Select only files ending with "aBeta.tsv"
# aBeta_files <- all_files[grepl("aBeta\\.tsv$", all_files)]

celltype_output_dfs <- list()

for (file in all_files){
  
  celltype <- strsplit(basename(file), "_")[[1]][[1]]
  trem2var <- strsplit(basename(file), "_")[[1]][[5]] |>
    strsplit(".tsv")
  contrast <- strsplit(basename(file), "_")[[1]][[4]]
  dream_data <- read.table(file, header=TRUE, sep = "\t")
  dream_data <- dream_data |>
    mutate("celltype"     = celltype) |>
    mutate("trem2"        = trem2var) |>
    mutate("significance" = if_else(adj.P.Val <= 0.05, TRUE, FALSE)) |>
    mutate("contrast_truncated"     = contrast)
  
  celltype_output_dfs <- c(celltype_output_dfs, list(dream_data))
}

combined_celltype_df <- do.call(rbind, celltype_output_dfs)
# truncate name to make visualisation more readable
combined_celltype_df <- combined_celltype_df |>
  mutate("contrast_truncated" = case_when(
    contrast_truncated == "NeuropathologicalDiagnosis" ~ "NPD",
    contrast_truncated == "pct4G8PositiveArea"         ~ "4G8",
    contrast_truncated == "pctPHF1PositiveArea"        ~ "PHF1"))

# heatmap_data <- combined_celltype_df |>
#   select(c(logFC, celltype)) |>
#   as.matrix()
# 
# heatmap_data <- reshape2::dcast(heatmap_data, celltype ~ gene, value.var = "logFC")
# rownames(heatmap_data) <- heatmap_data$celltype
# heatmap_data <- heatmap_data[, -1]  # Remove the 'celltype' column for the heatmap
# 
# heatmap(t(heatmap_data))
# pheatmap::pheatmap(heatmap_data, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, 
#                    color = colorRampPalette(c("blue", "white", "red"))(50))


CV_df <- combined_celltype_df |>
  filter(trem2 == "CV")

ggplot(CV_df, aes(x = celltype, y = geneset, fill = z.std, color = significance )) +
  geom_point(aes(size = abs(z.std)), shape = 21) +
  xlab("") +
  ylab("") +
  scale_x_discrete(position = "top") +
  facet_grid(contrast_truncated ~ trem) +
  scale_fill_distiller(name = "Z-score", 
                       palette = "RdYlBu", direction = -1) +
  scale_size(range = c(2, 7)) +
  guides(size = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 14,
                                   angle = 90, hjust = 0, vjust = 0),
        axis.text.y = element_text(colour = "black", size = 16), 
  ) 


