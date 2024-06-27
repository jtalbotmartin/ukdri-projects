##########################################################################
# function to plot individual TREM2 variant enrichment scores
##########################################################################

plot_facet_trem2var_segregated <- function(trem2var_stratified_df, title_label){

    # plot geneset enrichment scores, by each celltype, faceted by each contrast
    # colour by z-score, outline by significance

    ggplot(trem2var_stratified_df, aes(x = contrast_truncated, y = geneset, fill = z.std, stroke = ifelse(significance, 0.75, 0))) +
    geom_point(aes(size = abs(z.std)), shape = 21) +
    xlab("") +
    ylab("") +
    ggtitle(title_label) +
    scale_x_discrete(position = "bottom") +
    facet_grid(cols = vars(celltype)) +
    scale_fill_distiller(name = "Z-score", 
                        palette = "RdYlBu", direction = -1) +
    scale_size(range = c(2, 7)) +
    guides(size = "none", color = "none") +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black", size = 11, angle = 90, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(colour = "black", size = 11), 
            strip.text.x = element_text(size = 12, angle = 45), 
            strip.background = element_blank(),
            panel.spacing = unit(0.1, "lines")
  )
}

##########################################################################
# function to process enrichment scores from GSEA analysis for use in 
# downstream plotting
##########################################################################

process_gsea_enrichment <- function(enrichment_results_toplevel_dir){

    # List all files in the directory
    all_files <- list.files(enrichment_results_toplevel_dir, full.names = TRUE, recursive = TRUE)

    celltype_output_dfs <- list()

    # process all files in enrichment result top-level directory (+ its subdirs)
    # extract key information from filename, read and transform data for visualisation

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

    # combine dfs for downstream mutation and filtering
    combined_celltype_df <- do.call(rbind, celltype_output_dfs)

    # truncate name to make visualisation more readable
    combined_celltype_df <- combined_celltype_df |>
    mutate("contrast_truncated" = case_when(
        contrast_truncated == "NeuropathologicalDiagnosis" ~ "NPD",
        contrast_truncated == "pct4G8PositiveArea"         ~ "4G8",
        contrast_truncated == "pctPHF1PositiveArea"        ~ "PHF1"))

    # return unfiltered object for visualisation
    return(combined_celltype_df)
}