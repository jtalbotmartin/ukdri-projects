
# - for each file, calculate proportion up/down regulated
# - store in table with AD / PHF1, model, filename, %up/down
# - grouped pathway analysis plot?

# Function to calculate differential gene expression statistics for each cell type
# Input: model_path - directory path containing subdirectories for each cell type, each containing TSV files with gene expression data
# Output: A dataframe containing statistics for each contrast (pairwise comparison) of gene expression within each cell type


dge_celltype_statistics <- function(model_path){
  # Get list of subdirectories (cell types) within the model_path directory
  celltypes <- dir(model_path)
  # Initialize an empty list to store results for all cell types
  res_all <- list()

  # Iterate over each cell type directory
  for( celltype in celltypes){

    # Construct full file path for the current cell type
    file_path <- sprintf("%s/%s", model_path, celltype)

    # Get list of TSV files in the current cell type directory
    file_list <- list.files(path = file_path, pattern =".tsv" , full.names = TRUE) #[-4]

    # Read each TSV file into a list of data frames
    dt_l <- lapply(file_list, read.delim)

    # Assign names to the list elements based on file names (without extension)
    names(dt_l) <- basename(file_list) %>% tools::file_path_sans_ext()

    # Calculate statistics for each contrast within the cell type
    res_l <- lapply(names(dt_l), function(x){

      # Count genes with significant upregulation according to padj and logFC thresholds
      de_up <- dt_l[[x]] %>%
        dplyr::filter(padj <= 0.05, logFC >= 0.25) %>%
        dplyr::pull(gene) %>%
        as.character() %>%
        length()

      # Count genes with significant downregulation according to padj and logFC thresholds
      de_down <- dt_l[[x]] %>%
        dplyr::filter(padj <= 0.05, logFC <= -0.25) %>%
        dplyr::pull(gene) %>%
        as.character() %>%
        length()

      # Create a dataframe with contrast statistics
      tmp_dt <- data.frame(contrast = x,
                           total_expressed_gene = nrow(dt_l[[x]]),
                           up = de_up,
                           down = de_down,
                           pct_up = de_up*100/nrow(dt_l[[x]]),
                           pct_down = de_down*100/nrow(dt_l[[x]]),
                           model= basename(model_path))
    })

    # Combine the list of dataframes into a single dataframe for the current cell type
    res <- do.call(rbind, res_l)

    # Store the results for the current cell type in the res_all list
    res_all[[celltype]] <- res

  }
  # Combine results for all cell types into a single dataframe
  res <- do.call(rbind, res_all)

  # Extract additional information from contrast names
  res <- res %>%
    mutate(celltype = purrr::map_chr(contrast, ~strsplit(., "_")[[1]][1]),
           trem2 = purrr::map_chr(contrast, ~strsplit(., "_")[[1]][4]), #from 5
           pct_down = -1*pct_down)

  # Reshape the dataframe + create labels for visualization

  res_l <- pivot_longer(res,
                        cols = c(pct_up,pct_down),
                        names_to = "pct_de",
                        values_to = "value"
  ) %>%
    mutate(label = ifelse(value > 0, paste0(round(value, 2), "%"), paste0(-1* round(value, 2), "%")),
           label_y = ifelse(value > 0, value + 1, value -1)
    )

  return(res_l)

}

#################

# Function to plot proportion of differentially expressed genes (DEGs) for each cell type
# Input:
#   - graph_data: dataframe containing data to be plotted, with columns 'trem2', 'value', 'celltype', 'label', and 'label_y'
#   - graph_title: title of the plot

plot_DE_prop_celltypes <- function (graph_data, graph_title){

  palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")


  ggplot(graph_data, aes(x = trem2, y = value, fill = trem2)) +
    geom_bar(stat = "identity", position = position_dodge(0.8)) +
    facet_wrap(~ celltype) +
    geom_hline(yintercept = 0) +
    scale_color_manual(name = "TREM2", values = palette_choice,
                       aesthetics = c("colour", "fill")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 16, colour = "black"),
      axis.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 16),
      legend.title = ggplot2::element_text(size = 16),
      strip.text = ggplot2::element_text(size = 16)
    ) +
    geom_text(aes(label = label, y = label_y)) +
    labs(x = "TREM2_Var", y = "% Up/Down DE Genes",
         title = graph_title)

}


plot_comparison <- function (graph_data, graph_title){

  palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")


  ggplot(graph_data, aes(x = TREM2Compare, y = value, fill = trem2)) +
    geom_bar(stat = "identity", position = position_dodge(0.8)) +
    facet_wrap(~ celltype) +
    geom_hline(yintercept = 0) +
    scale_color_manual(name = "TREM2", values = palette_choice,
                       aesthetics = c("colour", "fill")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 16, colour = "black"),
      axis.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 16),
      legend.title = ggplot2::element_text(size = 16),
      strip.text = ggplot2::element_text(size = 16)
    ) +
    geom_text(aes(label = label, y = label_y)) +
    labs(x = "TREM2_Var (0 = Previous Run, 1 = New Run)", y = "% Up/Down DE Genes",
         title = graph_title)

}

append_comparison_label_to_df <- function(df, comparison_label){
  df <- df |>
    mutate("comparison_label" = comparison_label)
  return(df)
}



plot_num_DEG_celltype_comparison <- function (graph_data, graph_title, x_label, y_label){

  palette_choice <- paletteer::paletteer_d("ggsci::nrc_npg")


  ggplot(graph_data, aes(x = TREM2Compare, y = value, fill = trem2)) +
    geom_bar(stat = "identity", position = position_dodge(0.8)) +
    facet_wrap(~ celltype) +
    geom_hline(yintercept = 0) +
    scale_color_manual(name = "TREM2", values = palette_choice,
                       aesthetics = c("colour", "fill")) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 14, colour = "black", angle = 90, vjust = 0.5, hjust=1),
      axis.title = ggplot2::element_text(size = 14, margin = margin(t = 20)),
      legend.text = ggplot2::element_text(size = 14),
      legend.title = ggplot2::element_text(size = 14),
      strip.text = ggplot2::element_text(size = 18)
    ) +
    geom_text(aes(label = label, y = label_y)) +
    labs(x = x_label, y = y_label,
         title = graph_title)

}

combine_dfs_for_comparison_by_TREM2 <- function(df1, df2){
  combined_dfs <- rbind(df1, df2)
  combined_dfs <- combined_dfs |>
    filter(trem2 %in% c("CV", "R47H", "R62H")) |>
    mutate(TREM2Compare = paste0(trem2, ":", comparison_label))
  return(combined_dfs)
}

rename_ENIN_ExcInh <- function(df){
    df <- df %>%
      mutate(celltype = ifelse(substr(celltype, 1, 2) == "EN", paste0("Exc", substr(celltype, 3, nchar(celltype))), celltype)) %>%
      mutate(celltype = ifelse(substr(celltype, 1, 2) == "IN", paste0("Inh", substr(celltype, 3, nchar(celltype))), celltype))
    return(df)


}

