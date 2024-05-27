
# - for each file, calculate proportion up/down regulated
# - store in table with AD / PHF1, model, filename, %up/down
# - grouped pathway analysis plot?

library(dplyr)
library(ggplot2)
library(tidyr)

setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/")


# directories:

ADvControl_dirs <-c("previous_results/de_results_diagnosis_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group_age_PMD",
                    "ADvControl/de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
  
PHF1_dirs  <-     c("previous_results/de_results_TREM2_stratified_PHF1_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group",
                    "PHF1byTREM2/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group",
                    "PHF1byTREM2/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")


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

# compress Ex and In Neurons together

compress_neuronal_celltypes <- function(dataframe){
 
  compressed_df <- dataframe %>%
    mutate(celltype = ifelse(substr(celltype, 1, 2) == "EN", "EN", celltype)) %>%
    mutate(celltype = ifelse(substr(celltype, 1, 2) == "IN", "IN", celltype)) %>%
    group_by(celltype, trem2, pct_de) %>%                   # Group by trem2 and pct_de columns
    summarise(
      total_expressed_gene = sum(total_expressed_gene),
      up = sum(up),
      down = sum(down)
    ) %>%
    mutate(
      value = ifelse(pct_de == "pct_up", (up / total_expressed_gene) * 100, (down / total_expressed_gene) * 100),
      label = paste0(round(value, 2), "%"),
      label_y = ifelse(value > 0, value + 1, value - 1)
    ) %>%
    ungroup()
  
  return(compressed_df)
  
}

compress_vasc_celltypes <- function(dataframe){
  
  compressed_df <- dataframe %>%
    mutate(celltype = ifelse(celltype == "Endo", "Vasc", celltype)) %>%
    mutate(celltype = ifelse(celltype == "FB", "Vasc", celltype)) %>%
    mutate(celltype = ifelse(celltype == "SMC", "Vasc", celltype)) %>%
    group_by(celltype, trem2, pct_de) %>%                   # Group by trem2 and pct_de columns
    summarise(
      total_expressed_gene = sum(total_expressed_gene),
      up = sum(up),
      down = sum(down)
    ) %>%
    mutate(
      value = ifelse(pct_de == "pct_up", (up / total_expressed_gene) * 100, (down / total_expressed_gene) * -100),
      label = paste0(round(value, 2), "%"),
      label_y = ifelse(value > 0, value + 1, value - 1)
    ) %>%
    ungroup()
  
  return(compressed_df)
  
}

append_comparison_label_to_df <- function(df, comparison_label){
 df <- df |>
    mutate("comparison_label" = comparison_label)
  return(df)
}

process_prev_df <- function(df){
  df <- df |>
    select("celltype", "trem2", "pct_de", "total_expressed_gene", "up", "down", "value", "label", "label_y", "comparison_label") |>
    filter(trem2 != "TREM2var")
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
      axis.text = ggplot2::element_text(size = 16, colour = "black"),
      axis.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 16),
      legend.title = ggplot2::element_text(size = 16),
      strip.text = ggplot2::element_text(size = 16)
    ) +
    geom_text(aes(label = label, y = label_y)) +
    labs(x = x_label, y = y_label, 
         title = graph_title)
  
}

########################################

# get plots of 

<model> <-
  dge_celltype_statistics(<dir>)
plot_DE_prop_celltypes(<model>, <label>)


########################################

#---

# generate statistics between previous and new results for ADvControl
prev_ADvControl <- dge_celltype_statistics("previous_results/de_results_diagnosis_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group_age_PMD")
prev_ADvControl <- append_run_label(prev_ADvControl, "0")
prev_ADvControl <- process_prev_df(prev_ADvControl) # doing this for processing sake 


new_ADvControl  <- dge_celltype_statistics("ADvControl/de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
test_ADvControl <- new_ADvControl
new_ADvControl <- compress_neuronal_celltypes(new_ADvControl)
new_ADvControl <- compress_vasc_celltypes(new_ADvControl)
new_ADvControl <- append_run_label(new_ADvControl, "1")

compare_ADvControl <- rbind(prev_ADvControl, new_ADvControl)
compare_ADvControl <- compare_ADvControl |>
  mutate(TREM2Compare = paste0(trem2, ":", comparison_label))

plot_comparison(compare_ADvControl, "ADvControl Comparison w/ Previous Run")

# sanity check
# > new_ADvControl |>
#   +     filter(celltype == "EN") |>
#   +     pull(total_expressed_gene) |>
#   +     sum()


##################################################


#---

prev_PHF1 <- dge_celltype_statistics("previous_results/de_results_TREM2_stratified_PHF1_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group")
prev_PHF1 <- process_prev_df(prev_PHF1) # doing this for processing sake 
prev_PHF1 <- compress_neuronal_celltypes(prev_PHF1)
prev_PHF1 <- append_run_label(prev_PHF1, "0")



new_PHF1 <- dge_celltype_statistics("PHF1byTREM2/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
new_PHF1 <- compress_neuronal_celltypes(new_PHF1)
new_PHF1 <- compress_vasc_celltypes(new_PHF1)
new_PHF1 <- append_run_label(new_PHF1, "1")


new_PHF1_noBR <- dge_celltype_statistics("PHF1byTREM2/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group")
new_PHF1_noBR <- compress_neuronal_celltypes(new_PHF1_noBR)
new_PHF1_noBR <- compress_vasc_celltypes(new_PHF1_noBR)
new_PHF1_noBR <- append_run_label(new_PHF1_noBR, "1")

compare_PHF1 <- rbind(prev_PHF1, new_PHF1)
compare_PHF1 <- compare_PHF1 |>
  mutate(TREM2Compare = paste0(trem2, ":", comparison_label))

plot_comparison(compare_PHF1, "PHF1byTREM2 Comparison w/ Previous Run")

compare_PHF1_noBR <- rbind(prev_PHF1, new_PHF1_noBR)
compare_PHF1_noBR <- compare_PHF1 |>
  mutate(TREM2Compare = paste0(trem2, ":", comparison_label))

plot_comparison(compare_PHF1, "PHF1byTREM2 Comparison w/ Previous Run (No Brain Region in Model)")

#---

ABeta <- dge_celltype_statistics("ABbyTREM2/de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group")


# plots for new 

plot_DE_prop_celltypes(new_ADvControl, "ADvControl")
plot_DE_prop_celltypes(new_PHF1, "PHF1")
plot_DE_prop_celltypes(new_PHF1_noBR, "PHF1")
plot_DE_prop_celltypes(ABeta, "Aβ")










# # Assuming your dataframe is called 'df'
# compressed_df <- new_ADvControl %>%
#   mutate(celltype = ifelse(substr(celltype, 1, 3) == "Exc", "Exc", celltype)) %>%
#   mutate(celltype = ifelse(substr(celltype, 1, 3) == "Inh", "Inh", celltype)) %>%
#   group_by(trem2, pct_de) %>%                   # Group by trem2 and pct_de columns
#   summarise_at(vars(total_expressed_gene, up, down), sum)  # Sum the values in specified columns

# If you want to keep the original rows along with the compressed rows, you can use bind_rows
final_df <- bind_rows(compressed_df, df %>% filter(substr(celltype, 1, 3) != "Exc"))



