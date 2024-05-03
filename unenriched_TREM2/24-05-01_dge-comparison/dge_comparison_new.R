
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


dge_celltype_statistics <- function(model_path){
  celltypes <- dir(model_path)
  res_all <- list()
  # outdir <- "bar_plot_DEG"
  # dir.create(outdir)
  
  for( celltype in celltypes){
    
    file_path <- sprintf("%s/%s", model_path, celltype)
    
    file_list <- list.files(path = file_path, pattern =".tsv" , full.names = TRUE) #[-4]
    
    dt_l <- lapply(file_list, read.delim)
    names(dt_l) <- basename(file_list) %>% tools::file_path_sans_ext()
    
    res_l <- lapply(names(dt_l), function(x){
      
      de_up <- dt_l[[x]] %>%
        dplyr::filter(padj <= 0.05, logFC >= 0.25) %>%
        dplyr::pull(gene) %>%
        as.character() %>%
        length()
      
      de_down <- dt_l[[x]] %>%
        dplyr::filter(padj <= 0.05, logFC <= -0.25) %>%
        dplyr::pull(gene) %>%
        as.character() %>%
        length()
      
      
      tmp_dt <- data.frame(contrast = x,
                           total_expressed_gene = nrow(dt_l[[x]]),
                           up = de_up,
                           down = de_down,
                           pct_up = de_up*100/nrow(dt_l[[x]]),
                           pct_down = de_down*100/nrow(dt_l[[x]]),
                           model= basename(model_path))
    })
    
    res <- do.call(rbind, res_l)
    
    res_all[[celltype]] <- res
    
  }
  
  res <- do.call(rbind, res_all)
  res <- res %>%
    mutate(celltype = purrr::map_chr(contrast, ~strsplit(., "_")[[1]][1]),
           trem2 = purrr::map_chr(contrast, ~strsplit(., "_")[[1]][4]), #from 5
           pct_down = -1*pct_down)
  
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

#####################

#---

prev_ADvControl <- dge_celltype_statistics("previous_results/de_results_diagnosis_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group_age_PMD")
new_ADvControl  <- dge_celltype_statistics("ADvControl/de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")

compare_ADvControl <- rbind(prev_ADvControl, new_ADvControl)

#---

prev_PHF1     <- dge_celltype_statistics("previous_results/de_results_TREM2_stratified_PHF1_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group")
new_PHF1      <- dge_celltype_statistics("PHF1byTREM2/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
new_PHF1_noBR <- dge_celltype_statistics("PHF1byTREM2/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group")


compare_PHF1_BR <- rbind(prev_PHF1, new_PHF1)
compare_PHF1_noBR  <- rbind(prev_PHF1, new_PHF1_noBR)
compare_PHF1_models  <- rbind(new_PHF1, new_PHF1_noBR)

# plots for new 

plot_DE_prop_celltypes(new_ADvControl, "ADvControl")
plot_DE_prop_celltypes(new_PHF1, "PHF1")
plot_DE_prop_celltypes(new_PHF1_noBR, "PHF1")


