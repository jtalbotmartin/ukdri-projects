
# - for each file, calculate proportion up/down regulated
# - store in table with AD / PHF1, model, filename, %up/down
# - grouped pathway analysis plot?

library(dplyr)
library(ggplot2)
library(tidyr)
library(scFlow)

setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/")


# directories:

ADvControl_dirs <-c("previous_results/de_results_diagnosis_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group_age_PMD",
                    "ADvControl/de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
  
PHF1_dirs  <-     c("previous_results/de_results_TREM2_stratified_PHF1_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group",
                    "PHF1byTREM2/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group",
                    "PHF1byTREM2/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")



extract_tables_to_celltype_list <- function(model_path){
  
  celltypes <- dir(model_path)

  # Initialize an empty list to store results for each celltype
  res_all <- list()
  
  list_res <- list()
  
  for( celltype in celltypes){
    
    file_path <- sprintf("%s/%s", model_path, celltype)
    
    # List all files with the extension ".tsv" within the current celltype directory
    file_list <- list.files(path = file_path, pattern =".tsv" , full.names = TRUE) #[-4]
    
    # Read each file into a data frame and store it in a list, with filenames as list names
    dt_l <- lapply(file_list, read.delim)
    names(dt_l) <- basename(file_list) %>% tools::file_path_sans_ext()
    
    # Combine the data frames into a single data frame and store it under the current celltype
    # res <- do.call(rbind, dt_l) 
    
    res_all[[celltype]] <- dt_l
    
  }
  
  # Combine all data frames stored in res_all into a single data frame
  # res <- do.call(rbind, res_all)
  return(res_all)
}

extract_tables_to_aggregated_df <- function(model_path){
  
  celltypes <- dir(model_path)
  
  # Initialize an empty list to store results for each celltype
  res_all <- list()
  
  list_res <-
    
    for( celltype in celltypes){
      
      file_path <- sprintf("%s/%s", model_path, celltype)
      
      # List all files with the extension ".tsv" within the current celltype directory
      file_list <- list.files(path = file_path, pattern =".tsv" , full.names = TRUE) #[-4]
      
      # Read each file into a data frame and store it in a list, with filenames as list names
      dt_l <- lapply(file_list, read.delim)
      names(dt_l) <- basename(file_list) %>% tools::file_path_sans_ext()
      
      # Combine the data frames into a single data frame and store it under the current celltype
      res <- do.call(rbind, dt_l)
      res_all[[celltype]] <- res
      
    }
  
  # Combine all data frames stored in res_all into a single data frame
  res <- do.call(rbind, res_all)
  return(res)
}

###########

# input one celltype, containing 3 stratified dfs - generates file with list of 
# sig deg genes for that celltype

get_sig_deg_celltype <- function(celltype_list_of_df, dir_path, padj_threshold = 0.05, logfc_threshold = 0.25){
  
  dir.create(dir_path)
  
  for (i in 1:length(celltype_list_of_df)){ 
    
    filename <- paste0("de_", names(celltype_list_of_df)[[i]], ".tsv")
    
    sig_de_genes <- celltype_list_of_df[[i]] %>%
      dplyr::filter(padj <= padj_threshold, abs(logFC) >= logfc_threshold) %>%
      dplyr::pull(gene) %>%
      as.character()
    
    sig_de_genes <- data.frame(value = unlist(sig_de_genes))
    
    file_path <- filename
    
    write.table(sig_de_genes, file = paste0(dir_path, "/", filename), sep = "\t", row.names = FALSE, col.names = FALSE)
    
  }
  
  # return(sig_de_genes)
}


# input one celltype, containing 3 stratified dfs - generates file outputs from
# IPA analysis

generate_ipa_celltype <- function(celltype_list_of_df, dir_path, padj_threshold = 0.05, logfc_threshold = 0.25){
  
  dir.create(dir_path)
  
  for (i in 1:length(celltype_list_of_df)){ 
    
    dir.create(paste0(dir_path,"/ipa_", names(celltype_list_of_df)[[i]]))

    sig_de_genes <- celltype_list_of_df[[i]] %>%
      dplyr::filter(padj <= padj_threshold, abs(logFC) >= logfc_threshold) %>%
      dplyr::pull(gene) %>%
      as.character()

    enrichment_result <- pathway_analysis_enrichr(
      sig_de_genes,
      enrichment_database = c("GO_Biological_Process_2021",
                              "KEGG_2021_Human"))

    write.table(enrichment_result$GO_Biological_Process_2021, file = paste0(dir_path,"/ipa_", names(celltype_list_of_df)[[i]],"/GO_Biological_Process_2021.tsv"), sep = "\t", row.names = FALSE)
    write.table(enrichment_result$KEGG_2021_Human, file = paste0(dir_path,"/ipa_", names(celltype_list_of_df)[[i]],"/KEGG_2021_Human.tsv"), sep = "\t", row.names = FALSE)

    if ("GO_Biological_Process_2021" %in% names(enrichment_result$plot)) {
      png(paste0(dir_path,"/ipa_", names(celltype_list_of_df)[[i]], "/GO_Biological_Process_2021.png"))
      plot(enrichment_result$plot$GO_Biological_Process_2021)
      dev.off()
    }

    if ("KEGG_2021_Human" %in% names(enrichment_result$plot)) {
      png(paste0(dir_path,"/ipa_", names(celltype_list_of_df)[[i]], "/KEGG_2021_Human.png"))
      plot(enrichment_result$plot$KEGG_2021_Human)
      dev.off()
    }

    metadata <- data.frame(value = unlist(enrichment_result$metadata))
    write.table(metadata, file = paste0(dir_path,"/ipa_", names(celltype_list_of_df)[[i]], "/metadata.tsv"), sep = "\t", row.names = TRUE, col.names = FALSE)

    qs::qsave(enrichment_result, file = paste0(dir_path,"/ipa_", names(celltype_list_of_df)[[i]], "/ipa.qs"))
    
  }
}




ADvControl_list  <- extract_tables_to_celltype_list("ADvControl/de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
PHF1_list        <- extract_tables_to_celltype_list("PHF1byTREM2/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
PHF1_noBR_list   <- extract_tables_to_celltype_list("PHF1byTREM2/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group")

# getting tables of DEGs for each celltype

setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/")


##### get DEGs for each celltype, stratified by TREM2Var

# alter path, and list to generate results for each celltype, after generating sig_deg directory

for (celltype in 1:length(PHF1_noBR_list)){
  dir_path = paste0("PHF1byTREM2/sig_deg/", names(PHF1_noBR_list)[[celltype]])
  get_sig_deg_celltype(PHF1_noBR_list[[celltype]], dir_path, 0.05, 0.25)
}


####### IPA results ##################

setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/ipa/")


##### get DEGs for each celltype, stratified by TREM2Var

# alter path, and list to generate results for each celltype, after generating sig_deg directory

for (celltype in 1:length(ADvControl_list)){
  dir_path = paste0("ADvControl/", names(ADvControl_list)[[celltype]])
  generate_ipa_celltype(ADvControl_list[[celltype]], dir_path, 0.05, 0.25)
}






