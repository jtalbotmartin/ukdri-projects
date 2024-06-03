
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

generate_ipa_celltype_separate_updown <- function(celltype_list_of_df, dir_path, padj_threshold = 0.05, logfc_threshold = 0.25){
  
  dir.create(dir_path)
  
  for (i in 1:length(celltype_list_of_df)){ 
    
    # split celltype model name, to make directories for each variant
    # outdir <- paste0(dir_path,"/PA_", strsplit(names(celltype_list_of_df)[[i]], "_")[[1]][[4]])
    # dir.create(outdir)
    
    contrast <- "CV"
    
    for(direction in c("up", "down")){
      
      if(direction == "up"){
        
        sig_de_genes <- celltype_list_of_df[[i]] %>%
          dplyr::filter(padj <= 0.05, logFC >= 0.25, gene_biotype == "protein_coding") %>%
          dplyr::pull(gene) %>%
          as.character()
        
        res_l <- pathway_analysis_enrichr(gene_file = sig_de_genes,
                                          enrichment_database = c(
                                            "GO_Molecular_Function_2023",
                                            "GO_Cellular_Component_2023",
                                            "GO_Biological_Process_2023",
                                            "WikiPathways_2021_Human",
                                            "KEGG_2021_Human"),
                                          is_output = FALSE)
        
        
        
        if(length(setdiff(names(res_l), c("plot", "metadata")))> 0){
          
          res_l$plot <- NULL
          res_l$metadata <- NULL
          
          res <- do.call(rbind, res_l)
          res$database <- rownames(res)
          res$database <- gsub(".[0-9]{1,3}$", "", fixed = F, res$database )
          contrast_name <- paste(contrast, direction, sep = "_")
          res$DEGs <- contrast_name
          write.table(res, file = paste0(dir_path, "/UP_", (names(celltype_list_of_df)[[i]]), ".tsv"), sep = "\t", row.names = F, quote = F)
        }
        
      } else if (direction == "down"){
        sig_de_genes <- celltype_list_of_df[[i]] %>%
          dplyr::filter(padj <= 0.05, logFC <= -0.25, gene_biotype == "protein_coding") %>%
          dplyr::pull(gene) %>%
          as.character()
        
        res_l <- pathway_analysis_enrichr(gene_file = sig_de_genes,
                                          enrichment_database = c(
                                            "GO_Molecular_Function_2023",
                                            "GO_Cellular_Component_2023",
                                            "GO_Biological_Process_2023",
                                            "WikiPathways_2021_Human",
                                            "KEGG_2021_Human"),
                                          is_output = FALSE)
        
        
        
        if(length(setdiff(names(res_l), c("plot", "metadata")))> 0){
          
          res_l$plot <- NULL
          res_l$metadata <- NULL
          
          res <- do.call(rbind, res_l)
          res$database <- rownames(res)
          res$database <- gsub(".[0-9]{1,3}$", "", fixed = F, res$database )
          contrast_name <- paste(contrast, direction, sep = "_")
          res$DEGs <- contrast_name
          write.table(res, file = paste0(dir_path, "/DOWN_", (names(celltype_list_of_df)[[i]]), ".tsv"), sep = "\t", row.names = F, quote = F)
                   
      }
    
        # names(celltype_list_of_df)[[i]])
    
      }
    } 
  }
}
##################### generate celltype_lists

DGE_tables_dir <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/Analyses/cellranger_mapping/Biogen_snRNAseq_analysis/TREM2_unenriched_snRNAseq/DGE/"

NPD_1_list <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
ABeta_1_list <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
ABeta_2_list <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_BrainRegion_APOEgroup_CD33Group"))
PHF1_1_list <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group"))
PHF1_2_list <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
PHF1_3_list <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_BrainRegion_APOEgroup_CD33Group"))

######################

####### IPA results ##################

setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/Analyses/cellranger_mapping/Biogen_snRNAseq_analysis/TREM2_unenriched_snRNAseq/Pathway_Analysis")


##### get IPA for each celltype, for each model

# NPD_1_list, "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group/"
# **ABeta_1_list, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group/"
# ABeta_2_list, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_BrainRegion_APOEgroup_CD33Group/"
# PHF1_1_list, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group/"
# **PHF1_2_list, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group/"
# PHF1_3_list, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_BrainRegion_APOEgroup_CD33Group/"


####### Running Baptiste's

NPD <- extract_tables_to_celltype_list("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/Analyses/cellranger_mapping/Biogen_snRNAseq_analysis/TREM2_enriched_snRNAseq/DGE/de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
ABeta <- extract_tables_to_celltype_list("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/Analyses/cellranger_mapping/Biogen_snRNAseq_analysis/TREM2_enriched_snRNAseq/DGE/de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
PHF1 <- extract_tables_to_celltype_list("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/Analyses/cellranger_mapping/Biogen_snRNAseq_analysis/TREM2_enriched_snRNAseq/DGE/de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")


# alter path, and list to generate results for each celltype, after generating sig_deg directory

for (celltype in 1:length(<list>)){
  dir_path = paste0("<path>", names(<list>)[[celltype]])
  generate_ipa_celltype_separate_updown(<list>[[celltype]], dir_path)
}

setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/Analyses/cellranger_mapping/Biogen_snRNAseq_analysis/enriched_updated_pathway_analysis")
source("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/DGE_scripts/pathway_enrichment_enrichR.r")


# setwd("de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
for (celltype in 1:length(PHF1)){
  dir.create("pathway_TREM2Variant_pctPHF1PositiveArea")
  dir_path = paste0("pathway_TREM2Variant_pctPHF1PositiveArea/", names(PHF1)[[celltype]])
  generate_ipa_celltype_separate_updown(PHF1[[celltype]], dir_path)
}





