
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