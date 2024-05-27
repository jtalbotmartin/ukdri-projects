# up and down separately
# combine results from database



dt <- read.delim("Micro_pct4G8PositiveArea_in_CV.tsv") %>%
  mutate(gene = as.character(gene))
contrast <- "CV"

for(direction in c("up", "down")){
  
  if(direction == "up"){
    
    de_gene <- dt %>%
      filter(padj <= 0.05, logFC >= 0.25, gene_biotype == "protein_coding") %>%
      pull(gene)
    
  } else {
    de_gene <- dt %>%
      filter(padj <= 0.05, logFC <= -0.25, gene_biotype == "protein_coding") %>%
      pull(gene) 
  }
  
  res_l <- pathway_analysis_enrichr(gene_file = de_gene,
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
    write.table(res, file = file_name, sep = "\t", row.names = F, quote = F)



