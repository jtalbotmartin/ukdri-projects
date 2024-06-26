library(dplyr)
library(ggplot2)

source("/Volumes/jmm17/home/ukdri_projects_code/toolkit/workflow/aucell_dream/gsea_vis.R")
enrichment_res_path <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/selective_enrichment/enrichment_res"

# return base object for downstream visualisation
combined_celltype_df <- process_gsea_enrichment(enrichment_res_path)

# for each variant, filter object and plot facet

CV <- combined_celltype_df |>
  filter(trem2 == "CV") 
plot_facet_trem2var_segregated(CV, "Celltype GSEA Enrichment Scores for TREM2 CV by Contrast")

R47H <- combined_celltype_df |>
  filter(trem2 == "R47H") 
plot_facet_trem2var_segregated(R47H, "Celltype GSEA Enrichment Scores for TREM2 R47H by Contrast")

R62H <- combined_celltype_df |>
  filter(trem2 == "R62H") 
plot_facet_trem2var_segregated(R62H, "Celltype GSEA Enrichment Scores for TREM2 R62H by Contrast")
