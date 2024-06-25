# 12/6/24 generating correlation plots for dge results 
# Can you please generate some correlation plot for the regressions in the initial cohort, with CD33 and without CD33 models?
# So this will be 4G8 and PHF1, the full model (w CD33) and full model (w/o CD33) for all celltypes.

library(dplyr)
library(ggplot2)
library(tidyr)

setwd("/Volumes/jmm17/home/ukdri_projects_code/unenriched_TREM2/24-06_dge-correlation")
source("../../toolkit/R/ipa.R")
source("../../toolkit/R/dge_correlation.R")
source("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/DGE_scripts/pathway_enrichment_enrichR.r")

DGE_tables_dir <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/DGE/"

# original models

NPD_1_list    <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
ABeta_1_list <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
# PHF1_1_list  <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group"))
PHF1_2_list  <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))

# no CD33

NPD_2_list    <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))
ABeta_2_list <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))
PHF1_3_list  <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))
# PHF1_2_list  <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))

setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/correlation_figures")

# for (celltype in 1:length(list)){
#   dir.create("pathway_")
#   dir_path = paste0("pathway_<>/", names(list)[[celltype]])
#   generate_ipa_celltype_separate_updown(list[[celltype]], dir_path)
# }

for (celltype in 1:length(NPD_1_list)){

  if (names(NPD_1_list)[celltype] == names(NPD_2_list)[celltype]){
    dir.create("NPD_comparison", showWarnings = FALSE)
    dir_path = paste0("NPD_comparison/", names(NPD_1_list)[[celltype]])
    dir.create(dir_path)
    dge_correlations(NPD_1_list[[celltype]], NPD_2_list[[celltype]], "NeuropathologicalDiagnosis", dir_path, "wCD33", "woCD33")
  }
}

for (celltype in 1:length(ABeta_1_list)){
  
  if (names(ABeta_1_list)[celltype] == names(ABeta_2_list)[celltype]){
    dir.create("4G8_comparison", showWarnings = FALSE)
    dir_path = paste0("4G8_comparison/", names(ABeta_1_list)[[celltype]])
    dir.create(dir_path)
    dge_correlations(ABeta_1_list[[celltype]], ABeta_2_list[[celltype]], "pct4G8PositiveArea", dir_path, "wCD33", "woCD33")
  }
}

for (celltype in 1:length(PHF1_2_list)){
  
  if (names(PHF1_2_list)[celltype] == names(PHF1_3_list)[celltype]){
    dir.create("PHF1_comparison", showWarnings = FALSE)
    dir_path = paste0("PHF1_comparison/", names(PHF1_2_list)[[celltype]])
    dir.create(dir_path)
    dge_correlations(PHF1_2_list[[celltype]], PHF1_3_list[[celltype]], "pctPHF1PositiveArea", dir_path, "wCD33", "woCD33")
  }
}