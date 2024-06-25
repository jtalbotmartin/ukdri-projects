# 3/6/24 running IPA analysis from DEGs of full cohort for Biogen (generated using common geneset)

library(dplyr)
library(ggplot2)
library(tidyr)
library(scFlow)

setwd("/Volumes/jmm17/home/ukdri_projects_code/unenriched_TREM2/24-05_IPA")
source("../../toolkit/R/ipa.R")
source("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/DGE_scripts/pathway_enrichment_enrichR.r")

DGE_tables_dir <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/full_cohort/DGE/"


# original models

NPD1_list    <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
ABeta_1_list <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
PHF1_1_list  <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group"))
PHF1_2_list  <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))

# no CD33 

NPD_2_list    <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))
ABeta_2_list <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))
PHF1_3_list  <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))
# PHF1_2_list  <- extract_tables_to_celltype_list(paste0(DGE_tables_dir, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))

setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/ipa/full_cohort/Pathway_Analysis_CommonSet")

# for (celltype in 1:length(list)){
#   dir.create("pathway_")
#   dir_path = paste0("pathway_<>/", names(list)[[celltype]])
#   generate_ipa_celltype_separate_updown(list[[celltype]], dir_path)
# }

# original runs


for (celltype in 1:length(NPD1_list)){
  dir.create("pathway_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
  dir_path = paste0("pathway_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group/", names(NPD1_list)[[celltype]])
  generate_ipa_celltype_separate_updown(NPD1_list[[celltype]], dir_path)
}

for (celltype in 1:length(ABeta_1_list)){
  dir.create("pathway_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
  dir_path = paste0("pathway_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group/", names(ABeta_1_list)[[celltype]])
  generate_ipa_celltype_separate_updown(ABeta_1_list[[celltype]], dir_path)
}

for (celltype in 1:length(PHF1_1_list)){
  dir.create("pathway_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group")
  dir_path = paste0("pathway_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group/", names(PHF1_1_list)[[celltype]])
  generate_ipa_celltype_separate_updown(PHF1_1_list[[celltype]], dir_path)
}

for (celltype in 1:length(PHF1_2_list)){
  dir.create("pathway_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
  dir_path = paste0("pathway_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group/", names(PHF1_2_list)[[celltype]])
  generate_ipa_celltype_separate_updown(PHF1_2_list[[celltype]], dir_path)
}


# no CD33


for (celltype in 1:length(NPD_2_list)){
  dir.create("pathway_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")
  dir_path = paste0("pathway_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup/", names(NPD_2_list)[[celltype]])
  generate_ipa_celltype_separate_updown(NPD_2_list[[celltype]], dir_path)
}

for (celltype in 1:length(ABeta_2_list)){
  dir.create("pathway_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")
  dir_path = paste0("pathway_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup/", names(ABeta_2_list)[[celltype]])
  generate_ipa_celltype_separate_updown(ABeta_2_list[[celltype]], dir_path)
}

for (celltype in 1:length(PHF1_3_list)){
  dir.create("pathway_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")
  dir_path = paste0("pathway_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup/", names(PHF1_3_list)[[celltype]])
  generate_ipa_celltype_separate_updown(PHF1_3_list[[celltype]], dir_path)
}

