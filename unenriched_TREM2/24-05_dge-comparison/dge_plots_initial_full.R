# 3/6/24: 
# plotting up/down % DGE for intitial and full cohorts after recalculation with common gene set
# plotting comparison figures for initial vs full cohort for each contrast

library(dplyr)
library(ggplot2)
library(tidyr)

setwd("/Volumes/jmm17/home/ukdri_projects_code/unenriched_TREM2/24-05_dge-comparison")
source("../../toolkit/R/dge_comparison.R")
setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/")

#################### initial cohort ##########################################


initial_cohort_base_path <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/DGE/"

NPD1_1_init <- dge_celltype_statistics(paste0(initial_cohort_base_path, 
                                  "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))

ABeta_1_init <- dge_celltype_statistics(paste0(initial_cohort_base_path, 
                                  "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))                                  

PHF1_1_init <- dge_celltype_statistics(paste0(initial_cohort_base_path, 
                                  "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group"))                                  
                                  
PHF1_2_init <- dge_celltype_statistics(paste0(initial_cohort_base_path, 
                                  "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))                                                                                                      



plot_DE_prop_celltypes(NPD1_1_init, "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
plot_DE_prop_celltypes(ABeta_1_init, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
plot_DE_prop_celltypes(PHF1_1_init, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group")
plot_DE_prop_celltypes(PHF1_2_init, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")

#################### full cohort ##########################################

full_cohort_base_path <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/full_cohort/DGE/"

NPD1_1_full <- dge_celltype_statistics(paste0(full_cohort_base_path, 
                                         "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))

ABeta_1_full <- dge_celltype_statistics(paste0(full_cohort_base_path, 
                                          "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))                                  

PHF1_1_full <- dge_celltype_statistics(paste0(full_cohort_base_path, 
                                         "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group"))                                  

PHF1_2_full <- dge_celltype_statistics(paste0(full_cohort_base_path, 
                                         "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))                                                                                                      



plot_DE_prop_celltypes(NPD1_1_full, "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
plot_DE_prop_celltypes(ABeta_1_full, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")
plot_DE_prop_celltypes(PHF1_1_full, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group")
plot_DE_prop_celltypes(PHF1_2_full, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")

################ initial:full cohort comparison#################################

# NPD_1

NPD1_1_init <- append_comparison_label_to_df(NPD1_1_init, "0")
NPD1_1_full <- append_comparison_label_to_df(NPD1_1_full, "1")

compare_NPD1_1 <- rbind(NPD1_1_init, NPD1_1_full)
compare_NPD1_1 <- compare_NPD1_1 |>
  mutate(TREM2Compare = paste0(trem2, ":", comparison_label))

plot_comparison(compare_NPD1_1, "NPDbyTREM2 Comparison of Initial and Full Cohorts - de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")

# 4G8_1

ABeta_1_init <- append_comparison_label_to_df(ABeta_1_init, "0")
ABeta_1_full <- append_comparison_label_to_df(ABeta_1_full, "1")

compare_ABeta_1 <- rbind(ABeta_1_init, ABeta_1_full)
compare_ABeta_1 <- compare_ABeta_1 |>
  mutate(TREM2Compare = paste0(trem2, ":", comparison_label))

plot_comparison(compare_ABeta_1, "4G8byTREM2 Comparison of Initial and Full Cohorts - de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")

# PHF1_1

PHF1_1_init <- append_comparison_label_to_df(PHF1_1_init, "0")
PHF1_1_full <- append_comparison_label_to_df(PHF1_1_full, "1")

compare_PHF1_1 <- rbind(PHF1_1_init, PHF1_1_full)
compare_PHF1_1 <- compare_PHF1_1 |>
  mutate(TREM2Compare = paste0(trem2, ":", comparison_label))

plot_comparison(compare_PHF1_1, "PHF1byTREM2 Comparison of Initial and Full Cohorts - de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group")

# PHF1_2

PHF1_2_init <- append_comparison_label_to_df(PHF1_2_init, "0")
PHF1_2_full <- append_comparison_label_to_df(PHF1_2_full, "1")

compare_PHF1_2 <- rbind(PHF1_2_init, PHF1_2_full)
compare_PHF1_2 <- compare_PHF1_2 |>
  mutate(TREM2Compare = paste0(trem2, ":", comparison_label))

plot_comparison(compare_PHF1_2, "PHF1byTREM2 Comparison of Initial and Full Cohorts - de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group")






