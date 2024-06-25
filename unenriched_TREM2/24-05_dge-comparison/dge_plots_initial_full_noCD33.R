# 11/6/24: 
# plotting up/down % DGE for intitial and full cohorts after recalculation with common gene set, using reduced model dropping CD33Group
# plotting comparison figures for initial vs full cohort for each contrast

library(dplyr)
library(ggplot2)
library(tidyr)

setwd("/Volumes/jmm17/home/ukdri_projects_code/unenriched_TREM2/24-05_dge-comparison")
source("../../toolkit/R/dge_comparison.R")
setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/")

#################### initial cohort ##########################################


initial_cohort_base_path <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/DGE/"

NPD1_2_init <- dge_celltype_statistics(paste0(initial_cohort_base_path, 
                                  "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))

ABeta_2_init <- dge_celltype_statistics(paste0(initial_cohort_base_path, 
                                  "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))                                  

PHF1_3_init <- dge_celltype_statistics(paste0(initial_cohort_base_path, 
                                  "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))                                  

plot_DE_prop_celltypes(NPD1_2_init, "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")
plot_DE_prop_celltypes(ABeta_2_init, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")
plot_DE_prop_celltypes(PHF1_3_init, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")

#################### full cohort ##########################################

full_cohort_base_path <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/full_cohort/DGE/"

NPD1_2_full <- dge_celltype_statistics(paste0(full_cohort_base_path, 
                                         "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))

ABeta_2_full <- dge_celltype_statistics(paste0(full_cohort_base_path, 
                                          "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))                                  

PHF1_3_full <- dge_celltype_statistics(paste0(full_cohort_base_path, 
                                         "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup"))                                  

plot_DE_prop_celltypes(NPD1_2_full, "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")
plot_DE_prop_celltypes(ABeta_2_full, "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")
plot_DE_prop_celltypes(PHF1_3_full, "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")

################ initial:full cohort comparison#################################

# NPD_2

NPD1_2_init <- append_comparison_label_to_df(NPD1_2_init, "0")
NPD1_2_full <- append_comparison_label_to_df(NPD1_2_full, "1")

compare_NPD1_2 <- rbind(NPD1_2_init, NPD1_2_full)
compare_NPD1_2 <- compare_NPD1_2 |>
  mutate(TREM2Compare = paste0(trem2, ":", comparison_label))

plot_comparison(compare_NPD1_2, "NPDbyTREM2 Comparison of Initial and Full Cohorts - de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")

# 4G8_2

ABeta_2_init <- append_comparison_label_to_df(ABeta_2_init, "0")
ABeta_2_full <- append_comparison_label_to_df(ABeta_2_full, "1")

compare_ABeta_2 <- rbind(ABeta_2_init, ABeta_2_full)
compare_ABeta_2 <- compare_ABeta_2 |>
  mutate(TREM2Compare = paste0(trem2, ":", comparison_label))

plot_comparison(compare_ABeta_2, "4G8byTREM2 Comparison of Initial and Full Cohorts - de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")

# PHF1_3

PHF1_3_init <- append_comparison_label_to_df(PHF1_3_init, "0")
PHF1_3_full <- append_comparison_label_to_df(PHF1_3_full, "1")

compare_PHF1_3 <- rbind(PHF1_3_init, PHF1_3_full)
compare_PHF1_3 <- compare_PHF1_3 |>
  mutate(TREM2Compare = paste0(trem2, ":", comparison_label))

plot_comparison(compare_PHF1_3, "PHF1byTREM2 Comparison of Initial and Full Cohorts - de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup")



