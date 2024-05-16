library(dplyr)
library(ggplot2)
library(tidyr)

source("../../toolkit/R/dge_comparison.R")

setwd("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/")

### read in all DGE results for comparison
##### previous runs of initial cohort

prev_runs_initial_base <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/archive/previous_results/"

prev_NPD  <- dge_celltype_statistics(paste0(prev_runs_initial_base, "de_results_diagnosis_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group_age_PMD"))
#~diagnosis+(1|manifest)+cngeneson+pc_mito+sex+brain_region+apoe+CD33_group + // age // + PMD
prev_4G8  <- dge_celltype_statistics(paste0(prev_runs_initial_base, "de_results_TREM2_stratified_amyloid_beta_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group"))
#~amyloid_beta+(1|manifest)+cngeneson+pc_mito+sex+brain_region+apoe+CD33_group
prev_PHF1 <- dge_celltype_statistics(paste0(prev_runs_initial_base, "de_results_TREM2_stratified_PHF1_cngeneson_pc_mito_sex_brain_region_apoe_CD33_group"))
#~PHF1+(1|manifest)+cngeneson+pc_mito+sex+brain_region+apoe+CD33_group



##### new runs of initial cohort

new_initial_runs_base <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/DGE/"

new_NPD_1_compare <- dge_celltype_statistics(paste0(new_initial_runs_base, 
                                  "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
# ~NeuropathologicalDiagnosis+(1|manifest)+cngeneson+pc_mito+Sex+Age+PostMortemInterval+BrainRegion+APOEgroup+CD33Group
new_4G8_1<- dge_celltype_statistics(paste0(new_initial_runs_base, 
                                  "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
#~pct4G8PositiveArea+(1|manifest)+cngeneson+pc_mito+Sex+Age+PostMortemInterval+BrainRegion+APOEgroup+CD33Group
new_4G8_2_compare<- dge_celltype_statistics(paste0(new_initial_runs_base, 
                                  "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_BrainRegion_APOEgroup_CD33Group"))
#~pct4G8PositiveArea+(1|manifest)+cngeneson+pc_mito+Sex+BrainRegion+APOEgroup+CD33Group
new_PHF1_1<- dge_celltype_statistics(paste0(new_initial_runs_base, 
                                  "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_APOEgroup_CD33Group"))
#~pctPHF1PositiveArea+(1|manifest)+cngeneson+pc_mito+Sex+Age+PostMortemInterval+APOEgroup+CD33Group
new_PHF1_2<- dge_celltype_statistics(paste0(new_initial_runs_base, 
                                  "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
#~pctPHF1PositiveArea+(1|manifest)+cngeneson+pc_mito+Sex+Age+PostMortemInterval+BrainRegion+APOEgroup+CD33Group
new_PHF1_3<- dge_celltype_statistics(paste0(new_initial_runs_base, 
                                  "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_APOEgroup_CD33Group"))
#~pctPHF1PositiveArea+(1|manifest)+cngeneson+pc_mito+Sex+APOEgroup+CD33Group
new_PHF1_4_compare<- dge_celltype_statistics(paste0(new_initial_runs_base, 
                                  "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_BrainRegion_APOEgroup_CD33Group"))
#~pctPHF1PositiveArea+(1|manifest)+cngeneson+pc_mito+Sex+BrainRegion+APOEgroup+CD33Group



##### new runs of full cohort (only for comparison)

new_full_runs_base <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/full_cohort/DGE/"

full_cohort_compare_NPD <- dge_celltype_statistics(paste0(new_full_runs_base, 
                                                         "de_TREM2Variant_NeuropathologicalDiagnosis_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
# ~NeuropathologicalDiagnosis+(1|manifest)+cngeneson+pc_mito+Sex+Age+PostMortemInterval+BrainRegion+APOEgroup+CD33Group
full_cohort_compare_4G8 <- dge_celltype_statistics(paste0(new_full_runs_base, 
                                                         "de_TREM2Variant_pct4G8PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
#~pct4G8PositiveArea+(1|manifest)+cngeneson+pc_mito+Sex+Age+PostMortemInterval+BrainRegion+APOEgroup+CD33Group
full_cohort_compare_PHF1 <- dge_celltype_statistics(paste0(new_full_runs_base, 
                                                           "de_TREM2Variant_pctPHF1PositiveArea_cngeneson_pc_mito_Sex_Age_PostMortemInterval_BrainRegion_APOEgroup_CD33Group"))
#~pctPHF1PositiveArea+(1|manifest)+cngeneson+pc_mito+Sex+Age+PostMortemInterval+BrainRegion+APOEgroup+CD33Group


########################################
########################################
########################################


### process for comparisons

### full vs initial cohort (new runs) comparisons

##### NPD

new_NPD_1_compare       <- append_comparison_label_to_df(new_NPD_1_compare, "0")
full_cohort_compare_NPD <- append_comparison_label_to_df(full_cohort_compare_NPD, "1")

compare_full_initial_NPD <- combine_dfs_for_comparison_by_TREM2(new_NPD_1_compare, full_cohort_compare_NPD)
plot_num_DEG_celltype_comparison(compare_full_initial_NPD, "NPD by TREM2 in Initial vs Full Cohorts", "TREM2_Var (0 =  Initial Cohort, 1 = Full Cohort)", "% Up/Down DE Genes")

##### 4G8

new_4G8_1               <- append_comparison_label_to_df(new_4G8_1, "0")
full_cohort_compare_4G8 <- append_comparison_label_to_df(full_cohort_compare_4G8, "1")

compare_full_initial_4G8 <- combine_dfs_for_comparison_by_TREM2(new_4G8_1, full_cohort_compare_4G8)
plot_num_DEG_celltype_comparison(compare_full_initial_4G8, "4G8 by TREM2 in Initial vs Full Cohorts", "TREM2_Var (0 =  Initial Cohort, 1 = Full Cohort)", "% Up/Down DE Genes")

##### PHF1

new_PHF1_2               <- append_comparison_label_to_df(new_PHF1_2, "0")
full_cohort_compare_PHF1 <- append_comparison_label_to_df(full_cohort_compare_PHF1, "1")

compare_full_initial_PHF1 <- combine_dfs_for_comparison_by_TREM2(new_PHF1_2, full_cohort_compare_PHF1)
plot_num_DEG_celltype_comparison(compare_full_initial_PHF1, "PHF1 by TREM2 in Initial vs Full Cohorts", "TREM2_Var (0 =  Initial Cohort, 1 = Full Cohort)", "% Up/Down DE Genes")


### prev runs vs current run comparison (initial cohort)

##### NPD
## highest-level: remap labels for EN and IN, drop unused celltypes

prev_NPD_basic_comparison <- rename_ENIN_ExcInh(prev_NPD) |>
  filter(celltype %in% c("Astro", "Micro", "Oligo", "OPC")) |>
  append_comparison_label_to_df("0")

new_NPD_1_basic_comparison <- new_NPD_1_compare |>
  filter(celltype %in% c("Astro", "Micro", "Oligo", "OPC")) |>
  append_comparison_label_to_df("1")

basic_compare_4G8 <- combine_dfs_for_comparison_by_TREM2(prev_NPD_basic_comparison, new_NPD_1_basic_comparison)
plot_num_DEG_celltype_comparison(basic_compare_4G8, "NPD by TREM2 in Initial Cohort: Previous and Updated Run ", "TREM2_Var (0 = Prev, 1 = Updated)", "% Up/Down DE Genes")

## compression of celltypes: remap + recalculate for Vasc and neuronal


##### 4G8
# note - needed to rename all files from amyloid_beta to abeta to ensure trem2 captured
# c("Astro", "Micro", "Oligo", "OPC", "EN-L2-3", "EN-L3-5", "EN-L4-6", "IN-LAMP5", "IN-PVALB", "IN-VIP")

prev_4G8_basic_comparison <- rename_ENIN_ExcInh(prev_4G8) |>
  mutate(celltype = gsub("^Exc-L4-6-I$", "Exc-L4-6", celltype)) |>
  filter(celltype %in% c("Astro", "Micro", "Oligo", "OPC", "Exc-L2-3", "Exc-L3-5", "Exc-L4-6", "Inh-LAMP5", "Inh-PVALB", "Inh-VIP")) |>
  append_comparison_label_to_df("0")

new_4G8_2_basic_comparison <- new_4G8_2_compare |>
  mutate(celltype = gsub("^Exc-L2-3-I$", "Exc-L2-3", celltype)) |>
  filter(celltype %in% c("Astro", "Micro", "Oligo", "OPC", "Exc-L2-3", "Exc-L3-5", "Exc-L4-6", "Inh-LAMP5", "Inh-PVALB", "Inh-VIP")) |>
  append_comparison_label_to_df("1")

basic_compare_4G8 <- combine_dfs_for_comparison_by_TREM2(prev_4G8_basic_comparison, new_4G8_2_basic_comparison)
plot_num_DEG_celltype_comparison(basic_compare_4G8, "4G8 by TREM2 in Initial Cohort: Previous and Updated Run ", "TREM2_Var (0 = Prev, 1 = Updated)", "% Up/Down DE Genes")


##### PHF1

prev_PHF1_basic_comparison <- rename_ENIN_ExcInh(prev_PHF1) |>
  mutate(celltype = gsub("^Exc-L4-6-I$", "Exc-L4-6", celltype)) |>
  filter(celltype %in% c("Astro", "Micro", "Oligo", "OPC", "Exc-L2-3", "Exc-L3-5", "Exc-L4-6", "Inh-LAMP5", "Inh-PVALB", "Inh-VIP")) |>
  append_comparison_label_to_df("0")

new_PHF1_4_basic_comparison <- new_PHF1_4_compare |>
  mutate(celltype = gsub("^Exc-L2-3-I$", "Exc-L2-3", celltype)) |>
  filter(celltype %in% c("Astro", "Micro", "Oligo", "OPC", "Exc-L2-3", "Exc-L3-5", "Exc-L4-6", "Inh-LAMP5", "Inh-PVALB", "Inh-VIP")) |>
  append_comparison_label_to_df("1")

basic_compare_PHF1 <- combine_dfs_for_comparison_by_TREM2(prev_PHF1_basic_comparison, new_PHF1_4_basic_comparison)
plot_num_DEG_celltype_comparison(basic_compare_PHF1, "PHF1 by TREM2 in Initial Cohort: Previous and Updated Run ", "TREM2_Var (0 = Prev, 1 = Updated)", "% Up/Down DE Genes")
