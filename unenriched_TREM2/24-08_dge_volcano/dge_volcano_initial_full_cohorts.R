library(dplyr)
library(ggplot2)

source("/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/additional_scripts/plot_volcano.r")


plot_save_volcano_celltypes_in_model <- function(base_dir, plots_dir, model){
  
  model_plots_dir <- paste0(plots_dir, model)
  dir.create(model_plots_dir)
  
  celltypes <-dir(paste0(base_dir, model))
  
  for( celltype in celltypes){
    
    model_celltype_plots_dir <- paste0(model_plots_dir, "/", celltype)
    dir.create(model_celltype_plots_dir)
    
    file_path <- sprintf("%s%s/%s", base_dir, model, celltype)
    file_list <- list.files(path = file_path, pattern =".tsv" , full.names = TRUE) #[-4]
    
    for (file in file_list){
      plot_name <- strsplit(basename(file), split = ".tsv")[[1]]
      # file <- read.table(file_list[[1]], header = TRUE, sep = "\t")
      file <- read.table(file, header = TRUE, sep = "\t")
      volcano_plot <- volcano_plot(file)
      ggsave(filename = paste0(plot_name, ".png"), plot = volcano_plot, path = paste0(model_celltype_plots_dir))
    }
  }
}

### generate plots for full cohort

base_dir <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/full_cohort/DGE/"
plots_dir <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/volcano_plots/full_cohort/"

#getting models in directory
model_full_cohort <- dir(base_dir, pattern = "^de_")

for (model in model_full_cohort){
  plot_save_volcano_celltypes_in_model(base_dir, plots_dir, model)
}

### generate plots for initial cohort

base_dir <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/DGE/"
plots_dir <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/volcano_plots/initial_cohort/"

#getting models in directory
model_initial_cohort <- dir(base_dir, pattern = "^de_")

for (model in model_initial_cohort){
  plot_save_volcano_celltypes_in_model(base_dir, plots_dir, model)
}


