library(tidyverse)
library(scFlow)
library(SummarizedExperiment)
library(SingleCellExperiment)

sce_path <- "/data/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/results/final/SCE/final_sce"
save_path <- "/data/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/objects/filter_sce"

sce <- read_sce(sce_path) # super slow when running on local machine
# sce <- readRDS("/data/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/results/final/sce_rds_obj/sce.rds")

cat(paste0("Study IDs present in the sce object pre-filtering ", unique(sce$StudyID, "\n")))

cat("Filtering SCE by StudyID")

sce <- sce[ , sce$StudyID %in% c("IGFQ000883", "IGFQ000955", "IGFQ001110")] # filter to only initial cohort
sce$StudyID <- droplevels(sce$StudyID)

cat(paste0("Study IDs present in the sce object post-filtering ", unique(sce$StudyID, "\n")))

# Save the sce into .qs obj for quicker loading
qs::qsave(sce, file = file.path(save_path, "sce.qs"))

