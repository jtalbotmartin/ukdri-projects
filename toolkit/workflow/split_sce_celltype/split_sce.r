library(tidyverse)
library(scFlow)
library(SummarizedExperiment)
library(SingleCellExperiment)

sce_path <- "/data/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/results/final/SCE/final_sce"

outdir <- "/data/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/objects/split_sce_celltype"

sce <- read_sce(sce_path)
# sce <- readRDS("/data/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/results/final/sce_rds_obj/sce.rds")

sce <- sce[ , sce$StudyID %in% c("IGFQ000883", "IGFQ000955", "IGFQ001110")] # filter to only initial cohort
sce$StudyID <- droplevels(sce$StudyID)

sce <- sce[!is.na(SummarizedExperiment::rowData(sce)$gene), ]
sce <- sce[!duplicated(SummarizedExperiment::rowData(sce)$gene), ]

sce$TREM2Variant <- as.character(sce$TREM2Variant)
sce$TREM2Variant <- ifelse(sce$TREM2Variant == "None", "CV", sce$TREM2Variant)
sce$TREM2Variant <- as.factor(sce$TREM2Variant)

sce$APOEgroup <- case_when(grepl("E4", sce$APOE) ~ "APOE4-pos", TRUE ~ "APOE4-neg")
sce$APOEgroup <- as.factor(sce$APOEgroup)

# Save the sce into .qs obj for quicker loading
qs::qsave(sce, file = file.path(outdir, "sce.qs"))

# sce <- qs::qread(file.path(outdir, "sce.qs"))

idx <- unique(sce$cluster_celltype)

dir.create(file.path(outdir, "celltype_sce"))

# Optional: If celltype clusters below a certain number needs to be excluded 
# idx <- names(table(sce$cluster_celltype)[(table(sce$cluster_celltype) >= 1000)])

for (celltype in idx) {
  sce_subset <- sce[, sce$cluster_celltype == celltype]
  sce_name <- file.path(outdir, "celltype_sce", paste0(celltype, "_sce.qs"))
  qs::qsave(sce_subset, file = sce_name)
}
