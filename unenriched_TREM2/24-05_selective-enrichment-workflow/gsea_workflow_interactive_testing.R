# Perform interactive AUcell ranking and gene-set enrichment 
# ...to investigate impact of log transformation on AUCell

##  ............................................................................
##  Load packages                                                           ####
library(AUCell)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(purrr)
library(argparse)
library(cli)
library(matrixStats)
library(doParallel)
library(variancePartition)
library(waldo)

sce_path     <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/sce_objects/initial_cohort/split_sce_celltype/celltype_sce/Micro_sce.qs"
# geneset_path <- "/Volumes/jmm17/home/ukdri_projects_data/unenriched_trem2/gazestani_et_al.qs"
geneset_path <- "/Volumes/jmm17/projects/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/assets/manuscript_geneset/Mancuso_2024.qs"
outdir       <- "/Volumes/jmm17/home/archive"
# geneset_name <- "gazestani"
geneset_name <- "mancuso"


# set working directory
setwd(outdir)
outdir <- "selective_enrichment_test_aucell_exprmatrix"
dir.create(outdir)

##  ............................................................................
##  Read Geneset                                                            ####

genesets_l <- qs::qread(geneset_path)

# alternative genesets
#database_dir <- "/rds/general/project/ukdrmultiomicsproject/live/Users/Nurun/gsea_databases"
#database_list <- list.files(database_dir, pattern = "GO_Biological_Process_2021.gmt", full.names = TRUE)
#genesets_l <- fgsea::gmtPathways(database_list)
#geneset_name <- "all_go_bp"

##  ............................................................................
##  Read and Process SCE object                                             ####

sce <- qs::qread(sce_path)

sce <- sce[!is.na(SummarizedExperiment::rowData(sce)$gene), ]
sce <- sce[!duplicated(SummarizedExperiment::rowData(sce)$gene), ]
rownames(sce) <- SummarizedExperiment::rowData(sce)$gene

celltype <- "micro"

##  ............................................................................
##  split celltype sce by TREM2Variant, and perform analysis                ####

variant <- c("CV", "R47H", "R62H")
  
sce_subset <- sce[ , sce[["TREM2Variant"]] %in% variant[[2]]]
sce_subset$manifest <- droplevels(sce_subset$manifest)

sce_subset <- sce_subset[!is.na(SummarizedExperiment::rowData(sce_subset)$gene), ]
sce_subset <- sce_subset[!duplicated(SummarizedExperiment::rowData(sce_subset)$gene), ]
rownames(sce_subset) <- SummarizedExperiment::rowData(sce_subset)$gene

#   expr_mat <- scater::calculateCPM(sce_subset, exprs_values = "counts")
#   expr_mat <- as.matrix(expr_mat, "dgCMatrix")
#   expr_mat <- as.matrix(log2(expr_mat + 1))

expr_mat_rawcounts <- SingleCellExperiment::counts(sce_subset)
expr_mat_rawcounts <- as.matrix(expr_mat_rawcounts)

expr_mat_cpm <- scater::calculateCPM(sce_subset, exprs_values = "counts")
expr_mat_cpm <- as.matrix(expr_mat_cpm)
# expr_mat_cpm <- as.matrix(log2(expr_mat_cpm + 1))

expr_mat_log2 <- scater::calculateCPM(sce_subset, exprs_values = "counts")
#   expr_mat_log2 <- as.matrix(expr_mat_log2, "dgCMatrix")
expr_mat_log2 <- as.matrix(log2(expr_mat_log2 + 1))

##  ............................................................................
##  get celltype name and process expression matrix                         ####
  # 
cells_rankings_rawcounts <- AUCell_buildRankings(exprMat = expr_mat_rawcounts,
                                       nCores= future::availableCores(),
                                       plotStats=FALSE)
                                       
cells_rankings_cpm <- AUCell_buildRankings(exprMat = expr_mat_cpm,
                                       nCores= future::availableCores(),
                                       plotStats=FALSE)

cells_rankings_log2 <- AUCell_buildRankings(exprMat = expr_mat_log2,
                                       nCores= future::availableCores(),
                                       plotStats=FALSE)
# 

#   aucell_cellrankings_dir <- sprintf("%s/%s_%s/%s", outdir, "aucell_cell_rankings", geneset_name, celltype)
#   dir.create(aucell_cellrankings_dir, recursive = TRUE, showWarnings = FALSE)
#   qs::qsave(cells_rankings, file = sprintf("%s/%s_cells_rankings_%s.qs", aucell_cellrankings_dir, celltype, variant))

# alt workflow: reading in pre-generated rankings
# cli::cli_text("Loading cells rankings")
# cells_rankings <- qs::qread(sprintf("aucell/%s_cells_rankings.qs", celltype))
  
set.seed(123)

cells_AUC_rawcounts <- AUCell_calcAUC(geneSets = genesets_l, 
                            rankings = cells_rankings_rawcounts,
                            nCores = future::availableCores(), 
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings_rawcounts)))

cells_AUC_cpm <- AUCell_calcAUC(geneSets = genesets_l, 
                                      rankings = cells_rankings_cpm,
                                      nCores = future::availableCores(), 
                                      aucMaxRank = ceiling(0.05 * nrow(cells_rankings_cpm)))

cells_AUC_log2 <- AUCell_calcAUC(geneSets = genesets_l, 
                                      rankings = cells_rankings_log2,
                                      nCores = future::availableCores(), 
                                      aucMaxRank = ceiling(0.05 * nrow(cells_rankings_log2)))
rm(cells_rankings_cpm)
rm(cells_rankings_rawcounts)
rm(cells_rankings_log2)
gc()

##  ............................................................................
##  dream: process, run, write                                              ####

AUscore_raw <- as.matrix(getAUC(cells_AUC_rawcounts))
AUscore_cpm <- as.matrix(getAUC(cells_AUC_cpm))
AUscore_log2 <- as.matrix(getAUC(cells_AUC_log2))

metadata <- colData(sce_subset) %>% as.data.frame()
metadata <- metadata %>% 
  mutate(total_features_by_counts = as.vector(scale(total_features_by_counts)),
         pc_mito = as.vector(scale(pc_mito)),
         Age = as.vector(scale(Age)),
         PostMortemInterval = as.vector(scale(PostMortemInterval)),
         NeuropathologicalDiagnosis = factor(NeuropathologicalDiagnosis, levels = c("Control", "AD")))

rm(sce_subset)
gc()

param = SnowParam(future::availableCores(), "SOCK", progressbar=TRUE)
register(param)

##  ............................................................................
##  pct4G8PositiveArea                                                      ####

form <- ~ pct4G8PositiveArea + (1| manifest) + total_features_by_counts + pc_mito + Sex + Age + PostMortemInterval + BrainRegion + CD33Group + APOEgroup

fitmm_raw = dream( AUscore_raw, form, metadata)
fitmm_cpm = dream( AUscore_cpm, form, metadata)
fitmm_log2 = dream( AUscore_log2, form, metadata)


# note: coef is the column name / number specifying interest group: use column 2, because for categorical variables, an error occurs due to
# all groups being added to the name of the variable
res_dream_raw <- topTable( fitmm_raw, coef=2, number=Inf )
res_dream_cpm <- topTable( fitmm_cpm, coef=2, number=Inf )
res_dream_log2 <- topTable( fitmm_log2, coef=2, number=Inf )


waldo::compare(res_dream_raw, res_dream_cpm, x_arg = "raw", y_arg = "cpm")
waldo::compare(res_dream_raw, res_dream_log2, x_arg = "raw", y_arg = "log2")
waldo::compare(res_dream_cpm, res_dream_log2, x_arg = "cpm", y_arg = "log2")