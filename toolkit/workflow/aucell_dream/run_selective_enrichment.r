#!/usr/bin/env Rscript

# Perform AUcell ranking and gene-set enrichment

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

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

# specify options, add and assign variables
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")

required$add_argument(
  "--sce_path",
  help = "path to sce object",
  metavar = "sce object",
  required = TRUE
)

required$add_argument(
  "--geneset_path",
  help = "path to geneset object",
  metavar = "geneset object",
  required = TRUE
)

required$add_argument(
  "--outdir",
  help = "path to output directory",
  metavar = "output directory ",
  required = TRUE
)

required$add_argument(
  "--geneset_name",
  help = "name of geneset analysed",
  metavar = "geneset name",
  required = TRUE
)

args <- parser$parse_args()

sce_path     <- args$sce_path
geneset_path <- args$geneset_path
outdir       <- args$outdir
geneset_name <- args$geneset_name

# set working directory
setwd(outdir)
outdir <- "selective_enrichment"
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

cli::cli_text("Loading  {args$sce_path} ")
sce <- qs::qread(sce_path)

sce <- sce[!is.na(SummarizedExperiment::rowData(sce)$gene), ]
sce <- sce[!duplicated(SummarizedExperiment::rowData(sce)$gene), ]
rownames(sce) <- SummarizedExperiment::rowData(sce)$gene


##  ............................................................................
##  get celltype name and process expression matrix                         ####

celltype <- gsub("_sce", "", tools::file_path_sans_ext(basename(args$sce)))

cli::cli_text("Normalising expression matrix")
 
expr_mat <- scater::calculateCPM(sce, exprs_values = "counts")
expr_mat <- as.matrix(log2(expr_mat + 1))
 
##  ............................................................................
##  get celltype name and process expression matrix                         ####

cli::cli_text("Building cells rankings")
# 
cells_rankings <- AUCell_buildRankings(exprMat = expr_mat,
                                       nCores= future::availableCores(),
                                       plotStats=FALSE)
# 
dir.create("aucell")
qs::qsave(cells_rankings, file = sprintf("aucell/%s_cells_rankings.qs", celltype))

# alt workflow: reading in pre-generated rankings
  # cli::cli_text("Loading cells rankings")
  # cells_rankings <- qs::qread(sprintf("aucell/%s_cells_rankings.qs", celltype))

cli::cli_text("Calculating AUC")

set.seed(123)
cells_AUC <- AUCell_calcAUC(geneSets = genesets_l, 
                            rankings = cells_rankings,
                            nCores = future::availableCores(), 
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

rm(cells_rankings)
gc()

##  ............................................................................
##  dream: process, run, write                                              ####

AUscore <- as.matrix(getAUC(cells_AUC))

# dir.create("AUscore")
# qs::qsave(AUscore, file = sprintf("AUscore/%s_auscore.qs", celltype))

min_val <- 1
AUscore_log <- log2(AUscore+min_val)

rm(cells_AUC)
gc()

metadata <- colData(sce) %>% as.data.frame()
metadata <- metadata %>% 
  mutate(total_features_by_counts = as.vector(scale(total_features_by_counts)),
         pc_mito = as.vector(scale(pc_mito)),
         Age = as.vector(scale(Age)),
         PostMortemInterval = as.vector(scale(PostMortemInterval)),
         NeuropathologicalDiagnosis = factor(NeuropathologicalDiagnosis, levels = c("Control", "AD")))

rm(sce)
gc()

outdir_name <- sprintf("%s/%s_%s", outdir, "enrichment_res", geneset_name)
dir.create(outdir_name)

param = SnowParam(future::availableCores(), "SOCK", progressbar=TRUE)
register(param)

##  ............................................................................
##  pct4G8PositiveArea                                                      ####

form <- ~ pct4G8PositiveArea + (1| manifest) + total_features_by_counts + pc_mito + BrainRegion + Sex + Age + PostMortemInterval

cli::cli_text(paste0("Calculating dream: ", form)

fitmm = dream( AUscore_log, form, metadata)

# note: coef is the column name / number specifying interest group: use column 2, because for categorical variables, an error occurs due to
# all groups being added to the name of the variable
res_dream <- topTable( fitmm, coef=2, number=Inf )

res_dream$geneset <- rownames(res_dream)
res_dream$form <- paste0(as.character(form), collapse = "")
res_dream$celltype <- celltype

write.table(x = res_dream, 
            file = sprintf("%s/%s_%s", outdir_name, celltype, "res_dream_pct4G8PositiveArea.tsv"),
            sep = "\t", 
            row.names = F) # row.names false to avoid issue with row and column name conflicts for excel users


##  ............................................................................
##  pctPHF1PositiveArea                                                     ####

form <- ~ pctPHF1PositiveArea + (1| manifest) + total_features_by_counts + pc_mito + BrainRegion + Sex + Age + PostMortemInterval
cli::cli_text(paste0("Calculating dream: ", form)

fitmm = dream( AUscore_log, form, metadata)
res_dream <- topTable( fitmm, coef=2, number=Inf )

res_dream$geneset <- rownames(res_dream)
res_dream$form <- paste0(as.character(form), collapse = "")
res_dream$celltype <- celltype

write.table(x = res_dream, 
            file = sprintf("%s/%s_%s", outdir_name, celltype, "res_dream_pctPHF1PositiveArea.tsv"),
            sep = "\t", 
            row.names = F)


##  ............................................................................
##  NPD                                                                     ####

form <- ~ NeuropathologicalDiagnosis + (1| manifest) + total_features_by_counts + pc_mito + BrainRegion + Sex + Age + PostMortemInterval
cli::cli_text(paste0("Calculating dream: ", form)

fitmm = dream( AUscore_log, form, metadata)
res_dream <- topTable( fitmm, coef=2, number=Inf )

res_dream$geneset <- rownames(res_dream)
res_dream$form <- paste0(as.character(form), collapse = "")
res_dream$celltype <- celltype

write.table(x = res_dream, 
            file = sprintf("%s/%s_%s", outdir_name, celltype, "res_dream_NeuropathologicalDiagnosis.tsv"),
            sep = "\t", 
            row.names = F)