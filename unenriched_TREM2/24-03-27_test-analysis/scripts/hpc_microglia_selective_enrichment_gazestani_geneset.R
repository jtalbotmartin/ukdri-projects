#!/usr/bin/env Rscript
# Perform AUcell ranking and gene-set enrichment of CS genes

#   ____________________________________________________________________________
#   Initialization                                                          ####

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

# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")

required$add_argument(
  "--sce_path",
  help = "path to sce.qs object",
  metavar = "/dir/sce/sce.qs",
  required = TRUE
)

required$add_argument(
  "--geneset_path",
  help = "path to sce.qs object",
  metavar = "/dir/sce/sce.qs",
  required = TRUE
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
args <- parser$parse_args()

sce_path <- args$sce_path
geneset_path <- args$geneset_path

# set working directory
setwd("~/overclockers_backup/projects/ad_unsorted/52_samples/analysis")

outdir <- "selective_enrichment"
dir.create(outdir)


# genesets_l <- qs::qread("selective_enrichment/selected_genesets.qs")
# geneset_name <- "selected"

genesets_l <- qs::qread(geneset_path)
geneset_name <- "gazestani"



#database_dir <- "/rds/general/project/ukdrmultiomicsproject/live/Users/Nurun/gsea_databases"
#database_list <- list.files(database_dir, pattern = "GO_Biological_Process_2021.gmt", full.names = TRUE)
#genesets_l <- fgsea::gmtPathways(database_list)
#geneset_name <- "all_go_bp"

cli::cli_text("Loading  {args$sce} ")

# sce <- qs::qread(args$sce)
sce <- readRDS(sce_path)



sce <- sce[!is.na(SummarizedExperiment::rowData(sce)$gene), ]
sce <- sce[!duplicated(SummarizedExperiment::rowData(sce)$gene), ]
rownames(sce) <- SummarizedExperiment::rowData(sce)$gene

celltype <- gsub("_sce", "", tools::file_path_sans_ext(basename(args$sce)))

cli::cli_text("Normalising expression matrix")
 
expr_mat <- scater::calculateCPM(sce, exprs_values = "counts")
expr_mat <- as.matrix(log2(expr_mat + 1))
 
cli::cli_text("Building cells rankings")
# 
cells_rankings <- AUCell_buildRankings(exprMat = expr_mat,
                                       nCores= future::availableCores(),
                                       plotStats=FALSE)
# 
# 
# qs::qsave(cells_rankings, file = sprintf("aucell/%s_cells_rankings.qs", celltype))

#cli::cli_text("Loading cells rankings")
#cells_rankings <- qs::qread(sprintf("aucell/%s_cells_rankings.qs", celltype))


cli::cli_text("Calculating AUC")
set.seed(123)
cells_AUC <- AUCell_calcAUC(geneSets = genesets_l, 
                            rankings = cells_rankings,
                            nCores = future::availableCores(), 
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

rm(cells_rankings)
gc()

####dream####
AUscore <- as.matrix(getAUC(cells_AUC))

#file_name <- sprintf("%s/AUscore/%s_auscore.qs", outdir, celltype)
#qs::qsave(AUscore, file = file_name)

min_val <- 1
AUscore_log <- log2(AUscore+min_val)

rm(cells_AUC)
gc()

metadata <- colData(sce) %>% as.data.frame()
metadata <- metadata %>% 
  mutate(total_features_by_counts = as.vector(scale(total_features_by_counts)),
         pc_mito = as.vector(scale(pc_mito)),
         age = as.vector(scale(age)),
         PMI = as.vector(scale(PMI)),
         diagnosis = factor(diagnosis, levels = c("Control", "AD")))

rm(sce)
gc()


##dream#

outdir_name <- sprintf("%s/%s_%s", outdir, "enrichment_res", geneset_name)
dir.create(outdir_name)

param = SnowParam(future::availableCores(), "SOCK", progressbar=TRUE)
register(param)

form <- ~ aBeta + (1| manifest) + total_features_by_counts + pc_mito + brain_region + sex + age + PMI
fitmm = dream( AUscore_log, form, metadata)
res_dream <- topTable( fitmm, coef='aBeta', number=Inf )

write.table(x = res_dream, 
            file = sprintf("%s/%s_%s", outdir_name, celltype, "res_dream_aBeta.tsv"),
            sep = "\t", 
            row.names = T)


form <- ~ pTau + (1| manifest) + total_features_by_counts + pc_mito + brain_region + sex + age + PMI
fitmm = dream( AUscore_log, form, metadata)
res_dream <- topTable( fitmm, coef='pTau', number=Inf )

write.table(x = res_dream, 
            file = sprintf("%s/%s_%s", outdir_name, celltype, "res_dream_pTau.tsv"),
            sep = "\t", 
            row.names = T)


form <- ~ diagnosis + (1| manifest) + total_features_by_counts + pc_mito + brain_region + sex + age + PMI
fitmm = dream( AUscore_log, form, metadata)
res_dream <- topTable( fitmm, coef='diagnosisAD', number=Inf )

write.table(x = res_dream, 
            file = sprintf("%s/%s_%s", outdir_name, celltype, "res_dream_diagnosis.tsv"),
            sep = "\t", 
            row.names = T)


