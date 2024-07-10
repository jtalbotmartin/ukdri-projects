#!/usr/bin/env Rscript

# Getting sce from first TREM2 unenriched run in sce format for interactive use

# load packages                                                           
library(scFlow)
library(argparse)

# parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

# specify options
required <- parser$add_argument_group("Required", "required arguments")

required$add_argument(
  "--sce_load_path",
  help = "directory path for saved sce files",
  metavar = "load_path",
  default = "~"               
)

required$add_argument(
  "--sce_save_path",
  help = "path to save ",
  metavar = "job id",
  default = "~"               
)

args <- parser$parse_args()

sce_load_path <- args$sce_load_path
sce_save_path <- args$sce_save_path

sce <- read_sce(sce_load_path)

# save to .rds object for interactive use
saveRDS(sce, sce_save_path)
