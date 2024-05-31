#!/usr/bin/env Rscript
# Perform differential expression analysis 

# helper function ####
filter_sce_genes_for_de <- function(sce,
                                    min_counts = 1L,
                                    min_cells_pc = 0.1) {
  n_genes_before <- dim(sce)[[1]]
  
  min_cells <- floor(dim(sce)[[2]] * min_cells_pc)
  
  keep_genes <- apply(
    SingleCellExperiment::counts(sce), 1,
    function(x) {
      length(x[x >= min_counts]) > min_cells
    }
  )
  
  sce <- sce[keep_genes, ]
  
  n_genes_after <- dim(sce)[[1]]
  
  cli::cli_text(c(
    "Selected {n_genes_after} from {n_genes_before} genes ",
    "with >{min_counts} count(s) in >{min_cells_pc * 100}% of cells."
  ))
  
  return(sce)
}

#   ____________________________________________________________________________

#   Initialization                                                          ####

##  ............................................................................
##  Load packages                                                           ####
library(apcluster)
library(scFlow)
library(tidyverse)
library(argparse)

##  ............................................................................
##  Parse command-line arguments                                            ####

# create parser object
parser <- ArgumentParser()

# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")

required$add_argument(
  "--sce",
  help = "path to sce.qs object",
  metavar = "/dir/sce/sce.qs",
  required = TRUE
)

required$add_argument(
  "--dependent_var",
  help = "dependent variable",
  metavar = "group",
  required = TRUE
)

required$add_argument(
  "--ref_class",
  help = "reference class within dependent variable",
  metavar = "NULL",
  required = TRUE
)

required$add_argument(
  "--confounding_vars",
  help = "confounding variables",
  metavar = "age,sex,pc_mito",
  required = TRUE
)

required$add_argument(
  "--stratification_var",
  help = "stratified by var",
  metavar = "TREM2",
  required = TRUE
)

required$add_argument(
  "--subset_var",
  help = "subset group for variance calculation",
  metavar = "Diagnosis",
  required = TRUE
)

required$add_argument(
  "--subset_class",
  help = "subset class for variance calculations",
  metavar = "Control",
  required = TRUE
)

required$add_argument(
  "--output_dir",
  help = "output dir path",
  metavar = "de_results",
  required = TRUE
)

required$add_argument(
  "--ensembl_mappings",
  help = "file path to enesembl mapping",
  metavar = "ensembl_mappings_human.tsv",
  required = TRUE
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
args <- parser$parse_args()

# set working directory

celltype <- gsub("_sce.qs", "", basename(args$sce))

args$confounding_vars <- strsplit(args$confounding_vars, ",")[[1]]
confounding_vars_default <- c("cngeneson", "pc_mito")
mod <- paste(c(confounding_vars_default, args$confounding_vars), collapse = "_")

if (args$ref_class == "NULL") {
  args$ref_class <- NULL
}

outdir <- sprintf("%s/de_%s_%s_%s/%s", args$output_dir, args$stratification_var, args$dependent_var, mod, celltype)
cli::cli_text(c(
  c("Creating ", cli::col_green(c(outdir, " \n")))
))
dir.create(outdir, recursive = TRUE)


#### preparing sce ####

cli::cli_text(c(
  c("Reading {.strong {celltype}} from ", cli::col_green(c(args$sce, " \n")))
))

sce <- qs::qread(args$sce)
rownames(sce) <- SingleCellExperiment::rowData(sce)$ensembl_gene_id


cli::cli_text("{.strong manifest} with less than {.strong 10} nuclei are excluded!")

idx <- names(table(sce$manifest)[(table(sce$manifest) >= 10)])
sce <- sce[ , sce$manifest %in% idx]
sce$manifest <- droplevels(sce$manifest)

# filter sce to have a common set of gene for all TREM2Variant
sce <- filter_sce_genes_for_de(sce, min_counts = 1, min_cells_pc = 0.1)


for(i in unique(sce[[args$stratification_var]])){
  cli::cli_h2("Subsetting {.strong {i}} from {.emph {args$stratification_var}}")
  
  sce_subset <- sce[ , sce[[args$stratification_var]] %in% i]

  exclude_vars <- c()

  for(vars in c(args$dependent_var, args$confounding_vars)){
  sce_subset[[vars]] <- as.character(sce_subset[[vars]])
  sce_subset[[vars]] <- ifelse(sce_subset[[vars]] == "N/A", NA, sce_subset[[vars]])
  sce_subset <- sce_subset[ , !is.na(sce_subset[[vars]])]
  is_numeric <- !any(is.na(as.numeric(sce_subset[[vars]])))
  if (is_numeric) {
    print(c(vars, "contains numeric values, converting it to numeric"))
    sce_subset[[vars]] <- as.numeric(sce_subset[[vars]])
  } else {
    print(c(vars, "does not contain numeric values, converting it to factor"))
    sce_subset[[vars]] <- as.factor(sce_subset[[vars]])
    n_lvl <- length(levels(sce_subset[[vars]]))
    if(n_lvl < 2) {exclude_vars <- c(exclude_vars, vars)}
  }
  }

  if(!is.null(exclude_vars)){
  cli::cli_text("Dropping {.emph {exclude_vars}} from confounding_vars!")
  confounding_vars <- setdiff(args$confounding_vars, exclude_vars)
  } else {
    confounding_vars <- args$confounding_vars
  }

  sce_subset$manifest <- droplevels(sce_subset$manifest)

  if(!is.numeric(args$dependent_var)){
    assertthat::assert_that(
      length(levels(sce_subset[[args$dependent_var]])) > 1,
      msg = "The length of the depending variable is less than 2."
    )
  }


  de_results <- perform_de(
    sce = sce_subset,
    de_method = "MASTZLM",
    ebayes = FALSE,
    mast_method = "glmer",
    min_counts = 0,
    min_cells_pc = 0.10,
    rescale_numerics = TRUE,
    dependent_var = args$dependent_var,
    ref_class = args$ref_class,
    confounding_vars = c(confounding_vars_default, confounding_vars),
    random_effects_var = "manifest",
    ensembl_mapping_file = args$ensembl_mappings,
    unique_id_var = "manifest",
    subset_var = args$subset_var,
    subset_class = args$subset_class,
    nAGQ = 0,
    parallel = TRUE
  )


file_name <- sprintf("%s/%s_%s_in_%s.tsv", outdir, celltype, args$dependent_var, i)
write.table(
    x = de_results[[1]], 
    file = file_name,
    col.names = TRUE, row.names = FALSE, sep = "\t")

}
