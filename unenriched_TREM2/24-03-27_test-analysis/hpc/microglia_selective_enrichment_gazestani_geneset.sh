#!/bin/bash

#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=10:mem=32gb
#PBS -N gazestani_selecive_enrichment_1
#PBS -o /rds/general/user/jmm17/home/ukdri-projects_code/unenriched_TREM2/24-03-27_test-analysis/logs/gazestani_selecive_enrichment_conda.out
#PBS -e /rds/general/user/jmm17/home/ukdri-projects_code/unenriched_TREM2/24-03-27_test-analysis/logs/gazestani_selecive_enrichment_conda.err

cd $PBS_O_WORKDIR

module load anaconda3/personal

sce_path="/rds/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/results/final/sce_rds_obj/sce.rds"
geneset_path="~/ukdri_projects_data/unenriched_trem2/gazestani_et_al.qs"

script_path="~/ukdri-projects_code/unenriched_TREM2/24-03-27_test-analysis/hpc_microglia_selective_enrichment_gazestani_geneset.R"

START=$(date)

echo job started at $START

Rscript ~/ukdri-projects_code/unenriched_TREM2/24-03-27_test-analysis/scripts/microglia_selective_enrichment_gazestani_geneset.R \
--sce_load_path=$sce_load_path \
--geneset_path=$geneset_path

END=$(date)
echo job ended at $END