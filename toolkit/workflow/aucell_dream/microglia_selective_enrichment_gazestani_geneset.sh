#!/bin/bash

#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=10:mem=64gb
#PBS -N gazestani_selecive_enrichment
#PBS -o /rds/general/user/jmm17/home/ukdri-projects_code/unenriched_TREM2/24-03-27_test-analysis/logs/gazestani_selecive_enrichment_docker_7.out
#PBS -e /rds/general/user/jmm17/home/ukdri-projects_code/unenriched_TREM2/24-03-27_test-analysis/logs/gazestani_selecive_enrichment_docker_7.err

cd $PBS_O_WORKDIR

singularity_container_path="/rds/general/user/jmm17/home/ukdri_projects_resources/images/snrnaseq-selective-enrichment_latest.img"
# sce_path="/rds/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/results/final/sce_rds_obj/sce.rds"
sce_path="/data/MAP_analysis/TREM2_unenriched_scflow/archive/results_pre_snomics_rerun/final/sce_rds_obj/sce.rds"
# geneset_path="~/ukdri_projects_data/unenriched_trem2/gazestani_et_al.qs"
geneset_path="/project/ukdri_projects_data/unenriched_trem2/gazestani_et_al.qs"
script_path="~/ukdri-projects_code/unenriched_TREM2/24-03-27_test-analysis/hpc/microglia_selective_enrichment_gazestani_geneset.R"

START=$(date)

echo job started at $START
singularity exec \
-B /rds/general/ephemeral/user/$USER/ephemeral:/tmp,/rds/general/ephemeral/user/$USER/ephemeral:/var/tmp,/rds/general/user/jmm17/home:/project,/rds/general/project/ukdrmultiomicsproject/live:/data \
$singularity_container_path \
Rscript ~/ukdri-projects_code/unenriched_TREM2/24-03-27_test-analysis/scripts/microglia_selective_enrichment_gazestani_geneset.R \
--sce_path=$sce_path \
--geneset_path=$geneset_path

END=$(date)
echo job ended at $END