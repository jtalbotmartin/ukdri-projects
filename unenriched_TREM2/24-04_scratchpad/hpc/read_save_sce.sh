#!/bin/bash

#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=10:mem=32gb
#PBS -N rds_sce
#PBS -o rds_sce.out
#PBS -e rds_sce.err

cd $PBS_O_WORKDIR

sce_load_path="/data/MAP_analysis/TREM2_unenriched_scflow/results/final/SCE/final_sce"
sce_save_path="/data/MAP_analysis/TREM2_unenriched_scflow/results/final/sce_rds_obj/sce.rds"

singularity_container_path="/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines/snRNAseq/singularity-cache/nfancy-scflow-0.7.4.img"

script_path="/rds/general/user/jmm17/home/ukdri-projects_code/unenriched_TREM2/24-03-27_test-analysis/scripts/read_save_sce.R"

START=$(date)

echo job started at $START

singularity exec \
-B /rds/general/ephemeral/user/$USER/ephemeral:/tmp,/rds/general/ephemeral/user/$USER/ephemeral:/var/tmp,/rds/general/user/jmm17/home:/project,/rds/general/project/ukdrmultiomicsproject/live:/data \
$singularity_container_path \
Rscript $script_path \
--sce_load_path=$sce_load_path \
--sce_save_path=$sce_save_path

# assign other parameters here

END=$(date)
echo job ended at $END