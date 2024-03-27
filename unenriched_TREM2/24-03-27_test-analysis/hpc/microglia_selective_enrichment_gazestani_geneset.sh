#!/bin/bash

#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=10:mem=32gb
#PBS -N rds_sce
#PBS -o rds_sce.out
#PBS -e rds_sce.err

cd $PBS_O_WORKDIR

sce_path="/data/MAP_analysis/TREM2_unenriched_scflow/results/final/sce_rds_obj/sce.rds"
geneset_path="/project/ukdri_projects_data/unenriched_trem2/gazestani_et_al.qs"

singularity_container_path="/rds/general/user/jmm17/home/ukdri_projects_resources/nfancy-scflow-0.7.2.img"

script_path="/rds/general/user/jmm17/home

/ukdri_projects_code/unenriched_TREM2/scripts/24-03-27_test-analysis/read_save_sce.R"
/rds/general/user/jmm17/home
/project/ukdri_projects_data/unenriched_trem2/gazestani_et_al.qs

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