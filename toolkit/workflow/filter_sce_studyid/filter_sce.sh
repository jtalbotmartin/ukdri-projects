#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:mem=80gb
#PBS -N filter_sce
#PBS -o filter.out
#PBS -e filter.err

cd /rds/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/objects/filter_sce

resource_dir=/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines

mount_dir=/rds

container_resource_dir=/data/general/project/ukdrmultiomicsproject/live/MAP_pipelines

START=$(date)
echo job started at $START

singularity exec \
-B /rds/general/user/$USER/ephemeral/tmp/:/tmp,/rds/general/user/$USER/ephemeral/tmp/:/var/tmp,$mount_dir:/data \
$resource_dir/snRNAseq/singularity-cache/nfancy-scflow-0.7.4.img \
Rscript /data/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/scripts/filter_sce/filter_sce.r
END=$(date)
echo job ended at $END
