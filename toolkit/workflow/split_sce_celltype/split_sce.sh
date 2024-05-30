#!/bin/bash
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:mem=80gb
#PBS -N split_sce
#PBS -o split.out
#PBS -e split.err

cd /rds/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/objects/split_sce_celltype

resource_dir=/rds/general/project/ukdrmultiomicsproject/live/MAP_pipelines

mount_dir=/rds

container_resource_dir=/data/general/project/ukdrmultiomicsproject/live/MAP_pipelines

START=$(date)
echo job started at $START

singularity exec \
-B /rds/general/user/$USER/ephemeral/tmp/:/tmp,/rds/general/user/$USER/ephemeral/tmp/:/var/tmp,$mount_dir:/data \
$resource_dir/snRNAseq/singularity-cache/nfancy-scflow-0.7.4.img \
Rscript /data/general/project/ukdrmultiomicsproject/live/MAP_analysis/TREM2_unenriched_scflow/dge/initial_cohort/scripts/split_sce/split_sce.r
END=$(date)
echo job ended at $END
