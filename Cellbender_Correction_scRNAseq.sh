#!/bin/bash

#SBATCH --account=marinapa99
#SBATCH --job-name='Cellbender_scRNAseq_%a'
#SBATCH --output=logs/Cellbender_scRNAseq_%a.log
#SBATCH --partition=gpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --array=1-48

## Setting up environment
echo -e ">>> Start time $(date) <<<"
start_time=$(date +%s)

source activate /nfs/turbo/umms-marinapa/Ahmed/conda/envs/cellbender/
outputDir=/nfs/turbo/umms-marinapa/Ahmed/Visium_Project/analysis/outputs/cellbender/scRNAseq
mkdir -p $outputDir

## Importing samples info
samples_info=/nfs/turbo/umms-marinapa/Ahmed/Visium_Project/data/sc_samples_manifest.txt
sample=$(cat $samples_info | cut -f2 | sed -n $[SLURM_ARRAY_TASK_ID]p)
path=$(cat $samples_info | cut -f5 | sed -n $[SLURM_ARRAY_TASK_ID]p)
outputDir=${outputDir}/${sample}
mkdir -p $outputDir
echo -e ">>> Processing $sample <<<"
echo -e ">>> Sample path $path <<<"
echo -e ">>> Output in $outputDir <<<"

## Running CellBender 
cd $outputDir
cellbender remove-background \
--cuda \
--input ${path}/raw_feature_bc_matrix.h5 \
--output ${outputDir}/${sample}.h5

## Reporting time
echo -e ">>> End time $(date) <<<"
end_time=$(date +%s)
runtime_seconds=$((end_time - start_time))
runtime_minutes=$((runtime_seconds / 60))
echo "Total runtime: $runtime_minutes minutes"