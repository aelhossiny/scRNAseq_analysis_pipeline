#!/bin/bash

#SBATCH --account=eicarpen99
#SBATCH --job-name='cellranger_hg38_%a'
#SBATCH --output=logs/cellranger_hg38_%a.log

#SBATCH --partition=standard
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mail-user=eicarpen@med.umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL   # When to send emails: BEGIN, END, FAIL, or ALL

#SBATCH --array=1

mkdir -p logs

# Import modules (Choose the version of cellranger that you would like to use)
## You check the available versions of cellranger by running the following command: module avail cellranger
ml cellranger/9.0.0

# Setting up variables
samples_manifest=samples_manifest.txt
ref=/nfs/turbo/umms-eicarpen/references/refdata-gex-GRCh38-2020-A
outputDir=/nfs/turbo/umms-eicarpen/longitudinal_PDAC_project/results/00_alignment/Realigned_CellRanger9
# mkdir -p $outputDir

# Extracting sample for the current job 
current_sample_name=$(cat $samples_manifest | cut -f1 | sed -n ${SLURM_ARRAY_TASK_ID}p)
current_sample_path=$(cat $samples_manifest | cut -f2 | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "Analyzing Sample $current_sample_name" 
echo -e "Fastq path $current_sample_path"

# Running cellranger
cellranger count --id=${current_sample_name} \
--create-bam true \
--transcriptome=$ref \
--fastqs=$current_sample_path \
--sample=$current_sample_name \
--localcores=16

mv ${current_sample_name} $outputDir

