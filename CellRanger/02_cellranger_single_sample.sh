#!/bin/bash

#SBATCH --account=eicarpen99
#SBATCH --job-name='cellranger_hg38'
#SBATCH --output=logs/cellranger_hg38.log

#SBATCH --partition=standard
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00

mkdir -p logs

# Import modules (Choose the version of cellranger that you would like to use)
## You check the available versions of cellranger by running the following command: module avail cellranger
ml cellranger/9.0.0

# Running cellranger
# sample id is going to be the name of the fastq file all the way before SXX_RX_XXX.fastq.gz
# for example if the fastq file is 12433-EC-1-GEX_S53_R1_001.fastq.gz, the sample id is going to be 12433-EC-1-GEX
sample_id='12433-EC-1-GEX'
ref='/nfs/turbo/umms-eicarpen/references/refdata-gex-GRCh38-2020-A'
fastqsDir='/nfs/turbo/umms-eicarpen/ScSequencing_Runs/PulmMet/12433-EC/fastqs_12433-EC'

cellranger count --id=$sample_id \
--create-bam true \
--transcriptome=$ref \
--fastqs=$fastqsDir \
--sample=$sample_id \
--localcores=16

