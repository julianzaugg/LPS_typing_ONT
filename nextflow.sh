#!/bin/bash
#
#SBATCH --time=20:00:00
#SBATCH --job-name=LPS_pipeline_ONT
#SBATCH --output=./s%j_job.pipeline_ONT_test.out
#SBATCH --error=./s%j_job.pipeline_ONT_test.error
#SBATCH --account=a_qcif_support
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

module load nextflow/23.04.2

#directory containing the nextflow.config file and the main.nf script
dir=/scratch/project_mnt/S0091/valentine/LPS/PIPELINE_ONT
cd ${dir}

#Samplesheet file
samplesheet=${dir}/samplesheet/samples_test.csv

#Directory that will be created to contain the output files
out_dir=${dir}/results_test

# Bunya Slurm account 
slurm_account='a_qcif_support'

#i) Basecalling and typing workflow
#directory containing the Nanopore raw pod5 files
#pod5_dir=${dir}/pod5
#nextflow main.nf --outdir ${out_dir} --pod5_dir ${pod5_dir} --samplesheet ${samplesheet} -resume --slurm_account ${slurm_account}

##ii) Typing workflow
#directory containing the Nanopore basecalled fastq files
fqdir=${dir}/fastq
nextflow main.nf --outdir ${out_dir} --fqdir ${fqdir} --samplesheet ${samplesheet} -resume --slurm_account ${slurm_account}



