#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --job-name=LPS_typing_ONT
#SBATCH --output=./%j_pipeline_ONT.out
#SBATCH --error=./%j_pipeline_ONT.err
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

# Path to samplesheet CSV (FASTQ paths are relative to the Nextflow launch directory)
SAMPLESHEET=/path/to/samplesheet/samples.csv

# Directory where results will be written
OUTDIR=/path/to/results

# Run the full pipeline (databases will be downloaded on first run)
nextflow run main.nf \
  -profile apptainer,slurm \
  --samplesheet ${SAMPLESHEET} \
  --outdir ${OUTDIR} \
  --slurm_account YOUR_ACCOUNT \
  -resume

# If databases have already been downloaded, add the skip flags:
# nextflow run main.nf \
#   -profile apptainer,slurm \
#   --samplesheet ${SAMPLESHEET} \
#   --outdir ${OUTDIR} \
#   --slurm_account YOUR_ACCOUNT \
#   --skip_download_checkm_db true \
#   --skip_download_sylph_db true \
#   --skip_download_bakta_db true \
#   --skip_download_amrfinder_db true \
#   -resume

# Typing-only run (skips most annotation steps):
# nextflow run main.nf \
#   -profile apptainer,slurm \
#   --samplesheet ${SAMPLESHEET} \
#   --outdir ${OUTDIR} \
#   --slurm_account YOUR_ACCOUNT \
#   --skip_download_checkm_db true \
#   --skip_download_sylph_db true \
#   --skip_download_bakta_db true \
#   --skip_download_amrfinder_db true \
#   --skip_polishing true \
#   --skip_mlst true \
#   --skip_quast true \
#   --skip_checkm true \
#   --skip_bakta true \
#   --skip_amrfinder true \
#   -resume
