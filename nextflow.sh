#!/bin/bash
#
#SBATCH --time=20:00:00
#SBATCH --job-name=LPS_pipeline_ONT
#SBATCH --output=./s%j_job.pipeline_ONT.out
#SBATCH --error=./s%j_job.pipeline_ONT.error
#SBATCH --account=a_uqds
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

module load nextflow/23.04.2

# Directory containing the nextflow.config file and the main.nf script
dir=/my/project/directory/LPS/ONT
cd ${dir}

# Samplesheet file
samplesheet=${dir}/samplesheet/samples_test.csv

# Directory that will be created to contain the output files
out_dir=${dir}/results

# Bunya Slurm account 
slurm_account='a_uqds'

##ii) Typing workflow
#directory containing the Nanopore basecalled fastq files
fqdir=${dir}/fastq
nextflow main.nf --outdir ${out_dir} --fqdir ${fqdir} --samplesheet ${samplesheet} -resume --slurm_account ${slurm_account}

# Add -resume to continue a run if interrupted
# include --skip_download_... if databases have already been downloaded
nextflow run main.nf -profile singularity -resume --outdir ${out_dir} --fqdir ${fqdir} --samplesheet ${samplesheet} -c nextflow_custom.config \
--skip_download_checkm_db true --skip_download_sylph_db true --skip_download_bakta_db true --skip_download_amrfinder_db true \
--slurm_account ${slurm_account}

# If you just want a (sub)type assignment, you can skip most of the annotation steps.
# You will need an assembly, but can skip polishing (though this may influence the subtyping)
# nextflow run main.nf -profile singularity -resume --outdir ${out_dir} --fqdir ${fqdir} --samplesheet ${samplesheet} -c nextflow_custom.config \
# --skip_download_checkm_db true --skip_download_sylph_db true --skip_download_bakta_db true --skip_download_amrfinder_db true \
# --skip_polishing true --skip_mlst true --skip_quast true --skip_checkm true --skip_bakta true --skip_amrfinder true \
# --slurm_account ${slurm_account}