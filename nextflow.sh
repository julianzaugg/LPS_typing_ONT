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

# Collect main results in one file for all samples
cd ${out_dir}
#5_checkm
echo -e  sampleID\\tMarker_lineage\\tNbGenomes\\tNbMarkers\\tNbMarkerSets\\t0\\t1\\t2\\t3\\t4\\t5+\\tCompleteness\\tContamination\\tStrain_heterogeneity > header_checkm
for file in `ls ${out_dir}/*/5_checkm/*checkm_lineage_wf_results.tsv`; do fileName=$(basename $file); sample=${fileName%%_checkm_lineage_wf_results.tsv}; grep -v Bin $file | sed s/^/${sample}_/  >> 5_checkm_lineage_wf_results.tsv.tmp; done
cat header_checkm 5_checkm_lineage_wf_results.tsv.tmp > 5_checkm_lineage_wf_results.tsv

#6_centrifuge
echo -e sampleID\\tname\\ttaxID\\ttaxRank\\tgenomeSize\\tnumReads\\tnumUniqueReads\\tabundance > header_centrifuge
for file in `ls ${out_dir}/*/6_centrifuge/*_centrifuge_report.tsv`; do fileName=$(basename $file); sample=${fileName%%_centrifuge_report.tsv}; grep multocida $file | grep "species"| grep -v subspecies | sed s/^/${sample}\\t/  >> 6_centrifuge_pasteurella_multocida_species_abundance.tsv.tmp; done
cat header_centrifuge 6_centrifuge_pasteurella_multocida_species_abundance.tsv.tmp > 6_centrifuge_pasteurella_multocida_species_abundance.tsv

#7_kaptive
echo -e sampleID\\tBest match locus\\tBest match type\\tMatch confidence\\tProblems\\tIdentity\\tCoverage\\tLength discrepancy\\tExpected genes in locus\\tExpected genes in locus, details\\tMissing expected genes\\tOther genes in locus\\tOther genes in locus, details\\tExpected genes outside locus\\tExpected genes outside locus, details\\tOther genes outside locus\\tOther genes outside locus, details\\tTruncated genes, details\\tExtra genes, details >  header_kaptive3
for file in `ls ${out_dir}/*/7_kaptive_v3/*_kaptive_results.tsv`; do fileName=$(basename $file); sample=${fileName%%_kaptive_results.tsv}; grep -v Assembly $file | sed s/^/${sample}_/  >> 7_kaptive_v3_9lps_results.tsv.tmp; done
cat header_kaptive3 7_kaptive_v3_9lps_results.tsv.tmp > 7_kaptive_v3_9lps_results.tsv

#8_clair3
echo -e  SAMPLEID\\tCHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tSAMPLE > header_clair3
for file in `ls ${out_dir}/*/8_clair3/*_clair3.snpeff.high_impact.vcf`; do fileName=$(basename $file); sample=${fileName%%_clair3.snpeff.high_impact.vcf}; grep -v "^#" $file | sed s/^/${sample}\\t/  >> 8_clair3_snpeff_high_impact.vcf.tmp; done
cat header_clair3 8_clair3_snpeff_high_impact.vcf.tmp > 8_clair3_snpeff_high_impact.vcf

#9_mlst
for file in `ls ${out_dir}/*/9_mlst/*_mlst_pmultocida_rirdc.csv`; do fileName=$(basename $file); sample=${fileName%%_mlst_pmultocida_rirdc.csv};  sed s/^/${sample}_/ $file >> 9_mlst_pmultocida_rirdc.csv; done

#clean folder $out_dir
rm $out_dir/*tmp  $out_dir/header*

