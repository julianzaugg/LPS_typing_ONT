# LPS_typing_ONT
Bioinformatics pipeline for Pasteurella multocida LPS typing using ONT sequencing data

  - [Overall pipeline](#Overall-pipeline)
  - [User guide](#Step-by-step-user-guide)
  - [Example data](#Example-data)
  - [Optional parameters](#Optional-parameters)
  - [Output files](#structure-of-the-output-folders)
  - [Assembly mode](#Running-the-workflow-in-assembly-mode-for-other-organisms)
  - [Troubleshooting](#Troubleshooting)
  - [Acknowledgements/citations/credits](#acknowledgements--citations--credits)
    
## Overall pipeline 

### 1. Basecalling

The raw Nanopore in pod5 format are basecalled using [Dorado](https://github.com/nanoporetech/dorado). The barcoding kit should be specified to enable read demultiplexing (e.g "SQK-NBD114-24" or "SQK-NBD114-96"). By default, the basecalling model selected is the super-accuracy ("sup") model (parameter "basecalling_model"). The folder containing the pod5 files is an input parameter to the pipeline (parameter "pod5_dir"). 

### 2. Nanopore reads quality metrics 

[Nanocomp](https://github.com/wdecoster/nanocomp) is used to compute Nanopore read metrics (e.g. Median Read Length, Read N50, Median Read Quality). Those metrics are computed for each barcode on the raw reads included in the basecalled fastq files.  

### 3. Flye assembly and polishing

- The Nanopore reads are assembled using the software [Flye](https://github.com/fenderglass/Flye). By default, Flye will use a reduced coverage for initial disjointig assembly to speed up the assemby (config parameter flye_args= "--asm-coverage 100").  
- The draft assemblies are subsequently polished using [Medaka](https://github.com/nanoporetech/medaka). The model parameter selected to run Medaka (e.g. r1041_e82_400bps_sup_v5.0.0) must correspond to the model used for the basecalling (e.g. dna_r10.4.1_e8.2_400bps_sup.cfg).  

### 4. 	Assembly quality assessment with QUAST

The software [QUAST](https://quast.sourceforge.net/quast.html) is used to compute genome assembly metrics on the polished assemblies.  

### 5. Assembly quality assessment with CheckM

The software [CheckM](https://github.com/Ecogenomics/CheckM) v1 (command [lineage_wf](https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow)) is used to compute genome assembly completeness and contamination, based on the presence or absence of marker genes. 

### 6. Centrifuge taxonomy classification

Nanopore reads are used as input to the taxonomy classifier [Centrifuge](https://ccb.jhu.edu/software/centrifuge/). The default database was downloaded from https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz. 

### 7. LPS typing using Kaptive

The LPS type of the sample is obtained using the software [Kaptive](https://kaptive.readthedocs.io/en/latest/) v3. The genomes assemblies are used as input to this tool. The 9 LPS database is used but can be modified (config parameter "kaptive_db_9lps").  

### 8. 	Variant calling using Clair3

- The reads are mapped to the reference LPS type sequence identified by Kaptive using the sequence alignment program [Minimap2](https://github.com/lh3/minimap2).   
- Then [clair3](https://github.com/HKU-BAL/Clair3/) is used to call variants in the reads as compared to the reference LPS type sequence. The clair3 model must be selected to match the ONT platform used and the basecalling model used (e.g. "r1041_e82_400bps_sup_v500"). The list of clair3 models is available at https://github.com/nanoporetech/rerio/tree/master/clair3_models.  
- Then [SnpEff](https://pcingola.github.io/SnpEff/) is used to annotate and predict the effects of the variants on genes and proteins (such as amino acid changes).  
- Finally, [SnpSift](https://pcingola.github.io/SnpEff/#snpsift) is used to extract the variants predicted to have a high impact on the protein (frameshift and stop_gained variants).    

### 9. 	MLST typing

The software [mlst](https://github.com/tseemann/mlst) is used to scan the genome assemblies against the  PubMLST typing scheme "pmultocida_2" by default (RIRDC). The typing scheme can be modified by specifying the parameter --mlst_scheme (e.g. --mlst_scheme "pmultocida").   

## Step by step user guide

Some files required to use the pipeline are provided to the user (see sections 1a, 1b and 2 below). Additional files must be created/modified by the user (see sections 1c, 3 and 4 below). 

**1) Clone the Github pipeline repository**

Navigate to a folder to your scratch space where you would like to run the pipeline (e.g. $raw_dir below) and clone the pipeline repository to import the required files: 
```
raw_dir=/scratch/project_mnt/SXXX/PIPELINE
cd $raw_dir
git clone https://github.com/vmurigneu/LPS_typing.git
```

It will create a repository called "LPS_typing". The following three files can be found in the pipeline repository:

- **a) Nextflow configuration file (nextflow.config)**  

When a Nexflow pipeline script is launched, Nextflow looks for a file named **nextflow.config** in the current directory. The configuration file defines default parameters values for the pipeline and cluster settings such as the executor (e.g. "slurm", "local") and queues to be used (https://www.nextflow.io/docs/latest/config.html).  

The pipeline uses separated Singularity containers for all processes. Nextflow will automatically pull the singularity images required to run the pipeline and cache those images in the singularity directory in the pipeline work directory by default or in the singularity.cacheDir specified in the nextflow.config file ([see documentation](https://www.nextflow.io/docs/latest/singularity.html)). Ensure that you have sufficient space in your assigned singularity directory as images can be large.   

An example configuration file can be found [here](https://github.com/vmurigneu/LPS_typing/blob/main/nextflow.config). 

- **b) Nextflow main script (main.nf)**

The main.nf script contains the pipeline code and is generally not user-modifiable. 

- **c) Nextflow execution bash script (nextflow.sh)**

This is the bash script used to launch the workflow on the HPC. The template slurm script provided can be used to launch the pipeline on UQ HPC Bunya and is available [here](https://github.com/vmurigneu/LPS_typing/blob/main/nextflow.sh). This file should be modified by the user to provide the path to the samplesheet file, Nanopore data files etc (see section "Step by step user guide" below). 

**2) Database files for Centrifuge, Kaptive, Minimap2 and CheckM**

Copy the databases folder from the RDM to the cloned pipeline repository on the scratch space (named "dir" below):
```
dir=/scratch/project_mnt/SXXX/PIPELINE/LPS_typing
cp -r /QRISdata/Q2313/Valentine/PIPELINES/databases ${dir}
```

**3) Prepare the samplesheet file (csv)**

The samplesheet file is a comma-separated values files that defines the names of the samples with their corresponding barcode and input fastq files. The header line should match the header line in the examples below. The samplesheet can be saved in a folder named samplesheet e.g. 
```
mkdir /scratch/project_mnt/SXXX/LPS_typing/samplesheet
vim /scratch/project_mnt/SXXX/LPS_typing/samplesheet/samples.csv
```

* **Basecalling and typing workflow** (soon)

The samplesheet contains one line for each sample with the following information: the barcode identifier from the Nanopore barcoding kit (column "barcode_id") and the sample identifier (column "sample_id").
```
barcode_id,sample_id
barcode17,PM1947
barcode18,PM1422
```

* **Typing workflow**
   
If basecalling is performed outside of the pipeline, the raw ONT pod5 files needs to be converted into basecalled reads using Dorado (fastq files). Then the user must copy the basecalled files in a directory (parameter "--fqdir") and specify the path to those files in the samplesheet file.   
The samplesheet contains one line for each sample with the following information: the sample identifier (column "sample_id") and the path to the corresponding basecalled read file (column "long_fastq"). File paths are given in relation to the workflow base directory, they are not absolute paths. 
```
sample_id,long_fastq
PM1947,fastq/barcode17.simplex_duplex.fastq.gz
PM1422,fastq/barcode18.simplex_duplex.fastq.gz
```

**4) Run the pipeline**

The pipeline will be launched on the HPC Bunya using the bash script nextflow.sh.   

* **Basecalling and typing workflow** (soon)

The raw ONT pod5 files must be copied in a directory (parameter "--pod5_dir").
```
pod5=/scratch/project_mnt/SXXX/LPS_typing_pipeline/pod5
mkdir $pod5
cp /path/to/pod5/files/ $pod5
```

Then the command to start the pipeline is:  
`nextflow main.nf --samplesheet /path/to/samples.csv --pod5_dir /path/to/pod5/directory/ --outdir /path/to/outdir/ --slurm_account 'account' `
```
--samplesheet: path to the samplesheet file
--outdir: path to the output directory to be created
--pod5_dir: path to the directory containing the Nanopore pod5 files
--slurm_account: name of the Bunya account (default='a_uqds') 
```

Once the nextflow.sh file is ready, the user can submit the pipeline on Bunya using the command:
```
sbatch nextflow.sh
```

* **Typing workflow**
  
The user must copy the basecalled files in a directory (parameter "--fqdir") and specify the path to those files in the samplesheet file (see above). 
```
fastq=/scratch/project_mnt/SXXX/LPS_typing_pipeline/fastq
mkdir $fastq
cp /path/to/fastq/files $fastq
```

Then the command to start the pipeline is:  
`nextflow main.nf --samplesheet /path/to/samples.csv --fqdir /path/to/fastq/directory/ --outdir /path/to/outdir/ --slurm_account 'account'`
```
--samplesheet: samplesheet file
--outdir: path to the output directory to be created
--fqdir: path to the directory containing the Nanopore basecalled fastq files
--slurm_account: name of the Bunya account (default='a_uqds') 
```

Note: to run the assembly and assembly metrics steps only (skip LPS typing and variant calling):  
`nextflow main.nf --samplesheet /path/to/samples.csv --fqdir /path/to/fastq/directory/ --outdir /path/to/outdir/ --slurm_account 'account' --skip_kaptive3 --skip_clair3`

Once the nextflow.sh file is ready, the user can submit the pipeline on Bunya using the command:
```
sbatch nextflow.sh
```

## Example data

To test the pipeline, we have provided some test data that contain a subset of 20,000 Nanopore reads for two samples. 

File | Description
---|---
[fastq/barcode17.simplex_duplex.test.fastq.gz](https://github.com/vmurigneu/LPS_typing/blob/main/fastq/barcode17.simplex_duplex.test.fastq.gz) | ONT fastq reads for barcode17
[fastq/barcode18.simplex_duplex.test.fastq.gz](https://github.com/vmurigneu/LPS_typing/blob/main/fastq/barcode18.simplex_duplex.test.fastq.gz) | ONT fastq reads for barcode18
[samplesheet/samples_test.csv](https://github.com/vmurigneu/LPS_typing/blob/main/samplesheet/samples_test.csv) | samplesheet for running the typing pipeline
[nextflow.sh](https://github.com/vmurigneu/LPS_typing/blob/main/nextflow.sh) | Nextflow execution bash script

In the nextflow.sh file, you must modify the directory "dir" line 16, the Bunya Slurm account name for your group line 7 and line 26. Then you can run the pipeline using:
`sbatch nextflow.sh`

## Optional parameters

Some parameters can be added to the command line in order to include or skip some steps and modify some parameters:

1. Basecalling:
* `--skip_basecalling`: skip the basecalling step (default=false)
* `--basecalling_model`: basecalling model (default="sup"), see [details](https://github.com/nanoporetech/dorado?tab=readme-ov-file#automatic-model-selection-complex)
* `--barcoding_kit`: name of the barcoding kit for Dorado (default="SQK-NBD114-24")

2. Nanopore reads quality metrics:
* `--skip_nanocomp`: skip the Nanocomp step (default=false)
* `--nanocomp_threads`: number of threads for the Nanocomp step (default=4)

3. Genome assembly and polishing:
* `--skip_assembly`: skip the assembly step (default=false). Note: it is not recommended to skip assembly as many steps in the downstream processing depends on the assembly results.   
* `--flye_args`: Flye optional parameters (default="--asm-coverage 100"), see [available parameters](https://github.com/mikolmogorov/Flye/blob/flye/docs/USAGE.md#-quick-usage)
* `--flye_threads`: number of threads for the assembly (default=4)
* `--genome_size`: estimated genome size (default="2.3M")
* `--skip_polishing`: skip the Medaka polishing step (default=false)
* `--medaka_threads`: number of threads for the polishing (default=8)
* `--medaka_model`: name of the Medaka model (default="r1041_e82_400bps_sup_v5.0.0"), see [details](https://github.com/nanoporetech/medaka#models)

4. Assembly quality assessment with QUAST:
* `--skip_quast`: skip the QUAST step (default=false)
* `--quast_threads`: number of threads for QUAST (default=4)

5. Assembly quality assessment with CheckM:
* `--skip_checkm`: skip the CheckM step (default=false)
* `--checkm_db`: path to the CheckM database folder (default="../../../databases/CheckM-1.2.2")

6. Centrifuge taxonomy classification:
* `--skip_centrifuge`: skip the Centrifuge classification step (default=false)
* `--skip_download_centrifuge_db`: skip the Centrifuge database downloading step (default=true)
* `--centrifuge_db_download_file`: link to the centrifuge database file to be downloaded (default= 'https://genome-idx.s3.amazonaws.com/centrifuge/nt_2018_3_3.tar.gz')
* `--centrifuge_db`: (default= '../../../databases/centrifuge/nt.*.cff')
* `--centrifuge_threads`: number of threads for Centrifuge classification step (default=6)

7. LPS typing using Kaptive:
* `--skip_kaptive3`: skip the Kaptive typing step (default=false). note: it will automatically skip the variant calling step.  
* `--kaptive_db_9lps`: path to the Kaptive database file (default="../../../databases/v1_kaptive3/9lps.gbk")

8. Variant calling using Clair3:
* `--skip_clair3`: skip the variant calling step (default=false)
* `--minimap_threads`: number of threads for the Minimap2 mapping step (default=6)
* `--clair3_threads`: number of threads for the Clair3 variant calling step (default=4)
* `--clair3_model`: path to the clair3 model folder (default="../../../databases/clair3_models/r1041_e82_400bps_sup_v500")
* `--clair3_args`: Clair3 optional parameters (default="--haploid_sensitive"), see [available parameters](https://github.com/HKU-BAL/Clair3?tab=readme-ov-file#options)
* `--skip_snpeff`: skip the variant annotation step (default=false)
* `--reference_LPS`: path to the file summarising the reference LPS sequence files (default="../../../databases/reference_LPS.txt")

9. MLST typing:
* `--skip_mlst`: skip the MLST typing step (default=false)
* `--mlst_scheme`: MLST typing scheme (default="pmultocida_2")

## Structure of the output folders

The pipeline will create several folders corresponding to the different steps of the pipeline. 
* **2_nanocomp:** Quality control of the Nanopore reads  
  * Full Nanocomp report in html format with plots and metrics table (NanoComp-report.html)    
  * Nanocomp summary text file (NanoStats.txt)
    
The main output folder (`--outdir`) will contain a folder per sample (the folder is named as in the column sample_id in the samplesheet file).

Each sample folder will contain the following folders:
* **3_assembly:** Flye assembly output files (.fasta, .gfa, .gv, .info.txt), see [details](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-flye-output). The final polished asssembly fasta file is sample_id_flye_polished.fasta.
* **4_quast:** QUAST output report file (sample_id_report.tsv).
* **5_checkm:** CheckM output file (sample_id_checkm_lineage_wf_results.tsv).  
* **6_centrifuge:**  Centrifuge taxonomy classification results for the Nanopore reads, see [details](https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-classification-output)
  * Centrifuge classification output: classification assignments for a read (sample_id_centrifuge_species_report.tsv)
  * Centrifuge summary output: classification summary for each genome or taxonomic unit (sample_id_centrifuge_report.tsv)
* **7_kaptive_v3:** Kaptive output files, see [details](https://kaptive.readthedocs.io/en/latest/Outputs.html)
    * LPS type results (sample_id_kaptive_results.tsv)
    * LPS sequence in fasta format (sample_id_flye_polished_kaptive_results.fna)
* **8_clair3:** Mapping files and variant calling results:
    * Minimap2 mapping file in bam format (sample_id_minimap2_mapped.bam and .bai index). Unmapped reads were excluded.
    * Minimap2 mapping statistics (sample_id_minimap2_flagstat.txt)  
    * Clair3 variants (sample_id_clair3.vcf)  
    * Clair3 variants annotated by SnpEff (sample_id_clair3.snpeff.vcf)  
    * Frameshift and stop_gained clair3 variants (sample_id_clair3.snpeff.high_impact.vcf). 
* **9_mlst:** MLST typing output file (sample_id_mlst.csv)
* **10_report:** Summary of results for all samples
    * Flye assembly statistics: assembly coverage, number of contigs, assembly size (3_flye_stats.tsv)  
    * QUAST combined report file (4_quast_report.tsv)  
    * Checkm results (5_checkm_lineage_wf_results.tsv)  
    * Centrifuge taxonomy results:
        - Abundance of P. multocida reads: 6_centrifuge_pasteurella_multocida_species_abundance.tsv
        - Information about the most abundant species identified: 6_centrifuge_most_abundant_species.tsv  
    * Kaptive results (7_kaptive_results.tsv)  
    * Clair3 variants results (8_clair3_snpeff_high_impact.vcf)  
    * MLST results (9_mlst.csv)  
    * Genotype results summarising the variants found in the genotype database (10_genotype_report.tsv)

## Running the workflow in assembly mode for other organisms

The default parameters are suited for Pasteurella multocida. The LPS typing and variant calling are specific to Pasteurella multocida. Here are the paraneters to use the workflow to assemble another species:  
* `--genome_size`: estimated genome size (default="2.3M")
* `--mlst_scheme`: MLST typing scheme (default="pmultocida_2")
* `--skip_kaptive3`: skip the Kaptive typing step (default=false). note: it will automatically skip the variant calling step.

 ## Troubleshooting

Here are a list of things to check if an error occurs:  

* read the error message in the log file
* the samplesheet contains the correct header line (cf https://github.com/vmurigneu/LPS_typing?tab=readme-ov-file#step-by-step-user-guide)
* the fastq files in my samplesheet are present in the fastq folder given in the nextflow.sh script (--fqdir).
* the /databases folder is present in the cloned pipeline repository and not empty (cf https://github.com/vmurigneu/LPS_typing?tab=readme-ov-file#step-by-step-user-guide)
