# LPS_typing_ONT
Bioinformatics pipeline for *Pasteurella multocida* LPS typing using Oxford Nanopore sequencing data

## Contents

- [Quick start](#quick-start)
- [Overall pipeline](#overall-pipeline)
- [Database setup](#database-setup)
- [User guide](#step-by-step-user-guide)
- [Example data](#example-data)
- [Optional parameters](#optional-parameters)
- [Output files](#structure-of-the-output-folders)
- [Assembly mode for other organisms](#running-the-workflow-in-assembly-mode-for-other-organisms)
- [Troubleshooting](#troubleshooting)
- [Acknowledgements / citations / credits](#acknowledgements--citations--credits)

---

## Quick start

**Requirements:**
- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [Singularity](https://sylabs.io/singularity/) or [Apptainer](https://apptainer.org/) (the pipeline uses containers for all tools)

**1. Clone the repository**
```bash
git clone https://github.com/vmurigneu/LPS_typing_ONT.git
cd LPS_typing_ONT
```

**2. Obtain databases** — see [Database setup](#database-setup)

**3. Create a samplesheet** — see [User guide](#step-by-step-user-guide)

**4. Run the pipeline**

Local / non-SLURM:
```bash
nextflow run main.nf \
  -profile apptainer \
  --samplesheet /path/to/samples.csv \
  --outdir /path/to/results \
  -resume
```

SLURM cluster:
```bash
nextflow run main.nf \
  -profile apptainer,slurm \
  --samplesheet /path/to/samples.csv \
  --outdir /path/to/results \
  --slurm_account YOUR_ACCOUNT \
  -resume
```

Or edit `nextflow.sh` (set `SAMPLESHEET`, `OUTDIR`, `YOUR_ACCOUNT`, `YOUR_PARTITION`) and submit:
```bash
sbatch nextflow.sh
```

---

## Overall pipeline

### 1. Basecalling

The basecalling and demultiplexing step are performed by the user outside of the pipeline.

### 2. Nanopore reads quality metrics

[Nanocomp](https://github.com/wdecoster/nanocomp) is used to compute Nanopore read metrics (e.g. Median Read Length, Read N50, Median Read Quality). Those metrics are computed for each barcode on the raw reads included in the basecalled fastq files.

### 3. Flye assembly and polishing

- The Nanopore reads are assembled using the software [Flye](https://github.com/fenderglass/Flye). By default, Flye will use a reduced coverage for initial disjointig assembly to speed up the assembly (config parameter `flye_args = "--asm-coverage 100"`).
- The draft assemblies are subsequently polished using [Medaka](https://github.com/nanoporetech/medaka). The model parameter selected to run Medaka (e.g. `r1041_e82_400bps_sup_v5.0.0`) must correspond to the model used for the basecalling (e.g. `dna_r10.4.1_e8.2_400bps_sup.cfg`).

### 4. Assembly quality assessment with QUAST

The software [QUAST](https://quast.sourceforge.net/quast.html) is used to compute genome assembly metrics on the polished assemblies.

### 5. Assembly quality assessment with CheckM

The software [CheckM](https://github.com/Ecogenomics/CheckM) v1 (command [lineage_wf](https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow)) is used to compute genome assembly completeness and contamination, based on the presence or absence of marker genes.

### 6. Sylph taxonomy classification

Nanopore reads are used as input to the taxonomy classifier [Sylph](https://sylph-docs.github.io). Both the GTDB R226 and RefSeq Fungi databases were downloaded from https://sylph-docs.github.io/pre%E2%80%90built-databases/.

### 7. LPS typing using Kaptive

The LPS type of the sample is obtained using the software [Kaptive](https://kaptive.readthedocs.io/en/latest/) v3. The genome assemblies are used as input to this tool. The 9 LPS database is used but can be modified (config parameter `kaptive_db_9lps`).

### 8. Variant calling using Clair3

- The reads are mapped to the reference LPS type sequence identified by Kaptive using the sequence alignment program [Minimap2](https://github.com/lh3/minimap2).
- Then [Clair3](https://github.com/HKU-BAL/Clair3/) is used to call variants in the reads as compared to the reference LPS type sequence. The Clair3 model must be selected to match the ONT platform and basecalling model used (e.g. `r1041_e82_400bps_sup_v500`). The list of Clair3 models is available at https://github.com/nanoporetech/rerio/tree/master/clair3_models.
- Then [SnpEff](https://pcingola.github.io/SnpEff/) is used to annotate and predict the effects of the variants on genes and proteins (such as amino acid changes).
- Finally, [SnpSift](https://pcingola.github.io/SnpEff/#snpsift) is used to extract the variants predicted to have a high impact on the protein (frameshift and stop_gained variants).

### 9. MLST typing

The software [mlst](https://github.com/tseemann/mlst) is used to scan the genome assemblies against the PubMLST typing scheme `pmultocida_2` by default (RIRDC). The typing scheme can be modified by specifying `--mlst_scheme` (e.g. `--mlst_scheme pmultocida`).

### 10. petG detection

The petG reference sequence is searched against the same assembly used for Kaptive typing using BLAST. A sample is reported as petG-positive when a hit has a genomic span greater than 1570 bp and at least 95% identity. Matching hit sequences are written to the sample output directory.

### 11. Subtype report

The pipeline generates a subtype report file (`10_ONT_subtype_report.tsv`) summarising the variants found in the subtype database. To be reported, the variant identified by Clair3 must be present in the subtype database with the following conditions:
- the variant must be identified at the same position in the reference sequence, and
- both the reference allele and the alternate allele must match their corresponding allele from the variant in the database.

If the subtype database includes `PHENOTYPE_DEFAULT` and `PHENOTYPE_MULTIPLE_SUBTYPES` columns after the `GENE` column, the report also assigns a phenotype for each subtype. Phenotype descriptions are added from `phenotype_lookup.tsv` when this file is present in the LPS reference database directory.

### 12. Genome annotation using Bakta

The software [Bakta](https://github.com/oschwengers/bakta) is used to annotate the genome assemblies. The default database is v6.0 from 2025-02-24, https://zenodo.org/records/14916843.

### 13. Antimicrobial resistance genes

The software [AMRFinderPlus](https://github.com/ncbi/amr) is used to identify AMR genes in the genome assemblies. The default database is version 2025-03-25.1 downloaded using `amrfinder_update`.

---

## Database setup

The pipeline uses two types of databases:

**Included in this repository** (`databases/` directory):

| Database | Contents | Notes |
|----------|----------|-------|
| LPS references (L1–L9) | `databases/LPS/` — FASTA/GenBank files, `reference_LPS.txt`, `LPS_subtype_database_v2.txt`, `phenotype_lookup.tsv`, `petG_X73_NZ_CM001580.fasta` | Required for LPS typing and variant calling |
| Kaptive3 LPS DB | `databases/kaptive3_LPS_db_v1/9lps.gbk`, `9lps.logic` | Required for LPS typing |

**Downloaded separately** — the pipeline can download these automatically on first run (set the corresponding `--skip_download_*_db false`), or you can download them manually and point to them with `--*_db`:

| Database | Required for | Default path | Tested version | Source |
|----------|-------------|-------------|---------------|--------|
| CheckM | Assembly quality | `<pipeline_dir>/databases/checkm_data_2015_01_16` | `checkm_data_2015_01_16` | [CheckM installation docs](https://github.com/Ecogenomics/CheckM/wiki/Installation) |
| Sylph GTDB + Fungi | Taxonomy classification | `<pipeline_dir>/databases/sylph/` | GTDB R226, Fungi RefSeq 2024-07-25 | [Sylph pre-built databases](https://sylph-docs.github.io/pre%E2%80%90built-databases/) |
| Clair3 model | Variant calling | `<pipeline_dir>/databases/clair3_models/r1041_e82_400bps_sup_v500` | `r1041_e82_400bps_sup_v500` | [Rerio models](https://github.com/nanoporetech/rerio/tree/master/clair3_models) — manual download |
| Bakta | Genome annotation | `<pipeline_dir>/databases/bakta/db` | v6.0 (2025-02-24) | [Zenodo 10.5281/zenodo.14916843](https://zenodo.org/records/14916843) |
| AMRFinderPlus | AMR gene identification | `<pipeline_dir>/databases/amrfinderplus/2025-03-25.1` | `2025-03-25.1` | [NCBI AMRFinderPlus](https://github.com/ncbi/amr/wiki/AMRFinderPlus-database) |

> **Reproducibility note:** The `--skip_download_*_db false` flags fetch the **latest** available version of each database, which may differ from the tested versions listed above and could produce different results. For reproducible analyses, use the tested versions by downloading them manually and setting the corresponding `--*_db` parameter.

---

## Step-by-step user guide

### 1. Clone the repository

```bash
git clone https://github.com/vmurigneu/LPS_typing_ONT.git
cd LPS_typing_ONT
```

The repository contains:
- **`nextflow.config`** — default parameters and container settings. An example is available [here](https://github.com/vmurigneu/LPS_typing_ONT/blob/main/nextflow.config).
- **`main.nf`** — pipeline code (not user-modifiable).
- **`nextflow.sh`** — SLURM submission template. Edit it to set your paths and account name before submitting.

### 2. Obtain databases

See [Database setup](#database-setup) above.

### 3. Prepare the samplesheet (CSV)

The samplesheet is a comma-separated file listing each sample and its FASTQ file. The header line must match the example below.

```
sample_id,long_fastq
PM1947,fastq/barcode17.simplex_duplex.fastq.gz
PM1422,fastq/barcode18.simplex_duplex.fastq.gz
```

`long_fastq` paths are resolved relative to the **Nextflow launch directory** (the directory where you run `nextflow run main.nf`), not relative to the samplesheet file's location.

### 4. Run the pipeline

**Local / non-SLURM:**
```bash
nextflow run main.nf \
  -profile apptainer \
  --samplesheet /path/to/samples.csv \
  --outdir /path/to/results \
  -resume
```

**SLURM cluster:**
```bash
nextflow run main.nf \
  -profile apptainer,slurm \
  --samplesheet /path/to/samples.csv \
  --outdir /path/to/results \
  --slurm_account YOUR_ACCOUNT \
  -resume
```

Or edit `nextflow.sh` (set `SAMPLESHEET`, `OUTDIR`, `YOUR_ACCOUNT`, `YOUR_PARTITION`) and submit:
```bash
sbatch nextflow.sh
```

**Typing-only run** (skip most annotation steps):
```bash
nextflow run main.nf \
  -profile apptainer,slurm \
  --samplesheet /path/to/samples.csv \
  --outdir /path/to/results \
  --slurm_account YOUR_ACCOUNT \
  --skip_polishing true \
  --skip_mlst true \
  --skip_quast true \
  --skip_checkm true \
  --skip_bakta true \
  --skip_amrfinder true \
  -resume
```

Use `-resume` to restart an interrupted run without recomputing completed steps.

---

## Example data

To test the pipeline, test data containing 20,000 Nanopore reads for two samples is provided:

File | Description
---|---
[fastq/barcode17.simplex_duplex.test.fastq.gz](https://github.com/vmurigneu/LPS_typing_ONT/blob/main/fastq/barcode17.simplex_duplex.test.fastq.gz) | ONT reads for barcode17
[fastq/barcode18.simplex_duplex.test.fastq.gz](https://github.com/vmurigneu/LPS_typing_ONT/blob/main/fastq/barcode18.simplex_duplex.test.fastq.gz) | ONT reads for barcode18
[samplesheet/samples_test.csv](https://github.com/vmurigneu/LPS_typing_ONT/blob/main/samplesheet/samples_test.csv) | Samplesheet for the test run
[nextflow.sh](https://github.com/vmurigneu/LPS_typing_ONT/blob/main/nextflow.sh) | SLURM submission template

Run the test data:
```bash
nextflow run main.nf \
  -profile apptainer \
  --samplesheet samplesheet/samples_test.csv \
  --outdir results_test \
  -resume
```

---

## Optional parameters

Some parameters can be added to the command line in order to include or skip some steps and modify some parameters:

2. Nanopore reads quality metrics:
* `--skip_nanocomp`: skip the Nanocomp step (default=false)
* `--nanocomp_threads`: number of threads for the Nanocomp step (default=4)

3. Genome assembly and polishing:
* `--skip_assembly`: skip the assembly step (default=false). Note: it is not recommended to skip assembly as many downstream steps depend on assembly results.
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
* `--skip_download_checkm_db`: skip downloading the CheckM database (default=true — assumes the database is already present at `--checkm_db`)
* `--checkm_db`: path to the CheckM database folder (default=`<pipeline_dir>/databases/checkm_data_2015_01_16`)

6. Sylph taxonomy classification:
* `--skip_sylph`: skip the Sylph classification step (default=false)
* `--skip_download_sylph_db`: skip downloading the Sylph databases (default=true — assumes databases are already present at `--sylph_db`)
* `--sylph_db_gtdb_file` and `--sylph_db_fungal_file`: URLs for Sylph GTDB and Fungi RefSeq database files to download
* `--sylph_tax_gtdb_metadata` and `--sylph_tax_fungal_metadata`: URLs for Sylph-tax metadata files to download (must match the database files)
* `--sylph_db`: glob pattern for pre-downloaded Sylph database files (default=`<pipeline_dir>/databases/sylph/*.syldb`)
* `--sylph_threads`: number of threads for the Sylph classification step (default=6)

7. LPS typing using Kaptive:
* `--skip_kaptive3`: skip the Kaptive typing step (default=false). Note: skipping Kaptive automatically skips variant calling.
* `--kaptive_db_9lps`: path to the Kaptive database file (default=`<pipeline_dir>/databases/kaptive3_LPS_db_v1/9lps.gbk`)

8. Variant calling using Clair3:
* `--skip_clair3`: skip the variant calling step (default=false)
* `--minimap_threads`: number of threads for the Minimap2 mapping step (default=6)
* `--clair3_threads`: number of threads for the Clair3 variant calling step (default=4)
* `--clair3_model`: path to the Clair3 model folder (default=`<pipeline_dir>/databases/clair3_models/r1041_e82_400bps_sup_v500`)
* `--clair3_args`: Clair3 optional parameters (default="--haploid_sensitive"), see [available parameters](https://github.com/HKU-BAL/Clair3?tab=readme-ov-file#options)
* `--skip_snpeff`: skip the variant annotation step (default=false)
* `--reference_LPS_directory`: directory containing the reference LPS sequence files (default=`<pipeline_dir>/databases/LPS`). This directory should contain `reference_LPS.txt`, `LPS_subtype_database_v2.txt`, `petG_X73_NZ_CM001580.fasta`, and optionally `phenotype_lookup.tsv`.

9. MLST typing:
* `--skip_mlst`: skip the MLST typing step (default=false)
* `--mlst_scheme`: MLST typing scheme (default="pmultocida_2")

10. petG detection:
* `--skip_petg`: skip petG detection with BLAST (default=false)
* `--petg_threads`: number of threads for the petG BLAST step (default=2)
* `--petg_min_length`: minimum genomic hit span for petG presence; hits must exceed this value (default=1570)
* `--petg_min_identity`: minimum percent identity for petG presence (default=95)

11. Report:
* `--reference_LPS_directory`: path to the subtype database directory (default=`<pipeline_dir>/databases/LPS`). This directory should contain `LPS_subtype_database_v2.txt` and, for phenotype descriptions, `phenotype_lookup.tsv`.

12. Genome annotation using Bakta:
* `--skip_bakta`: skip the genome annotation step (default=false)
* `--bakta_threads`: number of threads for the Bakta step (default=8)
* `--skip_download_bakta_db`: skip downloading the Bakta database (default=true)
* `--bakta_db`: path to the Bakta database files (default=`<pipeline_dir>/databases/bakta/db`)
* `--bakta_args`: Bakta optional parameters (default includes trusted protein sequences from the LPS database)

13. AMR gene identification using AMRFinderPlus:
* `--skip_amrfinder`: skip the AMR gene identification step (default=false)
* `--skip_download_amrfinder_db`: skip downloading the AMRFinderPlus database (default=true)
* `--amrfinder_db`: path to the AMRFinderPlus database files (default=`<pipeline_dir>/databases/amrfinderplus/2025-03-25.1`)
* `--amrfinder_args`: AMRFinderPlus optional parameters (default="")

---

## Structure of the output folders

The pipeline will create several folders corresponding to the different steps of the pipeline.
* **2_nanocomp:** Quality control of the Nanopore reads
  * Full Nanocomp report in HTML format with plots and metrics table (NanoComp-report.html)
  * Nanocomp summary text file (NanoStats.txt)

The main output folder (`--outdir`) will contain one folder per sample (named by the `sample_id` column in the samplesheet).

Each sample folder will contain the following folders:
* **3_assembly:** Flye assembly output files (.fasta, .gfa, .gv, .info.txt), see [details](https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#-flye-output). The final polished assembly fasta file is `sample_id_flye_polished.fasta`.
* **4_quast:** QUAST output report file (`sample_id_report.tsv`).
* **5_checkm:** CheckM output file (`sample_id_checkm_lineage_wf_results.tsv`).
* **6_sylph:** Sylph taxonomy classification results for the Nanopore reads, see [details](https://sylph-docs.github.io/Output-format/) and [more details](https://sylph-docs.github.io/sylph-tax-output-format/)
  * Sylph list of genomes detected: `sample_id_sylph_profile.tsv`
  * Sylph-tax individual sequence abundances combined: `sample_id_merged_sequence_abundance.tsv`
  * Sylph-tax individual taxonomic abundances combined: `sample_id_merged_taxonomic_abundance.tsv`
* **7_kaptive_v3:** Kaptive output files, see [details](https://kaptive.readthedocs.io/en/latest/Outputs.html)
  * LPS type results (`sample_id_kaptive_results.tsv`)
  * LPS sequence in fasta format (`sample_id_flye_polished_kaptive_results.fna`)
* **8_clair3:** Mapping files and variant calling results:
  * Minimap2 mapping file in bam format (`sample_id_minimap2_mapped.bam` and `.bai` index). Unmapped reads were excluded.
  * Minimap2 mapping statistics (`sample_id_minimap2_flagstat.txt`)
  * Clair3 variants (`sample_id_clair3.vcf`)
  * Clair3 variants annotated by SnpEff (`sample_id_clair3.snpeff.vcf`)
  * Frameshift and stop_gained Clair3 variants (`sample_id_clair3.snpeff.high_impact.vcf`)
* **9_mlst:** MLST typing output file (`sample_id_mlst.csv`)
* **13_petG:** petG BLAST output files
  * Accepted petG hit sequences (`sample_id_petG_hits.fasta`)
  * BLAST tabular output for all hits (`sample_id_petG_blast.tsv`)
  * BLAST tabular output for accepted hits (`sample_id_petG_blast.filtered.tsv`)
  * petG presence summary (`sample_id_petG_summary.tsv`)
* **10_report:** Summary of results for all samples
  * Flye assembly statistics: assembly coverage, number of contigs, assembly size (`3_ONT_flye_stats.tsv`)
  * QUAST combined report file (`4_ONT_quast_report.tsv`)
  * CheckM results (`5_ONT_checkm_lineage_wf_results.tsv`)
  * Sylph taxonomy results:
    - Abundance of *P. multocida* reads, and information on the most abundant species (if not *P. multocida*): `6_ONT_sylph_summary.tsv`
  * Kaptive results (`7_ONT_kaptive_results.tsv`)
  * Clair3 variants results:
    - All variants: `8_ONT_clair3_snpeff.vcf`
    - Only variants predicted to have a high impact on the protein: `8_ONT_clair3_snpeff_high_impact.vcf`
  * MLST results (`9_ONT_mlst.csv`)
  * Subtype results summarising the variants found in the subtype database (`10_ONT_subtype_report.tsv`). Columns:
    - SAMPLE: sample identifier
    - MLST: MLST sequence type number
    - TYPE: LPS type assigned by Kaptive when confidence is Typeable, otherwise untypeable
    - SUBTYPE: LPS subtype assigned by the pipeline (using the subtype database)
    - VARTYPE: description of the variant
    - ISOLATE_DATABASE: reference isolate from the subtype database that contained that variant
    - CHROM: name of the reference sequence for the LPS locus type
    - POS: variant position in the reference sequence for the LPS locus type
    - REF: reference allele sequence present in the LPS reference sequence
    - ALT: alternate allele sequence identified in the sample
    - GENE: gene containing the variant
    - PHENOTYPE: phenotype label assigned from the subtype database, if available
    - PHENOTYPE_DESCRIPTION: phenotype description from `phenotype_lookup.tsv`, if available
    - PETG_PRESENT: "yes" when a qualifying petG BLAST hit is present, blank otherwise
    - NOTE: note from the subtype database
  * AMRFinderPlus results (`12_ONT_amrfinder.tsv`)
* **11_bakta:** Bakta genome annotation output files, see [details](https://github.com/oschwengers/bakta?tab=readme-ov-file#output).
  * Annotations & sequences in (multi) GenBank format (`sample_id_bakta.gbff`)
  * Inference metrics as TSV (`sample_id_bakta.inference.tsv`)
  * Annotation summary in text format (`sample_id_bakta.txt`)
* **12_amrfinder:** AMRFinderPlus output file (`sample_id_amrfinder.tsv`), see [details](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#output-format).

---

## Running the workflow in assembly mode for other organisms

The default parameters are suited for *Pasteurella multocida*. The LPS typing and variant calling are specific to *Pasteurella multocida*. To use the workflow to assemble another species:
* `--genome_size`: estimated genome size (default="2.3M")
* `--mlst_scheme`: MLST typing scheme (default="pmultocida_2")
* `--skip_kaptive3 true`: skip LPS typing (this also skips variant calling automatically)

---

## Troubleshooting

Common things to check if an error occurs:
* Read the error message in the log file
* The samplesheet has the correct header line (see [User guide](#step-by-step-user-guide))
* The FASTQ files in the samplesheet are present at the paths listed, resolved relative to the Nextflow launch directory
* The `databases/` folder is present in the cloned pipeline repository and not empty (see [Database setup](#database-setup))
