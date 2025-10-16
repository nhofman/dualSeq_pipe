# Dual RNA-Seq analysis pipeline

## Description
In a dual RNA-Seq approach human HuH7 cells were infected with different pathogenic viruses. Samples were taken after predefined times of infection.
The pipeline was build to analyse data from (these) infection experiments where samples are a mix of host and pathogenic reads. A python script is used to build and execute a customized snakemake workflow based on the analysis modules and parameters given by the user. 
Available analysis modules are:
- QC
	- fastqc
- Preprocessing
	- none
	- fastp
- Mapping
	- combined {STAR}
- Gene_expression_analysis
	- count_table {featureCounts}
	- deseq2
- Report
	- bamCoverage {bigwig}
	- multiqc
- Variant_analysis
	- preprocessing
	- variant_calling

## Prerequisites
- Mamba/Conda
- Snakemake >= 7
- Python > 3.10
- Pyyaml

The required tools need to be installed before running the pipeline. If conda or mamba is preexisting on the system, the file dualseq.yaml can be used to create a conda/mamba environment containing all necessary tools and dependencies.

## Execution
The wrapper script can be executed with the following command:
```
usage: dualSeq.py --data DATA_FILE --pipeline-config PIPELINE_CONFIG --outdir OUTPUT_FOLDER [--profile PROFILE] [-t CORES] [-j JOBS] [--use-conda] [--conda-frontend {mamba,conda}]
                  [--conda-prefix CONDA_PREFIX] [--conda-create-envs-only] [--latency-wait LATENCY] [--other [OTHER ...]] [-v] [--verbose] [-h]

Dual RNA-Seq analysis pipeline

Required arguments:
  --data DATA_FILE      Path to comma- or tab-separated file that lists all input data.
  --pipeline-config PIPELINE_CONFIG
                        Path to YAML file that defines rules and rule specific parameters.
  --outdir OUTPUT_FOLDER, -o OUTPUT_FOLDER
                        Path to output folder.

Other arguments:
  --profile PROFILE     Path to folder containing profile 'config.yaml' for snakemake configuration. Can be used to set default values for command line options, e.g. cluster submission command. See
                        also: https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
  -t CORES, --cores CORES
                        Number of threads/cores (Default: 1).
  -j JOBS, --jobs JOBS  Maximal number of parallel jobs send to the cluster (Default: 1). Only used in cluster mode.
  --use-conda           Run job in conda environment, if defined in rule.
  --conda-frontend {mamba,conda}
                        Choose frontend for installing environmnents ['conda', 'mamba']. (Default: mamba)
  --conda-prefix CONDA_PREFIX
                        Path to folder where conda directories are created, can be a path relative to invocation dir or an absolute path.
  --conda-create-envs-only
                        Only create job-specific environments and exit. --use-conda has to be set.
  --latency-wait LATENCY
                        Seconds to wait before checking if all files of a rule were created (Default: 3). Should be increased if using cluster mode.
  --other [OTHER ...]   Add additional snakemake command line options, e.g. 'dry-run' ('--' is automatically placed in front.)
  -v, --version         Show program's version number and exit
  --verbose             Print debugging output
  -h, --help            Show this help message and exit

```

## Input

### Required:

#### Data file (--data): 

A comma- or tab-separated file that describes all input data. The samples are listed line by line.

- **name**: unique sample name 
- **reads**: path to file containing reads; *if single reads* 
- **forward_reads**: path to file containing forward reads; *if paired end*
- **reverse_reads**: path to file containing reverse reads; *if paired end*
- **condition**: sample condition, e.g. treatment
- **virus_genome**: path to virus genome file (fasta)
- **virus_genome_annotation**: path to virus annotation file (gff/gtf)

Example:

name | reads | condition | virus_genome | virus_genome_annotation
-----|-------|-----------|--------------|------------------------
EBOV_12h_1 | /path/to/reads.fq | EBOV_12h | /path/to/EBOV.fasta | /path/to/EBOV.gtf
H1N1_24h_1 | /path/to/reads.fq | H1N1_24h | /path/to/H1N1.fasta | /path/to/H1N1.gtf
Mock_12h_1 | /path/to/reads.fq | Mock_12h | | 
Mock_12h_2 | /path/to/reads.fq | Mock_12h | | 
Mock_24h_1 | /path/to/reads.fq | Mock_24h | | 

#### Pipeline config file (--pipeline-config): 

A YAML file that defines the modules that will be executed and the rule specific parameters. For each module,there are module-specific YAML files that define mandatory and optional parameters. An example pipeline_config.yaml can be found in the test folder.
	
Example:
```
pipeline:
  paired_end: false

QC:
  module: "fastqc"

preprocessing:
  module: "fastp"

mapping:
  module: "star"
  
  star:
    host_genome: "/path/to/human/genome.fasta"
    host_gtf: "/path/to/human/genome.gtf"
    star_index_dir: "/path/to/human/star_idx/"

gene_expression_analysis:
  modules: ["count_table", "deseq2"]

  count_table:
    gff_feature_type_human: "gene"
    gff_feature_name_human: "gene_id"
    gff_human: "/path/to/human/genome.gtf"

  deseq2:
    comparisons_human: "/path/to/comparisons.txt"
    color: "/path/to/color.txt"
    rRNA: "/path/to/rRNA_genes.txt"

```
**comparisons file**: A tab-separated file that defines the comparisons for the differential gene expression analysis. The 1st column specifies the name of the numerator (e.g. treated) while the 2nd column gives the name of the denominator (e.g. untreated) for the log2 fold change calculation. The file path is specified in the config file.

**host genome**: Host reference genome {fasta}. File path needs to be specified in config file.

**host annotation**: Host reference annotation {gff/gtf}. File path needs to be specified in config file.

**color file**: A tab-separated file that is used to specify custom colors for each pathogen (pathogen	color). The file pah is set for the deseq2 module in the config file.

### Optional:

#### Snakemake profile: 
A YAML file to define default values for snakemake command line parameters. The profile is used to specify the cluster command and resources.

#### Conda/Mamba

The pipeline uses the Snakemake conda integration to provide the necessary software packages, if conda or mamba is installed on the system. The environments are defined as YAML files and are set via the conda directive for each rule. Predefined environment files are provided in the directory `conda_envs`.
To use the conda integration, the flag `--use-conda` has to be set. The conda frontend can be defined with the argument `--conda-frontend <{mamba, conda}>`, available options are mamba ('default') and conda. The conda environments are created in the `.snakemake` directory within the current working directory per default. The option `--conda-prefix <directory>` can be used to set a user-defined directory.

## Results
- **QC**
	- FastQC result for each sample (html)
- **Preprocessing**
	- Fastp result for each sample (html)
- **Mapping**
	- Alignment files for each sample 
		- Host mapping results (bam)
		- Pathogen mapping resluts (bam)
	- Mapping statistics
		- table (csv, xlsx)
		- plot (pdf, png)
- **Gene expression analysis**
	- featureCounts
		- Count table (.txt)
		- featureCounts statistic plot (pdf, png)
	- Differential gene expression analysis
		- Boxplots of raw and normalized read counts (pdf, png)
		- Sample correlation heatmap (pdf, png)
		- PCA (pdf, png)
		- DESeq2 result tables (csv, xlsx)
		- Volcano plot / MAplot (pdf, png)
		- Summary of differentially expressed genes (csv, xlsx, pdf, png)
	- Functional analysis
		- GSEA (csv, pdf, svg)
- **Report**
	- MultiQC report (html)
	- bamCoverage files (bigwig)
- **Variant analysis**
	- Variant call files (vcf)
	- VCF statistic plots (pdf, png)
	- Plot of variant positions on genome (pdf, png) 

## Citation

**Sequence data availability:**

The raw sequence data of the time-resolved virus infection study were deposited in the Short Read Archive (SRA) of the National Center for Biotechnology Information (NCBI) under the BioProject ID [PRJNA1074963](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1074963/).

**Paper on infection study:**

Hofmann N, Bartkuhn M, Becker S, Biedenkopf N, Böttcher-Friebertshäuser E, Brinkrolf K, Dietzel E, Fehling SK, Goesmann A, Heindl MR, Hoffmann S, Karl N, Maisner A, Mostafa A, Kornecki L, Müller-Kräuter H, Müller-Ruttloff C, Nist A, Pleschka S, Sauerhering L, Stiewe T, Strecker T, Wilhelm J, Wuerth JD, Ziebuhr J, Weber F, Schmitz ML.2024.Distinct negative-sense RNA viruses induce a common set of transcripts encoding proteins forming an extensive network. J Virol98:e00935-24.https://doi.org/10.1128/jvi.00935-24

## Contact

Nina Hofmann<br />
nina.hofmann@computational.bio<br />

Bioinformatics and Systems Biology<br />
Justus Liebig University Giessen<br />
Heinrich-Buff-Ring 58<br />
35392 Giessen<br />
Germany
