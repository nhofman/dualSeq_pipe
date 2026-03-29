# Dual RNA-Seq analysis pipeline

## Idea
In a dual RNA-Seq approach human HuH7 cells were infected with different pathogenic viruses. Samples were taken after predefined times of infection. Total RNA was extracted and sequenced on an Illumina HiSeq4000.

A pipeline was developed for the analysis of data from such dual RNA-Seq infection experiments, in which samples consist of a mixture of host and pathogenic reads. The objective of this study was to analyze the human cell response to viral infections and to conduct a mutation analysis of the viral sequences. Consequently, the pipeline is divided into three segments: joint analyses, host-specific analyses, and virus-specific analyses. A Python script is used to build and execute a customized Snakemake workflow based on the analysis modules and parameters specified by the user. A workflow comprises a series of modules and submodules.

Available analysis modules are:

#### Joint analyses
---
- QC
	- fastqc
- Preprocessing
	- none
	- fastp
- Mapping
	- combined {STAR, joint host + virus reference}
  - consecutive {1. STAR, whole sample (host reference), 2. Bowtie2, unmapped reads from STAR (virus reference)}
- Report
	- bamCoverage {bigwig}
	- multiqc

#### Host-specific analyses
---
- Gene_expression_analysis
	- count_table {featureCounts}
	- deseq2
	- functional_analysis

#### Virus-specific analyses
---
- Variant_analysis
	- preprocessing
	- preprocessing_dedup
	- variant_calling


## Prerequisites
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) >= 7.32
- Python >= 3.10
- Pyyaml
- ([Mamba](https://mamba.readthedocs.io/en/latest/index.html)/[Conda](https://docs.conda.io/en/latest/))

The required tools need to be installed before running the pipeline. If conda or mamba is preexisting on the system, the file dualseq.yaml can be used to create a conda/mamba environment containing the required tools and dependencies.

```
mamba env create -f dualseq.yaml 
mamba activate dualseq
```
or 
```
mamba create -n <env_name> -f dualseq.yaml
mamba activate <env_name>
```

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

#### 1. Data file (--data): 

A comma- or tab-separated file that describes all input data. The samples are listed line by line. File paths can either be absolute paths or paths relative to this data file. In case of uninfected samples no virus genome and annotation are specified.

- **name**: unique sample name 
- **reads**: path to the file containing the sequencing reads; *if single reads* 
- **forward_reads**: path to the file containing forward reads; *if paired end*
- **reverse_reads**: path to the file containing reverse reads; *if paired end*
- **condition**: sample condition, e.g. treatment or control
- **virus_genome**: path to virus genome file (fasta)
- **virus_genome_annotation**: path to virus annotation file (gff/gtf)

###### Example for single reads:

name | reads | condition | virus_genome | virus_genome_annotation
-----|-------|-----------|--------------|------------------------
EBOV_12h_1 | /path/to/reads.fq | EBOV_12h | /path/to/EBOV.fasta | /path/to/EBOV.gtf
H1N1_24h_1 | /path/to/reads.fq | H1N1_24h | /path/to/H1N1.fasta | /path/to/H1N1.gtf
Mock_12h_1 | /path/to/reads.fq | Mock_12h | | 
Mock_12h_2 | /path/to/reads.fq | Mock_12h | | 
Mock_24h_1 | /path/to/reads.fq | Mock_24h | | 

#### 2. Pipeline config file (--pipeline-config): 

The pipeline config file is a YAML file that defines the modules that will be executed and the rule specific parameters. For each available module, there are module-specific YAML files that define mandatory and optional parameters. These can be found at `snakefiles/<module_name>/<submodule_name/<submodule_name>.yaml`. The parameter `paired-end` is used to define the sequencing strategy (paired-end or single read), which will trigger the use of specific module-files.
	
###### Example:
```
pipeline:
  paired_end: false

QC:
  module: "fastqc"

preprocessing:
  module: "fastp"

mapping:
  module: "combined"
  
  combined:
    host_genome: "/path/to/human/genome.fasta"
    host_gtf: "/path/to/human/genome.gtf"
    star_index_dir: "/path/to/human/star_idx/"

gene_expression_analysis:
  modules: ["count_table", "deseq2"]

  count_table:
    gff_feature_type_human: "exon"
    gff_feature_name_human: "gene_id"
    gff_human: "/path/to/human/genome.gtf"

  deseq2:
    comparisons_human: "/path/to/comparisons.txt"
    color: "/path/to/color.txt"
    rRNA: "/path/to/rRNA_genes.txt"

```

**host_genome** (required): Path to host reference genome {fasta}. The file path can be specified as an absolute path or relative to the pipeline config file.

**host_gtf** (required): Path to host reference annotation {gff/gtf}. The file path can be specified as an absolute path or relative to the pipeline config file.

**comparisons_human file** (required): A tab-separated file that defines the comparisons for the differential gene expression analysis. The 1st column specifies the name of the numerator (e.g. treated) while the 2nd column gives the name of the denominator (e.g. untreated) for the log2 fold change calculation. The file path can be specified as an absolute path or relative to the pipeline config file.

**color file** (optional): A tab-separated file that is used to specify custom colors for each pathogen (pathogen	color). The file path is set for the deseq2 module in the config file. Colors can either be name or HEX code. The file path can be specified as an absolute path or relative to the pipeline config file.

**rRNA file** (optional): A file containing a list of gene names that should be removed prior to the normalization of feature counts, e.g., rRNA genes. The file path can be specified as an absolute path or relative to the pipeline config file.

### Optional:

#### Snakemake profile: 

Snakemake allows users to define workflow profiles that specify the default values for Snakemake command-line parameters, such as the cluster command and default resources (https://snakemake.readthedocs.io/en/stable/executing/cli.html#executing-profiles). The profile is defined in a YAML file, that has to be named config.yaml or config.vX+.yaml, where X is the minimum supported Snakemake version. The location of the config file is specified with the option `--profile <profile_folder>`.

###### Example:

```
cluster: sbatch --mem={resources.mem_mb} -c {threads} -p {resources.partition} -o .snakemake/log/%j.slurm.out  -e .snakemake/log/%j.slurm.err
default-resources: [mem_mb=500, partition=bcf]
printshellcmds: True
```

### Conda/mamba usage:

The pipeline uses the Snakemake conda integration to provide the necessary rule-specific software packages, as long as conda or mamba is installed on the system. The environments are defined as YAML files and are set via the conda directive for each rule. Predefined environment files are provided in the directory `conda_envs`.
To use the conda integration, the flag `--use-conda` has to be set. The conda frontend can be defined with the argument `--conda-frontend <{mamba, conda}>`, available options are mamba ('default') and conda. The conda environments are created in the `.snakemake` directory within the current working directory per default. The option `--conda-prefix <directory>` can be used to set a user-defined directory. If the conda environments should be created without running the pipeline, the flag `--conda-create-envs-only` can be specified.

In the event that the conda integration is not utilized, it is up to the user to ensure that the necessary analysis tools are accessible to the pipeline. The required tools are specified in the respective YAML files, which are stored in the directory `conda_envs`.

### Test:

The pipeline can be tested with the exemplary input data and pipeline files provided in the `test\test_single` directory. The folder contains a pipeline configuration file, sample overview file, the comparison.txt and a snakemake profile. A virus reference genome and annotation can be found in `references`, while the human reference files have to be prepared by the user. The pipeline can be executed with the file `run_pipeline.sh`.

## Results
All analysis results can be found in the output directory specified by the user (`--outdir`). For each submodule, a separate folder is created within the respective module folder. Available results include:
- **QC**
	- FastQC result for each sample (html)
- **Preprocessing**
	- Preprocessed reads (fq.gz)
	- Fastp result for each sample (html)
- **Mapping**
	- Alignment files for each sample 
		- Host mapping results (bam)
		- Pathogen mapping results (bam)
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

## Custom scripts

The folder `custom_scripts` contains custom R scripts that were used to compare the differential gene expression between time points and viruses. The necessary R packages can be found in the respective scripts. The comparisons are based on the DESeq2 analyses of the dualSeq pipeline. 

- `compare_time_points.R`:
	- Compare DEGs between time points by generating UpSet plots (for each virus)
	- Compare DEGs for *virus_24h vs Mock 24h*, *virus_24h vs virus_BPL*, *virus_BPL vs Mock_BPL*
- `compare_viruses.R`:
	- Compare DEGs across specified viruses (*common host reponse*)
		- UpSet plot
	- ComplexHeatmap of common host response
	- Overrepresentation analysis of common host response

## Interactive visualization 

An interactive visualization of the analysis results of the virus infection study can be viewed on the Analysis Dashboard for Virus-Induced CEll Response based on RNA-Seq data (**ADVICER**). ADVICER is an R Shiny app that is publicly available under https://advicer.computational.bio/. Link to GitHub Repo: https://github.com/nhofman/ADVICER


## Citation

**Sequence data availability:**

The raw sequence data of the time-resolved virus infection study were deposited in the Short Read Archive (SRA) of the National Center for Biotechnology Information (NCBI) under the BioProject ID [PRJNA1074963](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1074963/).

**Paper on infection study:**

> Hofmann N, Bartkuhn M, Becker S, Biedenkopf N, Böttcher-Friebertshäuser E, Brinkrolf K, Dietzel E, Fehling SK, Goesmann A, Heindl MR, Hoffmann S, Karl N, Maisner A, Mostafa A, Kornecki L, Müller-Kräuter H, Müller-Ruttloff C, Nist A, Pleschka S, Sauerhering L, Stiewe T, Strecker T, Wilhelm J, Wuerth JD, Ziebuhr J, Weber F, Schmitz ML.2024.Distinct negative-sense RNA viruses induce a common set of transcripts encoding proteins forming an extensive network. J Virol98:e00935-24.https://doi.org/10.1128/jvi.00935-24

## Contact

Nina Hofmann<br />
advicer@computational.bio<br />

Bioinformatics and Systems Biology<br />
Justus Liebig University Giessen<br />
35392 Giessen<br />
Germany
