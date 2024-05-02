# Dual RNA-Seq analysis pipeline

## Description
In a dual RNA-Seq approach human HuH7 cells were infected with different pathogenic viruses. Samples were taken after predefined times of infection.
The pipeline was build to analyse data from (these) infection experiments where samples are a mix of host and pathogenic reads. A python script is used to build and execute a customized snakemake workflow based on the analysis modules and parameters given by the user. 
Available analysis modules are:
- QC
	- fastqc
- Preprocessing
	- fastp
- Mapping
	- combined {STAR}
- Analyses
	- count_table {featureCounts}
	- deseq2
- Summary
	- bamCoverage {bigwig}
	- multiqc
- Variant_analyses
	- Preprocessing
	- variant_calling

## Execution
```
usage: dualSeq.py --groups GROUPS_FILE --config CONFIG_FILE --output OUTPUT_FOLDER [--profile PROFILE] [-t CORES] [--cluster-nodes CLUSTER_NODES] [--use-conda] [--conda-frontend {mamba,conda}] [--conda-create-envs-only] [--latency-wait LATENCY] [-v] [--other OTHER] [--verbose] [-h]

Differential gene expression pipeline for dual RNA-Seq analysis

Required arguments:
  --groups GROUPS_FILE  Path to comma- or tab-separated file that lists all input data.
  --config CONFIG_FILE  Path to yaml file that defines rules and rule specific parameters.
  --output OUTPUT_FOLDER
                        Path to output folder.

Other arguments:
  --profile PROFILE     Path to folder containing profile 'config.yaml' for
                        snakemake configuration. Can be used to set default
                        values for command line options, e.g. cluster
                        submission command. See also: https://snakemake.readth
                        edocs.io/en/stable/executing/cli.html#profiles
  -t CORES, --cores CORES
                        Number of threads/cores (Default: 1). Defines locale
                        cores in cluster mode
  --cluster-nodes CLUSTER_NODES
                        Maximal number of parallel jobs send to the cluster
                        (Default: 1). Only used in cluster mode.
  --use-conda           Run job in conda environment, if defined in rule.
  --conda-frontend {mamba,conda}
                        Choose frontend for installing environmnents ['conda',
                        'mamba']. (Default: mamba)
  --conda-create-envs-only
                        Only create job-specific environments and exit. --use-
                        conda has to be set.
  --latency-wait LATENCY
                        Seconds to wait before checking if all files of a rule
                        were created (Default: 3). Should be increased if
                        using cluster mode.
  -v, --version         Show program's version number and exit
  --other OTHER         Add additional snakemake command line options, e.g.
                        'dry-run'. '--' is automatically placed in front.
  --verbose             Print debugging output
  -h, --help            Show this help message and exit


```

## Prerequisites
- Snakemake
- Mamba/Conda
- Python

The pipeline uses the Snakemake conda integration to provide the necessary software packages, if conda or mamba is installed on the system. The environments are defined as yaml files and can be set via the conda directive for each rule.

## Input

#### Required:
- **groups file**: A comma- or tab-separated file that describes all input data. The samples are listed line by line.

| Column names | Description |
|--------------|--------------|
| name | unique sample name |
| reads | path to file containing reads; *if single reads* |
| forward_reads | path to file containing forward reads; *if paired end* |
| reverse_reads | path to file containing reverse reads; *if paired end* |
| condition | sample condition, e.g. treatment |
| virus_genome | path to virus genome file (fasta) |
| virus_genome_annotation | path to virus annotation file (gff/gtf) |

- **config file**: A yaml file that defines the rules that will be executed and the rule specific parameters. An example pipeline_config.yaml can be found in the test folder.
	- **comparisons file**: A tab-separated file that defines the comparisons for the differential gene expression analysis. The 1st column specifies the name of the numerator (e.g. treated) while the 2nd column gives the name of the denominator (e.g. untreated) for the log2 fold change calculation. The file path is specified in the config file.
	- **host genome**: Host reference genome {fasta}. File path needs to be specified in config file.
	- **host annotation**: Host reference annotation {gff/gtf}. File path needs to be specified in config file.
}
 
#### Optional:
- **snakemake profile**: A yaml file to define default values for snakemake command line parameters. The profile is used to specify the cluster command and resources.
- **color file**: A tab-separated file that is used to specify custom colors for each pathogen (pathogen	color). The file pah is set for the deseq2 module in the config file.


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
		- table (tsv)
		- plot (svg, png)
- **Analyses**
	- Count table (.txt)
	- DESeq2 result tables (csv, xlsx)
	- Sample correlation heatmap (pdf)
	- PCA (pdf, svg)
	- featureCounts statistic plot (pdf)
	- Volcano plot / MAplot (pdf)
	- Summary of differentially expressed genes (tsv, png, svg)
	- GSEA (csv, pdf, svg)
- **Summary**
	- MultiQC (html)
	- bamCoverage files (bigwig)
- **Variant analysis**
	- variant call files (vcf)

## Contact

Nina Hofmann<br />
nina.hofmann@computational.bio<br />

Bioinformatics and Systems Biology<br />
Justus Liebig University Giessen<br />
Heinrich-Buff-Ring 58<br />
35392 Giessen<br />
Germany
