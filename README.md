# VDJPuzzle2 - README #

TCR and BCR reconstruction from scRNA-seq data

## Setup

Download this repository and move into the VDJPuzzle directory. 

To create the conda environment run `conda create --name vdjpuzzle --file environment-explicit.txt`. If you don't have conda installed, you can find it [here](https://conda.io/docs/user-guide/install/index.html). 

Type `source activate vdjpuzzle` to activate the environment. 

Symlink or copy bin/vdjpuzzle into your path by copying the following command in the .bashrc file in your home directory and substituting path_to_vdjpuzzle_dir with the absolute VDJPuzzle directory

`export PATH=/path_to_vdjpuzzle_dir/bin:$PATH`

VDJPuzzle requires the [Ensembl reference genome](https://ccb.jhu.edu/software/tophat/igenomes.shtml). Other reference genome can be utilized (see details below). This contains the bowtie index and genome annotation required to run VDJPuzzle.

run an example with `nohup vdjpuzzle Example --bowtie-index=path_to_bt2_index/genome --gtf=path_to_gene_annotations.gtf > LOG.txt &` from the VDJPuzzle directory, you can run it on a different directory but make sure that "Example" is pointing to the Example directory in this repository.

This command will take approximaly 30 minutes to complete. You will find the output in the summary_corrected directory.

## Run VDJPuzzle with a different reference genome
VDJPuzzle uses a BED file to locate the position of the VDJ genes in the genome. The BED files provided are built for the Ensembl reference genome.
If you would like to use a different reference genome, you can generate a new BED file using this [Python script](https://bitbucket.org/kirbyvisp/marmo/src/7cfeada825fb9a00d07ebe89a7e8599550b709f1/scripts/extract_receptors.py?at=master&fileviewer=file-view-default).

## Execution and parameters

Usage: `vdjpuzzle rna_seq_directory_name --bowtie-index=path_to_bt2_index/genome_prefix --gtf=path_to_gene_annotations.gtf [option]`

Note that --bowtie-index and --gtf parameters are mandatory. 

rna_seq_directory_name contains the fastq files organized by single cell (i.e. one sub-directory for each cell that include the fastq files from that cell, check the structure of the Example directory). All fastq files need to be zipped e.g. fastq.gz and paired data needs to be specified using \_1 and \_2 in the file name.

|parameter|description|
| ------------- |-------------|
|--help|show the help|
|--qsub|executes on the cluster (overrides --CPU flag)|
|--type=(t⎮b)|specifies tcell (default) or bcell analysis|
|--CPU=n|runs with n processes|
|--THR=n|runs bowtie and tophat with n threads (default 8)|
|--species=(human⎮mouse)|specified human (default) or mouse organism|
|--no-err-corr|Do not perform final error correction on consensus sequence|
|--only-statistics|Executes only summary statistics script|
|--no-statistics|Do not execute summary statistics script|
|--transcriptomic|Enable cuffquant/cuffnorm gene quantification measurements|
|--trim|Trim reads using Trimmomatic|
|--counts|Enable featureCounts gene quantification measurements|
|--bowtie-index=path\_to\_bt2\_index|Location of the bowtie index files including the prefix (e.g. /path/to/bt2/genome)|
|--gtf=path\_to\_gtf\|Location of the GTF annotation file|

An additional script to plot gene expression as an heatmap annotated with mutation rates and other phenotype data is provided in scripts/mutation\_gene\_expression\_analysis.R

|parameter|description|
| ------------- |-------------|
|--help|show the help|
|-g file|Gene annotation in the CuffNorm directory|
|-f file|CuffNorm FPKM matrix|
|-a file|Annotation file for each cell. First column contains cell ID|

## Run VDJPuzzle on a cluster
VDJPuzzle support the execution on a system with PBS scheduler by adding the --qsub option. Every system has different parameters, thus make sure to change these parameters at the beginning of the .sh files in the script directory. 


## Citation

Simone Rizzetto, David NP Koppstein, Jerome Samir, Mandeep Singh, Joanne H Reed, Curtis H Cai, Andrew R Lloyd, Auda A Eltahla, Christopher C Goodnow, and Fabio Luciani. B-cell receptor reconstruction from single-cell RNA-seq with VDJPuzzle. Bioinformatics, 2018

## Contact

If you have any problem running VDJPuzzle please contact Fabio Luciani (luciani@unsw.edu.au)