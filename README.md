# hic-scaffolding
This repo contains the scripts that can be used for Hi-C scaffolding with Arima Genomics Hi-C libraries 

## Installation

### This installation guide will provide instructions to install this pipeline on a Compute Canada (Aliance) cluster

First you will need to clone the repo
```
git clone https://github.com/pubudumanoj/hic-scaffolding.git
cd hic-scaffolding
```

Next, either download and install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html) or load the NextFlow module installed in the CC file system.
_e.g_ 
```
module load StdEnv/2020 nextflow/21.04.3
```
Then specify the directory path of the Hi-C fastq files in **in_dir** param. Make sure to add a "/" at the end of the path.

#### How to name the fastq files

Currently this pipeline only supports for paired end fastq files and you must have both forward and reverse reads in order to use the pipeline. If you have fastq files from multiple samples you can name samples as follows

e.g
sample 1
HiC05-01_Afraterculus_RAG1429A1-1_S1_**L001**_**R1**_001.fastq.gz
HiC05-01_Afraterculus_RAG1429A1-1_S1_**L001**_**R2**_001.fastq.gz

sample 2
HiC05-01_Afraterculus_RAG1429A1-1_S1_**L002**_**R1**_001.fastq.gz
HiC05-01_Afraterculus_RAG1429A1-1_S1_**L002**_**R1**_001.fastq.gz

There should be a common part for all the names of the samples and sample can be uniquily identified by a sample ID (L002 in above example). This should followed by the read type and the rest should be similar

