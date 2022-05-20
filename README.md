[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)

# hic-scaffolding
hid-scaffolding is a pipeline that can be used for `Hi-C scaffolding` with `Arima Genomics` Hi-C libraries. The pipeline uses [SALSA2](https://github.com/marbl/SALSA) for the scaffolding and it is currently designed for [Compute Canada Architecture](https://status.computecanada.ca/)

## Installation

### This installation guide will provide instructions to install this pipeline on a Compute Canada (Aliance) cluster

First you will need to clone the repo
```
git clone https://github.com/pubudumanoj/hic-scaffolding.git
cd hic-scaffolding
```

Next, either download and install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html) `>=21.10.3` or load the NextFlow module installed in the CC file system.
_e.g_ 
```
module load StdEnv/2020 nextflow/21.04.3
```

## Usage

After installing, specify the directory path of the Hi-C fastq files that you want to use for the scaffolding process in the `in_dir` param. Make sure to add a "/" at the end of the path.

#### How to name the fastq files

Currently this pipeline only supports for paired end fastq files and you must have both forward and reverse reads in order to use the pipeline. If you have fastq files from multiple samples you can name samples as follows

_e.g_ <br />
**sample 1** <br />
`HiC_Afraterculus_L001_R1_001.fastq.gz`  <br />
`HiC_Afraterculus_L001_R2_001.fastq.gz`

**sample 2** <br />
`HiC_Afraterculus_L002_R1_001.fastq.gz`  <br />
`HiC_Afraterculus_L002_R2_001.fastq.gz`

There should be a common part for all the names of the samples and sample can be uniquily identified by a sample ID (`L001 and L002` in above example). This should followed by the read type `(R1 and R2)` and the rest should be similar.

After correctly formatting fastq file names you should change the `fastq` param accordingly to match the below REGEX pattern <br />
`'*R{1,2}_001.fastq.gz'`

Then you should specify the path for the contigs assembly (reference fasta file) using `REF` param in the config file. Also make sure to clean the scaffold names in the fasta file. If you use simpler form the output files will be smaller in size and easy to process.
e.g
`>scaffold_1`

Optionally, you can modifiy each parameter defined in the config file accordingly.

To run the pipeline use this code.

```
nextflow run main2.nf -resume
```


