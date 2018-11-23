# K9-WGS-Pipeline
[![Travis status][travis-badge]][travis-link]

Nextflow pipeline for standardised mapping and variant calls on canine genomes.
The pipeline take fastqs (or bams) and outputs, bams, and gvcfs for joint
genotyping, followed by hard filtering.


## Usage

```bash
# nextflow run main.nf --help
N E X T F L O W  ~  version 0.31.1
Launching `main.nf` [shrivelled_fermi] - revision: 131a72393f
    Usage:
        nextflow run main.nf --fastqDir <directory>
    Options:
        --help
           Print this help message
        --fastqDir <Dir>
           Directory containing fastq samples (.fq.gz)
        --bamDir <Dir>
           Instead of --fastqDir, directory containing bam and indexed bam sample files (.bam, .bam.bai)
        --reference <file>
           Genome reference file (has to be indexed)
        --known <file>
           File with known sites for quality calibration.
        --outdir <dir>
           Directory for output files
        --onlyMap
           Only run the mapping steps
        --project
           Slurm project to run with
```


## Installation


The recommended way is to clone it from github:

```bash
# git clone https://github.com/NBISweden/K9-WGS-Pipeline.git
# cd K9-WGS-Pipeline
```


## Prerequisites

### Software

The pipeline requires [nextflow](https://www.nextflow.io/) and
[singularity](https://www.sylabs.io/singularity/) on the target system. These
are often pre-installed on HPC systems.

It is recommended that you pre-pull all the singularity images required by the
workflow, there is a script in the workflow directory to help you with this,
just run:

```bash
# scripts/pull_singularity.sh
```

### Data

The pipeline requries a `bwa` indexed reference genome and a `picard` genomic
dictionary. These can be created like this:

```bash
# bwa index ref.fasta
# java -Xmx4g /picard/CreateSequenceDictionary.jar REFERENCE=ref.fasta OUTPUT=ref.dict
```

You can also do this directly through the prepulled singularity images like so:

```bash
# singularity exec singularity/NBISweden-K9-WGS-Pipeline-bwa-0.7.12.img \
      bwa index ref.fasta
# singularity exec singularity/NBISweden-K9-WGS-Pipeline-picard-2.10.6.img \
      picard CreateSequenceDictionary REFERENCE=ref.fasta OUTPUT=ref.dict
```

## How to run


```
# nextflow run main.nf <options>
```

When running the mapping the workflow expects a `--fastqDir <dir>` parameter on
the command line. This directory should contain, gzipped, paired fastq files
with R1 and R2 in their filenames respectively. Specifically the filenames of
the readfiles should match either of these glob patterns.

 - `*R{1,2}*.fq.gz`
 - `*R{1,2}*.fastq.gz`

The known file is a bed file for quality calibration using BaseRecalibrator
from the GATK toolkit.

If specifying `onlyMap` no genotyping will be done.

If you already have created your mapping you can use `--bamDir` instead of
`--fastqDir` to specify a directory with bam files to run from.

### Running on HPC

## Output

The `--outdir` will have the following layout

```
<outdir>
out-tiny-fastq/
├── bam
│   ├── Sample.bai
│   └── Sample.bam
├── genotype
│   ├── all_filtered_indels_1.vcf
│   ├── all_filtered_snps_1.vcf
│   ├── all_pass_INDEL_1.recode.vcf
│   ├── all_pass_SNP_1.recode.vcf
│   ├── all_raw_INDEL_1.vcf
│   ├── all_raw_SNP_1.vcf
│   ├── all.vcf.gz
│   └── all.vcf.gz.tbi
├── haplotypeCaller
│   ├── Sample.g.vcf.gz
│   └── Sample.g.vcf.gz.tbi
└── reports
    ├── k9_dag.dot
    ├── k9_report.html
    ├── k9_timeline.html
    ├── k9_trace.txt
    ├── <Sample>.flagstat
    ├── <Sample>.marked.metrics
    ├── <Sample>.post_recal_data.table
    ├── <Sample>_R1_fastqc.html
    ├── <Sample>_R1_fastqc.zip
    ├── <Sample>_R2_fastqc.html
    ├── <Sample>_R2_fastqc.zip
    ├── <Sample>.recal_data.table
    ├── <Sample>.stats
    └── <Sample>.wgs_metrics
```

Most of this is fairly selfexplanatory, except for the reports directory.

The `k9_*` files are information from the workflow engine about how the whole
workflow went with timings and such. Then there are one set of `<Sample>*`
files for each pair of fastq files that the workflow has processed with
information on how the mapping went.


## Run on test data

First setup the testdata with `scripts/setup_testdata.sh` and then you can run
tests with the `scripts/test-one.sh` script.

```bash
# scripts/setup_testdata.sh
# scripts/test-one.sh singularity tiny fastq chr38
```

The `test-one.sh` script is mostly for testing on travis but it is very
convenient to use for local tests.

## Authors

- [Johan Viklund](https://github.com/viklund)
- [Jorrit Boekel](https://github.com/glormph)


[travis-badge]: https://api.travis-ci.org/NBISweden/K9-WGS-Pipeline.svg
[travis-link]: https://travis-ci.org/NBISweden/K9-WGS-Pipeline
