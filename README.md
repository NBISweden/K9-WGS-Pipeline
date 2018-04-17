# K9-WGS-Pipeline
[![Travis status][travis-badge]][travis-link]

Nextflow pipeline for standardised mapping and variant calls on canine genomes

## Usage

```bash
    Usage:
        nextflow run main.nf --fastqDir <directory>
    Options:
        --fastqDir <Dir>
           Directory containing fastq samples (.fq.gz)
        --bamDir <Dir>
           Instead of --fastqDir, directory containing bam and indexed bam sample files (.bam, .bam.bai)
        --reference <file>
           Genome reference file (has to be indexed)
        --known <file>
           Bed file with known sites for quality calibration.
        --outdir <dir>
           Directory for output files
        --onlyMap
           Only run the mapping steps
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

```bash
$ scripts/setup_testdata.sh
$ scripts/test-one.sh singularity tiny fastq chr38
```

The `test-one.sh` script is mostly for testing on travis but it is very
convenient to use for local tests.

## Authors

- [Johan Viklund](https://github.com/viklund)
- [Jorrit Boekel](https://github.com/glormph)


[travis-badge]: https://api.travis-ci.org/NBISweden/K9-WGS-Pipeline.svg
[travis-link]: https://travis-ci.org/NBISweden/K9-WGS-Pipeline
