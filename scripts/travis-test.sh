#!/bin/bash
set -xeuo pipefail

if [ -z "$PROFILE" ]; then
    echo "No profile specified"
    exit 1
fi

nextflow run main.nf \
    -profile singularity \
    --verbose \
    --fastqDir  test-data/test-data-tiny \
    --reference test-data/test-data-tiny/reference_chr38-1000000-1030000.fa \
    --known     test-data/test-data-tiny/known_chr38-1000000-1030000.bed \
    --out       out-tiny

# Very simple check to make sure that we have some output
outbytes=$(zcat out-tiny/all_samples_genotyping.vcf.gz | grep -v '^#' | wc -c)

if [ $outbytes -lt 1000 ]; then
    exit 1
fi

nextflow run main.nf \
    -profile singularity \
    --verbose \
    --fastqDir  test-data/test-data-small/ \
    --reference test-data/test-data-small/reference_subset.fa \
    --known     test-data/test-data-small/known_subset.bed \
    --out       out-small

# Very simple check to make sure that we have some output
outbytes=$(zcat out-small/all_samples_genotyping.vcf.gz | grep -v '^#' | wc -c)

if [ $outbytes -lt 1000 ]; then
    exit 1
fi
