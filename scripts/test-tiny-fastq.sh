#!/bin/bash

nextflow run main.nf \
    -profile singularity \
    --verbose \
    --fastqDir  test-data/test-data-tiny/ \
    --reference test-data/test-data-tiny/reference.fa \
    --known     test-data/test-data-tiny/known.bed \
    --out       out-tiny-fastq
