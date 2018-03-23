#!/bin/bash

nextflow run main.nf \
    -profile singularity \
    --verbose \
    --bamDir    test-data/test-data-small/ \
    --reference test-data/test-data-small/reference.fa \
    --known     test-data/test-data-small/known.bed \
    --out       out-small-bam
