#!/bin/bash
set -xeuo pipefail

if [ -z "$PROFILE" ]; then
    echo "No profile specified"
    exit 1
fi

nextflow run main.nf -profile $PROFILE --verbose --fastqDir test-data/test-data-tiny

# Very simple check to make sure that we have some output
outbytes=$(wc -c out/file.bam | awk '{print $1}')

if [ $outbytes -lt 100 ]; then
    exit 1
fi
