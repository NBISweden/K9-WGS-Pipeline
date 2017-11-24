#!/bin/bash
set -xeuo pipefail

nextflow run main.nf --verbose

# Very simple check to make sure that we have some output
outbytes=$(wc -c out/file.bam | awk '{print $1}')

if [ $outbytes -lt 100 ]; then
    exit 1
fi
