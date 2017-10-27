#!/bin/bash

for dir in {bwa,gatk,htslib,picard,samtools}; do
    # Not more than 4 jobs at the same time
    while [ $(jobs -p | wc -l) -ge 4 ]; do
        sleep 1
    done
    cd $dir
    docker build --tag viklund/k9-${dir}:latest . &
    cd ..
done

wait
