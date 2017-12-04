#!/bin/bash

## Needed for all, no need to download one thing in parallel
docker pull conda/miniconda3:latest

for dir in $(ls */Dockerfile | sed 's,/Dockerfile,,'); do
    # Not more than 4 jobs at the same time
    while [ $(jobs -p | wc -l) -ge 4 ]; do
        sleep 1
    done
    cd $dir
    docker build --tag viklund/k9-${dir}:latest . &
    cd ..
done

wait
