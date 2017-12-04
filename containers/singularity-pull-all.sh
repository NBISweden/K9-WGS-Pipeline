#!/bin/bash

for dir in $(ls */Dockerfile | sed 's,/Dockerfile,,'); do
    # Not more than 4 jobs at the same time
    while [ $(jobs -p | wc -l) -ge 4 ]; do
        sleep 1
    done
    singularity pull docker://viklund/k9-${dir}:latest &
done

wait
