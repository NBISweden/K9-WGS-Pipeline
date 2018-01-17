#!/bin/bash

for dir in $(ls */Singularity.* | sed 's,.*/Singularity.,,'); do
    singularity pull shub://NBISweden/K9-WGS-Pipeline:${dir}
done
