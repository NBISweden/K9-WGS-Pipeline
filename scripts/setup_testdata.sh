#!/usr/bin/env bash

singularity pull --name bwa.simg shub://NBISweden/K9-WGS-Pipeline:bwa-0.7.12 &
singularity pull --name picard.simg shub://NBISweden/K9-WGS-Pipeline:picard-1.97 &

wait

for DIR in test-data/*; do
    if [ ! -d "$DIR" ]; then
        continue
    fi

    cd $DIR
    echo "Processing $DIR"

    GENOME_REFERENCE=(reference*.fa)

    singularity exec ../../bwa.simg bwa index $GENOME_REFERENCE
    singularity exec ../../bwa.simg samtools faidx $GENOME_REFERENCE

    singularity exec ../../picard.simg picard CreateSequenceDictionary.jar \
        R=$GENOME_REFERENCE \
        O=${GENOME_REFERENCE/%.fa/.dict} || echo "Picard fail"

    cd ../..
done

rm bwa.simg
rm picard.simg
