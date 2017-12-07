#!/usr/bin/env bash

for DIR in test-data/*; do
    if [ ! -d "$DIR" ]; then
        continue
    fi

    cd $DIR
    echo "Processing $DIR"

    GENOME_REFERENCE=(reference*.fa)

    singularity exec docker://viklund/k9-bwa bwa index $GENOME_REFERENCE
    singularity exec docker://viklund/k9-bwa samtools faidx $GENOME_REFERENCE

    singularity exec docker://viklund/k9-picard picard CreateSequenceDictionary.jar \
        R=$GENOME_REFERENCE \
        O=${GENOME_REFERENCE/%.fa/.dict}

    cd ../..
done
