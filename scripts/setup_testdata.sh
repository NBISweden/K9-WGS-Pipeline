#!/usr/bin/env bash

singularity pull --name bwa.simg shub://NBISweden/K9-WGS-Pipeline:bwa-0.7.12 &
singularity pull --name picard.simg shub://NBISweden/K9-WGS-Pipeline:picard-1.97 &

wait

shopt -s nullglob

for dirpath in test-data/*/; do
    cd "$dirpath" || exit 1
    printf 'Processing "%s"\n' "$dirpath"

    for genome_reference in reference*.fa; do
        singularity exec ../../bwa.simg bwa index "$genome_reference"
        singularity exec ../../bwa.simg samtools faidx "$genome_reference"

        if ! singularity exec ../../picard.simg picard CreateSequenceDictionary.jar \
            R="$genome_reference" \
            O="${genome_reference%.fa}.dict"
        then
             printf 'Picard fail for "%s"\n' "$genome_reference" >&2
        fi
    done

    for bamfile in *.bam; do
        singularity exec ../../bwa.simg samtools index "$bamfile"
    done

    cd "$OLDPWD"
done

rm -f bwa.simg
rm -f picard.simg
