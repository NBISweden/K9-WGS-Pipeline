#!/bin/bash
set -xeuo pipefail

function chrom_for() {
    if [ "$1" == "tiny" ]; then
        echo "chr38"
    elif [ "$1" == "small" ]; then
        echo "chr36,chr37,chr38,chrX"
    fi
}

if [ -z "$PROFILE" ]; then
    echo "No profile specified"
    exit 1
fi

FAIL=0

for D in tiny small; do
    C=$(chrom_for $D)
    for T in fastq bam; do
        if ! ./scripts/test-one.sh $PROFILE $D $T $C; then
            echo "FAIL: $PROFILE $D $T"
            FAIL=1
        fi
    done
done

if [ $FAIL -eq 1 ]; then
    exit 1
fi
