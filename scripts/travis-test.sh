#!/bin/bash
set -xeuo pipefail

if [ -z "$PROFILE" ]; then
    echo "No profile specified"
    exit 1
fi

FAIL=0

for D in tiny small; do
    for T in bam fastq; do
        if ! ./scripts/test-one.sh $PROFILE $D $T; then
            echo "FAIL: $PROFILE $D $T"
            FAIL=1
        fi
    done
done

if [ $FAIL -eq 1 ]; then
    exit 1
fi
