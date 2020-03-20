#!/bin/bash
set -xeuo pipefail

function chrom_for() {
    case $1 in
	tiny)  echo chr38 ;;
	small) echo chr36,chr37,chr38,chrX
    esac
}

if [ -z "$PROFILE" ]; then
    echo 'No profile specified' >&2
    exit 1
fi

FAIL=0

for D in tiny small; do
    C=$(chrom_for $D)
    for T in fastq bam; do
        if ! ./scripts/test-one.sh "$PROFILE" "$D" "$T" "$C"; then
            printf 'FAIL: %s %s %s\n' "$PROFILE" "$D" "$T" >&2
            FAIL=1
        fi
    done
done

exit "$FAIL"
