#!/usr/bin/env bash

dirpath=singularity
if ! mkdir -p "$dirpath"; then
    printf 'Could not create directory "%s"\n' "$dirpath" >&2
    exit 1
fi

printf 'Will store singularity images in "%s"\n' "$dirpath"

awk '/shub:/ { gsub("\"", ""); print $NF }' conf/singularity.config |
while read -r image; do
    store=${image#shub://}
    store=${store//[:\/]/-}.img
    if [ ! -f "$dirpath/$store" ]; then
        printf '%s\n' "$store"
        singularity pull --name "$dirpath/$store" "$image"
    fi
done
