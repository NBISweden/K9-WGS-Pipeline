#!/usr/bin/env bash

# NBISweden-K9-WGS-Pipeline-bwa-0.7.12.img
for IMG in $( cat conf/singularity.config | grep shub | awk '{print $NF}' ); do
    IMG=$( echo $IMG | sed 's/"//g')
    STORE=$( echo $IMG | sed 's.shub://..;s./.-.g;s.:.-.g;s/$/.img/' )
    echo $STORE
    singularity pull --name "$STORE" "$IMG"
done
