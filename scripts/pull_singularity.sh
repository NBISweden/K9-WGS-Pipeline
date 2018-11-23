#!/bin/sh

DIR=singularity
mkdir $DIR

for IMG in $( cat conf/singularity.config | grep shub | awk '{print $NF}' ); do
    IMG=$( echo $IMG | sed 's/"//g')
    STORE=$( echo $IMG | sed 's.shub://..;s./.-.g;s.:.-.g;s/$/.img/' )
    if [ -f "$DIR/$STORE" ]; then
        continue
    fi
    echo $STORE
    singularity pull --name "$DIR/$STORE" "$IMG"
done
