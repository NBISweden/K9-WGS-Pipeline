#!/bin/bash

# Install Nextflow
cd $HOME
curl -fsSL https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Install Singularity, if needed
if [ ! -z "$SGT_VER" ]; then
    cd $HOME
    wget https://github.com/singularityware/singularity/releases/download/$SGT_VER/singularity-$SGT_VER.tar.gz
    tar xvf singularity-$SGT_VER.tar.gz
    cd singularity-$SGT_VER
    ./configure --prefix=/usr/local
    make
    sudo make install
    cd ..
    rm -rf singularity-$SGT_VER*
fi
