#!/bin/bash

# Install Nextflow
cd $HOME
curl -fsSL https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

singularity_url() {
    if [ "$SINGULARITY_VER" == "3.0.1" ]; then
        echo "https://github.com/sylabs/singularity/releases/download/v3.0.1/singularity-3.0.1.tar.gz"
    elif [ "$SINGULARITY_VER" == "3.0.0" ]; then
        echo "https://github.com/sylabs/singularity/releases/download/v3.0.0/singularity-v3.0.0.tar.gz"
    else
        echo "https://github.com/singularityware/singularity/releases/download/$SGT_VER/singularity-$SGT_VER.tar.gz"
    fi
}

# Install Singularity, if needed
if [ ! -z "$SGT_VER" ]; then
    cd $HOME
    URL=singularity_url
    wget "$URL"
    tar xvf singularity-*$SGT_VER.tar.gz # * since sometimes there's a v, sometimes not...
    cd singularity-$SGT_VER
    ./configure --prefix=/usr/local
    make
    sudo make install
    cd ..
    rm -rf singularity-$SGT_VER*
fi
