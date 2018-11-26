#!/bin/bash

# Install Singularity, if needed
if [ ! -z "$SGT_VER" ]; then
    #sudo apt-get update
    #sudo apt-get -y install build-essential curl git sudo man vim autoconf libtool \
    #    python-minimal python3 openjdk-8-jre linux-image-extra-$(uname -r) \
    #    linux-image-extra-virtual apt-transport-https ca-certificates \
    #    software-properties-common libssl-dev uuid-dev  \
    #    pkg-config # squashfs-tools libseccomp-dev libgpgme11-dev
    sudo apt-get update && sudo apt-get install -y wget git \
                                                    build-essential \
                                                    libtool \
                                                    autotools-dev \
                                                    libarchive-dev \
                                                    automake \
                                                    autoconf \
                                                    uuid-dev \
                                                    libssl-dev


    echo "Installing GO"
    export VERSION=1.11 OS=linux ARCH=amd64
    cd /tmp
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz >/dev/null 2>&1
    sudo tar -C /usr/local -xzf go$VERSION.$OS-$ARCH.tar.gz
    export GOPATH=${HOME}/go
    export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin
    echo "FULL PATH VARIABLE:"
    echo $PATH | tr ':' "\n"

    echo "Checking go version!"
    go version

    mkdir -p $GOPATH/src/github.com/sylabs
    cd $GOPATH/src/github.com/sylabs
    #git clone https://github.com/sylabs/singularity.git >/dev/null 2>&1
    git clone https://github.com/sylabs/singularity.git
    cd singularity
    git checkout "$SGT_VER"

    ./mconfig
    cd ./builddir
    echo "GREPPING GREPPING"
    grep "CC" Makefile
    echo "GREPPING GREPPING"
    grep "CFLAGS" Makefile
    cd ../
    echo "GREPPING GREPPING"
    find . -type f -exec grep -nH -- '-std' {} \;
    cd builddir
    make
    sudo make install
fi

# Install Nextflow
cd $HOME
curl -fsSL https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
