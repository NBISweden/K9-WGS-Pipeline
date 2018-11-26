#!/bin/bash

# Install Nextflow
cd $HOME
curl -fsSL https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Install Singularity, if needed
if [ ! -z "$SGT_VER" ]; then
    sudo apt-get -y install build-essential curl git sudo man vim autoconf libtool \
        python-minimal python3 openjdk-8-jre linux-image-extra-$(uname -r) \
        linux-image-extra-virtual apt-transport-https ca-certificates \
        software-properties-common libssl-dev uuid-dev libgpgme11-dev \
        squashfs-tools libseccomp-dev pkg-config

    export VERSION=1.11 OS=linux ARCH=amd64
    cd /tmp
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz >/dev/null 2>&1
    sudo tar -C /usr/local -xzf go$VERSION.$OS-$ARCH.tar.gz
    export GOPATH=${HOME}/go
    export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin

    mkdir -p $GOPATH/src/github.com/sylabs
    cd $GOPATH/src/github.com/sylabs
    #git clone https://github.com/sylabs/singularity.git >/dev/null 2>&1
    git clone https://github.com/sylabs/singularity.git
	cd singularity
    git checkout "$SGT_VER"

	./mconfig
	cd ./builddir
	make
	sudo make install
fi
