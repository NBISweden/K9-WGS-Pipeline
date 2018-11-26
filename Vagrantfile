# -*- mode: ruby -*-
# vi: set ft=ruby :

Vagrant.configure("2") do |config|
  config.vm.box = "ubuntu/xenial64"

  config.vm.provider "virtualbox" do |vb|
    vb.memory = "8128"
    vb.cpus = 4
  end

  config.vm.provision "shell", inline: <<-SHELL
    # Dependencies
    apt-get update
    apt-get -y install build-essential curl git sudo man vim autoconf libtool \
        python-minimal python3 openjdk-8-jre linux-image-extra-$(uname -r) \
        linux-image-extra-virtual apt-transport-https ca-certificates \
        software-properties-common libssl-dev uuid-dev libgpgme11-dev \
        squashfs-tools libseccomp-dev pkg-config

    # Go
    export VERSION=1.11 OS=linux ARCH=amd64
    cd /tmp
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz >/dev/null 2>&1
    sudo tar -C /usr/local -xzf go$VERSION.$OS-$ARCH.tar.gz
    export GOPATH=${HOME}/go
    export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin

    # Singularity
    mkdir -p $GOPATH/src/github.com/sylabs
    cd $GOPATH/src/github.com/sylabs
    git clone https://github.com/sylabs/singularity.git >/dev/null 2>&1
	cd singularity
	./mconfig
	cd ./builddir
	make
	sudo make install


    # Nextflow
    curl -s https://get.nextflow.io 2>/dev/null | bash
    mkdir -p /home/ubuntu/bin
    mv nextflow /home/ubuntu/bin
    chown -R ubuntu:ubuntu /home/ubuntu/bin

    su ubuntu -c /home/ubuntu/bin/nextflow # Downloads all dependencies

    # Docker
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg 2>/dev/null | sudo apt-key add -
    sudo add-apt-repository \
        "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
        $(lsb_release -cs) \
        stable"
    apt-get update
    apt-get -y install docker-ce

    usermod -a -G docker ubuntu
  SHELL
end
