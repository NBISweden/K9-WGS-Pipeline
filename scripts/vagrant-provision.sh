#!/usr/bin/env bash

# Dependencies
apt-get update
apt-get -y install build-essential curl git sudo man vim autoconf libtool \
    python-minimal python3 openjdk-8-jre linux-image-extra-$(uname -r) \
    linux-image-extra-virtual apt-transport-https ca-certificates \
    software-properties-common

# Singularity
git clone https://github.com/singularityware/singularity.git
cd singularity
./autogen.sh
./configure --prefix=/usr/local
make
make install
cd ..
echo 'bind path = /scratch' >> /usr/local/etc/singularity/singularity.conf

# Nextflow
curl -s https://get.nextflow.io | bash
mkdir -p /home/ubuntu/bin
mv nextflow /home/ubuntu/bin
chown -R ubuntu:ubuntu /home/ubuntu/bin

su ubuntu -c /home/ubuntu/bin/nextflow # Downloads all dependencies

# Docker
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository \
    "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
    $(lsb_release -cs) \
    stable"
apt-get update
apt-get -y install docker-ce

