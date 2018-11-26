#!/bin/bash

# Install Singularity, if needed
./scripts/travis-install-singularity.sh

# Install Nextflow
cd $HOME
curl -fsSL https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
