#!/bin/bash

echo "INSTALL DEPS"
sudo apt-get install squashfs-tools
./scripts/travis-install.sh

echo "SETUP TEST DATA"
./scripts/setup_testdata.sh

echo "RUN TRAVIS TEST"
./scripts/travis-test.sh

