#!/bin/bash

if [ ! -z "$FAKE" ]; then
    echo "FAKE!!"
    exit 1
fi

echo "INSTALL DEPS"
sudo apt-get install squashfs-tools
if ! ./scripts/travis-install.sh; then
    echo "Failed in install"
    exit 1;
fi

echo "SETUP TEST DATA"
if ! ./scripts/setup_testdata.sh; then
    echo "Failed to setup test data"
    exit 1;
fi

echo "RUN TRAVIS TEST"
if ! ./scripts/travis-test.sh; then
    echo "Test failure"
    exit 1
fi

