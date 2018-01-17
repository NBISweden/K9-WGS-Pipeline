#!/bin/bash

for dir in $(ls */Dockerfile | sed 's,/Dockerfile,,'); do
    docker pull viklund/k9-${dir}:latest
done
