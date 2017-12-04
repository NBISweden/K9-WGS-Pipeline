#!/bin/bash

for dir in $(ls */Dockerfile | sed 's,/Dockerfile,,'); do
    docker push viklund/k9-${dir}:latest
done
