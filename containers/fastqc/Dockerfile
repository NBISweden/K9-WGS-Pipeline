FROM openjdk:8

# This is copied from https://github.com/SciLifeLab/CAW
LABEL \
  author="Johan Viklund" \
  description="FastQC image" \
  maintainer="johan.viklund@gmail.com"

# Install libraries
RUN \
  apt-get update && apt-get install -y --no-install-recommends \
    wget \
  && rm -rf /var/lib/apt/lists/*

# Setup ENV variables
ENV FASTQC_VERSION=0.11.5

# Install FastQC
RUN \
  wget --quiet -O fastqc_v${FASTQC_VERSION}.zip \
    http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip \
  && unzip fastqc_v${FASTQC_VERSION}.zip -d /opt/ \
  && chmod 755 /opt/FastQC/fastqc \
  && ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc \
  && rm fastqc_v${FASTQC_VERSION}.zip

# Create UPPMAX directories
RUN mkdir /pica /proj /scratch /sw
