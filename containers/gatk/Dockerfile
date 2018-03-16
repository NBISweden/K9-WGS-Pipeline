FROM broadinstitute/gatk3:3.5-0

LABEL author="Johan Viklund" \
      maintainer="johan.viklund@nbis.se"

RUN apt-get -y install pkg-config && \
    curl -LO https://github.com/vcftools/vcftools/releases/download/v0.1.15/vcftools-0.1.15.tar.gz && \
    tar xzf vcftools-0.1.15.tar.gz && \
    cd vcftools-0.1.15 && \
    ./configure && \
    make && \
    make install

RUN mkdir /pica /proj /scratch /sw
