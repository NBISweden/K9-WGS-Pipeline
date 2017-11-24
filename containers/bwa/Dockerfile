FROM conda/miniconda3:latest

LABEL author="Johan Viklund" \
      maintainer="johan.viklund@nbis.se"

RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda install bwa=0.7.12 \
                     samtools=1.5

RUN mkdir /pica /proj /scratch /sw /vagrant
