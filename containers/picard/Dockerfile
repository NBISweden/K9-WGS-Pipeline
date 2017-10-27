FROM conda/miniconda3:latest

LABEL author="Johan Viklund" \
      maintainer="johan.viklund@nbis.se"

RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda install picard=1.97

RUN mkdir /pica /proj /scratch /sw
