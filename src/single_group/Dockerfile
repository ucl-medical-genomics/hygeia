FROM rocker/r-ver:latest

# Needed by Nextflow
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("tidyverse", "argparser", "Rcpp", "RcppArmadillo"))'

WORKDIR /hygeia/single_group
COPY ./src/ ./src/
RUN Rscript ./src/r/build.R --src_root_dir ./src

COPY ./bin ./bin
ENV HYGEIA_SINGLE_GRP_BASE_PATH="/hygeia/single_group/"
ENV PATH="/hygeia/single_group/bin:${PATH}"
