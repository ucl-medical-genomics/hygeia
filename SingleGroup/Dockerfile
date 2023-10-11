FROM r-base:4.3.1

RUN Rscript -e 'install.packages(c("tidyverse", "argparser", "Rcpp", "RcppArmadillo"))'

COPY . /hygeia/SingleGroup

RUN Rscript /hygeia/SingleGroup/src/r/build.R --src_root_dir /hygeia/SingleGroup/src

