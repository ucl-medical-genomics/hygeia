# Local Development

## Requirements

The following dependencies are required

* GCC compiler
* R 4.3.1
* Python

```bash
Rscript -e 'install.packages(c("tidyverse", "argparser", "Rcpp", "RcppArmadillo"))'

Rscript src/single_group/src/r/build.R

```