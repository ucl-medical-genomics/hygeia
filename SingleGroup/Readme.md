# Single Group Analysis

Hygeia can be used for single group analyses and is useful for the following use-cases:

* TODO
* TODO
* TODO

You can find out more by reading our paper: "A Bayesian framework for genome-wide inference of DNA methylation levels", available at: https://arxiv.org/abs/2211.07311.

## Usage

### Local Installation

After cloning the repository, you can run the following tools:

### Setup

```bash
# Ensure GCC compiler is available.

Rscript -e 'install.packages(c("tidyverse", "argparser", "Rcpp", "RcppArmadillo"))'
```

#### Build the codebase

```bash
Rscript src/r/build.R
```

#### Parameter Setup

First you need to setup your parameters:

```bash
mkdir ./example/params;
Rscript specify_parameters.R \
  --mu_csv_file ./example/params/mu.csv \
  --sigma_csv_file ./example/params/sigma.csv \
  --omega_csv_file ./example/params/omega.csv \
  --kappa_csv_file ./example/params/kappa.csv \
  --u_csv_file ./example/params/u.csv \
  --p_csv_file ./example/params/p.csv
```

Note that you can change the default params through the command line arguments.

### Simulating exemplar data (optional)

It is possible to simulate data from the model in case you want to test the
algorithm without providing your own data set.

```bash
mkdir example/simulated_data;
Rscript simulate_data.R \
  --mu_csv_file ./example/params/mu.csv \
  --sigma_csv_file ./example/params/sigma.csv \
  --omega_csv_file ./example/params/omega.csv \
  --kappa_csv_file ./example/params/kappa.csv \
  --u_csv_file ./example/params/u.csv \
  --p_csv_file ./example/params/p.csv \

  --regimes_csv_file ./example/simulated_data/regimes.csv \
  --n_methylated_reads_csv_file ./example/simulated_data/n_methylated_reads.csv \
  --genomic_positions_csv_file ./example/simulated_data/genomic_positions.csv \
  --n_total_reads_csv_file ./example/simulated_data/n_total_reads.csv \

  --number_of_samples 2 \
  --number_of_cpg_sites 250
```

### Running Hygeia Analysis

The script will run the single-group algorithm to estimate regime probabilities
and potentially the unknown model parameters. This tool require three input files
which need to match the format of the exemplar files:

1. genomic_positions.csv: holds the genomic positions of the CpG sites (in ascending order) in the first column, e.g., cell (2, 1) holds the genomic position of the second CpG site.
2. n_total_reads.csv: holds the total number of reads for each CpG site (rows) and each sample (columns), e.g. the number in cell (2, 3) is the total number of reads at CpG Site 2 in Sample 3.
3. n_methylated_reads.csv: holds the number of methylated reads for each CpG site (rows) and each sample (columns), e.g. the number in cell (2, 3) is the number of methylated reads at CpG Site 2 in Sample 3

```bash
mkdir ./example2/final_output
Rscript estimate_parameters_and_regimes.R \
  --mu_csv_file ./example2/params/mu.csv \
  --sigma_csv_file ./example2/params/sigma.csv \
  --omega_csv_file ./example2/params/omega.csv \
  --kappa_csv_file ./example2/params/kappa.csv \
  --u_csv_file ./example2/params/u.csv \
  --p_csv_file ./example2/params/p.csv \
  --n_methylated_reads_csv_file ./example2/simulated_data/n_methylated_reads.csv \
  --genomic_positions_csv_file ./example2/simulated_data/genomic_positions.csv \
  --n_total_reads_csv_file ./example2/simulated_data/n_total_reads.csv \
  --regime_probabilities_csv_file ./example2/final_output/regimes.csv \
  --theta_trace_csv_file ./example2/final_output/theta_trace.csv \
  --estimate_regime_probabilities \
  --estimate_parameters
```
