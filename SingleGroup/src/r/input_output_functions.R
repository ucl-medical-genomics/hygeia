suppressMessages(library("tidyverse"))

## Auxiliary mathematical functions
## ========================================================================== ##

## Auxiliary functions for writing simulated data and algorithm output to CSV
## files
## ========================================================================== ##

# Writes the genomic_positions to a csv file.
write_genomic_positions_to_csv_file <- function(
  genomic_positions, # numeric (integer) vector of genomic positions in ascending order
  file
) {
  tibble(genomic_positions = genomic_positions) %>%
    write_csv(file = file.path(file))
}
# Writes the total number of reads for each sample/CpG site to a csv file.
write_n_total_reads_to_csv_file <- function(
  n_total_reads, # (n_samples, n_cpg_sites)-numeric (integer) matrix containing the total number of reads for each sample/CpG site
  file
) {
  n_samples <- nrow(n_total_reads)
  n_total_reads %>%
    t() %>%
    as_tibble(.name_repair = ~ paste0("sample_", 1:n_samples)) %>%
    write_csv(file = file.path(file))
}
# Writes the number of methylated reads for each sample/CpG site to a csv file.
write_n_methylated_reads_to_csv_file <- function(
  n_methylated_reads, # (n_samples, n_cpg_sites)-numeric (integer) matrix containing the number of methylated reads for each sample/CpG site
  file
) {
  n_samples <- nrow(n_methylated_reads)
  n_methylated_reads %>%
    t() %>%
    as_tibble(.name_repair = ~ paste0("sample_", 1:n_samples)) %>%
    write_csv(file = file.path(file))
}
# Writes the (simulated) regimes to a csv file.
write_regimes_to_csv_file <- function(
  regimes, # numeric (integer) vector of regimes at each CpG site
  file
) {
  tibble(regime = regimes) %>%
    write_csv(file = file.path(file))
}
# Writes the regime-transition matrix $P$ to a CSV file.
write_p_to_csv_file <- function(
  p, # an (n_methylation_regimes, n_methylation_regimes)-stochastic matrix whose element $(i, j)$ is the probability of switching from Regime $i$ to Regime $j$ at a change point. Must have "0" on the diagonal.
  file
) {
  p %>%
    as_tibble(.name_repair = ~ paste0("regime_", 1:nrow(p))) %>%
    write_csv(file = file.path(file))
}
# Writes the parameter vector $\mu$ to a CSV file.
write_mu_to_csv_file <- function(
  mu, # a vector of length n_methylation_regimes
  file
) {
  tibble(mu = mu) %>% write_csv(file.path(file))
}
# Writes the parameter vector $\sigma$ to a CSV file.
write_sigma_to_csv_file <- function(
  sigma, # a vector of length n_methylation_regimes
  file
) {
  tibble(sigma = sigma) %>% write_csv(file.path(file))
}
# Writes the parameter vector $\kappa$ to a CSV file.
write_kappa_to_csv_file <- function(
  kappa, # a vector of length n_methylation_regimes
  file
) {
  tibble(kappa = kappa) %>% write_csv(file.path(file))
}
# Writes the parameter vector $\omega$ to a CSV file.
write_omega_to_csv_file <- function(
  omega, # a vector of length n_methylation_regimes
  file
) {
  tibble(omega = omega) %>% write_csv(file.path(file))
}
# Writes the parameter $u$ to a CSV file.
write_u_to_csv_file <- function(
  u, # an integer > 1
  file
) {
  tibble(u = u) %>% write_csv(file.path(file))
}

## Auxiliary functions for reading input data and some of the
# vector- or matrix-valued input parameters from CSV files
## ========================================================================== ##

# Reads the genomic positions from a CSV file and converts to a vector,
# ready to pass to the algorithm.
read_genomic_positions_from_csv_file <- function(
  file # CSV file with a single column (and n_cpg_sites rows) containing the genomic positions in ascending order
) {
  read_csv(file = file, show_col_types = FALSE) %>% pull(1)
}
# Reads the total number of reads from a CSV file and converts to a matrix,
# ready to pass to the algorithm.
read_n_total_reads_from_csv_file <- function(
  file # CSV file with n_samples columns and n_cpg_sites rows containing the total number of reads for each sample/CpG site
) {
  read_csv(file = file, show_col_types = FALSE) %>% as.matrix() %>% t()
}
# Reads the total number of reads from a CSV file and converts to a matrix,
# ready to pass to the algorithm.
read_n_methylated_reads_from_csv_file <- function(
  file # CSV file with n_samples columns and n_cpg_sites rows containing the number of methylated reads for each sample/CpG site
) {
  read_csv(file = file, show_col_types = FALSE) %>% as.matrix() %>% t()
}
# Reads the parameter vector $\mu$ from a CSV file.
read_mu_from_csv_file <- function(
  file # must contain a file "mu.csv" with n_methylation_regimes entries in the first column.
) {
  read_csv(file = file.path(file), show_col_types = FALSE) %>% pull(1)
}
# Reads the parameter vector $\sigma$ from a CSV file.
read_sigma_from_csv_file <- function(
  file # must contain a file "sigma.csv" with n_methylation_regimes entries in the first column.
) {
  read_csv(file = file.path(file), show_col_types = FALSE) %>% pull(1)
}
# Reads the parameter vector $\kappa$ from a CSV file.
read_kappa_from_csv_file <- function(
  file # must contain a file "kappa.csv" with n_methylation_regimes entries in the first column.
) {
  read_csv(file = file.path(file), show_col_types = FALSE) %>% pull(1)
}
# Reads the parameter vector $\omega$ from a CSV file.
read_omega_from_csv_file <- function(
  file # must contain a file "omega.csv" with n_methylation_regimes entries in the first column.
) {
  read_csv(file = file.path(file), show_col_types = FALSE) %>% pull(1)
}
# Reads the regime-transition matrix $P$ from a CSV file.
read_p_from_csv_file <- function(
  file # must contain a file "p.csv" with n_methylation_regimes rows and columns.
) {
  read_csv(file = file.path(file), show_col_types = FALSE) %>%
    as.matrix() %>%
    `dimnames<-`(NULL)
}
# Reads the parameter $u$ from a CSV file.
read_u_from_csv_file <- function(
  file # must contain a file "u.csv" with a single integer-valued cell in the first row/column
) {
  read_csv(file = file.path(file), show_col_types = FALSE) %>% pull(1)
}

