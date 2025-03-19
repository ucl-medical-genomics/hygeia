suppressMessages(library("tidyverse"))

## Generalized Function for Writing to CSV
write_to_csv_file <- function(data, file) {
  tibble(data) %>%
    write_csv(file = file.path(file))
}

## Generalized Function for Reading from CSV
read_from_csv_file <- function(file, as_matrix = FALSE, transpose = FALSE) {
  message(paste0("Reading data from file: ", file))
  data <- read_csv(file = file, show_col_types = FALSE)
  if (as_matrix) {
    data <- as.matrix(data)
    if (transpose) {
      data <- t(data)
    }
  }
  return(data)
}

## Utility Function for Creating Directories
create_dirs_for_file <- function(file_path) {
  if (!is.character(file_path)) return()
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

## Specific Functions Using the Generalized Functions

# Writing Functions
write_genomic_positions_to_csv_file <- function(genomic_positions, file) {
  write_to_csv_file(tibble(genomic_positions = genomic_positions), file)
}

write_n_total_reads_to_csv_file <- function(n_total_reads, file) {
  n_samples <- nrow(n_total_reads)
  data <- as_tibble(t(n_total_reads), .name_repair = ~ paste0("sample_", 1:n_samples))
  write_to_csv_file(data, file)
}

write_n_methylated_reads_to_csv_file <- function(n_methylated_reads, file) {
  n_samples <- nrow(n_methylated_reads)
  data <- as_tibble(t(n_methylated_reads), .name_repair = ~ paste0("sample_", 1:n_samples))
  write_to_csv_file(data, file)
}

write_regimes_to_csv_file <- function(regimes, file) {
  write_to_csv_file(tibble(regime = regimes), file)
}

write_p_to_csv_file <- function(p, file) {
  data <- as_tibble(p, .name_repair = ~ paste0("regime_", 1:nrow(p)))
  write_to_csv_file(data, file)
}

write_mu_to_csv_file <- function(mu, file) {
  write_to_csv_file(tibble(mu = mu), file)
}

write_sigma_to_csv_file <- function(sigma, file) {
  write_to_csv_file(tibble(sigma = sigma), file)
}

write_kappa_to_csv_file <- function(kappa, file) {
  write_to_csv_file(tibble(kappa = kappa), file)
}

write_omega_to_csv_file <- function(omega, file) {
  write_to_csv_file(tibble(omega = omega), file)
}

write_u_to_csv_file <- function(u, file) {
  write_to_csv_file(tibble(u = u), file)
}

# Reading Functions
read_genomic_positions_from_csv_file <- function(file) {
  read_from_csv_file(file) %>% pull(1)
}

read_n_total_reads_from_csv_file <- function(file) {
  read_from_csv_file(file, as_matrix = TRUE, transpose = TRUE)
}

read_n_methylated_reads_from_csv_file <- function(file) {
  read_from_csv_file(file, as_matrix = TRUE, transpose = TRUE)
}

read_mu_from_csv_file <- function(file) {
  read_from_csv_file(file) %>% pull(1)
}

read_sigma_from_csv_file <- function(file) {
  read_from_csv_file(file) %>% pull(1)
}

read_kappa_from_csv_file <- function(file) {
  read_from_csv_file(file) %>% pull(1)
}

read_omega_from_csv_file <- function(file) {
  read_from_csv_file(file) %>% pull(1)
}

read_p_from_csv_file <- function(file) {
  read_from_csv_file(file, as_matrix = TRUE)
}

read_u_from_csv_file <- function(file) {
  read_from_csv_file(file) %>% pull(1)
}
