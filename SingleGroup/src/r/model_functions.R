suppressMessages(library("tidyverse"))

## Auxiliary mathematical functions
## ========================================================================== ##

# The logit, i.e. inverse of the logistic function.
logit <- function(x) {
  return(log(x) - log(1.0 - x))
}
# The inverse-logit, i.e. logistic function.
inverse_logit <- function(x) {
  return(1.0 / (1.0 + exp(-x)))
}
# Normalises weights in log-space.
normalise_exp <- function(log_w) {
  log_w_max <- max(log_w)
  log_z <- log_w_max + log(sum(exp(log_w - log_w_max)))
  return(log_w - log_z)
}

## General model parameters and related functions
## ========================================================================== ##

# Returns the number of regimes.
get_n_methylation_regimes <- function(vartheta) {
  return(vartheta[2])
}
# Returns whether or not kappa is fixed.
get_is_kappa_fixed <- function(vartheta) {
  n_methylation_regimes <- vartheta[2]
  is_kappa_fixed <- vartheta[2 * n_methylation_regimes + 3]
  return(is_kappa_fixed)
}
# Returns a number of model parameters which are known/fixed, i.e. which are
# not estimated as part of the online parameter estimation algorithm.
get_known_parameters <- function(
    mu,
    sigma,
    is_kappa_fixed = TRUE, # are the parameters kappa fixed/known (otherwise, they are estimated alongside the other unknown parameters)
    kappa = rep(2, times = length(mu)),
    u # minimum distance between change points
)
{
  n_methylation_regimes <- length(mu)
  nu  <- mu * (1 - mu) / sigma^2 - 1
  alpha <- mu * nu
  beta <- (1 - mu) * nu
  if (is_kappa_fixed) {
    vartheta <- c(u, n_methylation_regimes, alpha, beta, is_kappa_fixed, kappa)
    dim_theta <- n_methylation_regimes^2
  } else {
    vartheta <- c(u, n_methylation_regimes, alpha, beta, is_kappa_fixed)
    dim_theta <- n_methylation_regimes * (n_methylation_regimes + 1)
  }

  return(
    list(
      vartheta = vartheta,
      dim_theta = dim_theta,
      n_methylation_regimes = n_methylation_regimes
    )
  )
}
# Converts the (potentially "unknown") model parameters to the vector $\theta$.
convert_model_parameters_to_theta <- function(
    vartheta,
    p,
    omega,
    kappa = rep(2, times = length(omega))
) {
  is_kappa_fixed <- get_is_kappa_fixed(vartheta)
  diag(p) <- -1
  theta <- c(log(p[p != -1]), logit(omega))
  if (!is_kappa_fixed) {
    theta <- c(theta, log(kappa))
  }
  return(theta)
}
# Converts the vector $\theta$ to the (potentially "unknown")
# model parameters.
convert_theta_to_model_parameters <- function(vartheta, theta) {

  n_methylation_regimes <- get_n_methylation_regimes(vartheta)
  RR <- n_methylation_regimes
  is_kappa_fixed <- get_is_kappa_fixed(vartheta)

  p_non_diag <- rep(NA, times = RR * (RR - 1))
  p <- matrix(0, nrow = RR, ncol = RR)
  for (rr in 1:RR) {
    idx_min <- (rr - 1) * (RR - 1) + 1
    idx_max <- rr * (RR - 1)
    p_non_diag[idx_min:idx_max] <- exp(normalise_exp(theta[idx_min:idx_max]))
    p[rr, -rr] <- p_non_diag[idx_min:idx_max]
  }

  idx_min <- RR * (RR - 1) + 1
  idx_max <- RR * RR
  omega <- inverse_logit(theta[idx_min:idx_max])

  if (!is_kappa_fixed) {
    idx_min <- RR * RR + 1
    idx_max <- RR * (RR + 1)
    kappa <- exp(theta[idx_min:idx_max])
  }

  if (is_kappa_fixed) {
    return(list(p = p, p_non_diag = p_non_diag, omega = omega))
  } else {
    return(list(p = p, p_non_diag = p_non_diag, omega = omega, kappa = kappa))
  }
}
# Simulates data (the total number of reads, the number of methylated reads,
# the change points and the associated regimes) from the model for a given
# number of CpG sites and given number of samples; writes the output to
# CSV files.
simulate_data <- function(
  regimes_csv_file,
  n_methylated_reads_csv_file,
  genomic_positions_csv_file,
  n_total_reads_csv_file,
  mu,
  sigma,
  p,
  omega,
  kappa,
  u,
  n_samples, # the number of (biological) samples
  n_cpg_sites, # number of observations/time steps/CpG sites
  lambda, # mean number of reads per CpG site/sample
  randomise_rng_seed, # should the seed for the random number generator be randomised?
  rng_seed # the seed of the random number generator
) {

  aux <- get_known_parameters(
    u = u,
    mu = mu,
    sigma = sigma,
    is_kappa_fixed = TRUE,
    kappa = kappa
  )

  vartheta <- aux$vartheta
  dim_theta <- aux$dim_theta

  theta <- convert_model_parameters_to_theta(
    vartheta,
    p = p,
    omega = omega,
    kappa = kappa
  )

  # The positions of the CpG sites within the genome in our simulated data:
  genomic_positions <- 1:n_cpg_sites

  # Simulates an (n_samples, n_cpg_sites) matrix containing the number of reads
  # for each CpG site and each sample. We simply draw the total number of reads
  # at each site from a Poisson distribution with mean lambda.
  n_total_reads <- matrix(
    rpois(n_samples * n_cpg_sites, lambda = lambda),
    nrow = n_samples,
    ncol = n_cpg_sites
  )

  aux <- simulateDataCpp(
    n_cpg_sites,
    vartheta,
    theta,
    genomic_positions,
    n_total_reads,
    randomise_rng_seed,
    rng_seed
  )

  matrix(
    unlist(aux$latentVariables),
    nrow = 2,
    ncol = n_cpg_sites
  ) %>%
    .[2, ] %>%
    write_regimes_to_csv_file(file = regimes_csv_file)

  matrix(
    unlist(aux$nMethylatedReads),
    nrow = n_samples,
    ncol = n_cpg_sites
  ) %>%
  write_n_methylated_reads_to_csv_file(file = n_methylated_reads_csv_file)

  genomic_positions %>%
    write_genomic_positions_to_csv_file(file = genomic_positions_csv_file)
  n_total_reads %>%
    write_n_total_reads_to_csv_file(file = n_total_reads_csv_file)

  return()
}





