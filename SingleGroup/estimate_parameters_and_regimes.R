################################################################################
#
# R script for estimating the single-group regimes
# as well as (potentially) the single-group model parameters
#
################################################################################

library("argparser")
parser <- arg_parser(description = "")
parser <- add_argument(
  parser,
  "--mu_csv_file",
  help = "CSV file containing parameter vector mu (see documentation)"
)
parser <- add_argument(
  parser,
  "--sigma_csv_file",
  help = "CSV file containing parameter vector sigma (see documentation)"
)
parser <- add_argument(
  parser,
  "--omega_csv_file",
  help = "CSV file containing parameter vector omega (see documentation)"
)
parser <- add_argument(
  parser,
  "--kappa_csv_file",
  help = "CSV file containing parameter vector kappa (see documentation)"
)
parser <- add_argument(
  parser,
  "--u_csv_file",
  help = "CSV file containing parameter vector u (see documentation)"
)
parser <- add_argument(
  parser,
  "--p_csv_file",
  help = "CSV file containing parameter vector p (see documentation)"
)

parser <- add_argument(
  parser,
  "--n_methylated_reads_csv_file",
  help = "path to the n_methylated_reads CSV file"
)
parser <- add_argument(
  parser,
  "--genomic_positions_csv_file",
  help = "path to the genomic_positions CSV file"
)
parser <- add_argument(
  parser,
  "--n_total_reads_csv_file",
  help = "path to the n_total_reads CSV file"
)
parser <- add_argument(
  parser,
  "--regime_probabilities_csv_file",
  help = "path to the output regime_probabilities CSV file"
)

parser <- add_argument(
  parser,
  "--theta_trace_csv_file",
  help = "path to the output theta trace CSV file"
)

parser <- add_argument(
  parser,
  "--is_kappa_fixed",
  type = "boolean",
  help = "should the parameters kappa governing the distance between changepoint be fixed to their initial values, i.e., not estimated from the data (only relevant if --estimate_parameters is set to TRUE)",
  default = TRUE
)
parser <- add_argument(
  parser,
  "--n_particles",
  type = "integer",
  help = "number of particles to use",
  default = 250
)
parser <- add_argument(
  parser,
  "--estimate_regime_probabilities",
  type = "boolean",
  flag = TRUE,
  help = "should the regimes be estimated?",
)
parser <- add_argument(
  parser,
  "--estimate_parameters",
  type = "boolean",
  flag = TRUE,
  help = "should the model parameters be updated throughout the algorithm (otherwise they are fixed at their initial values)?",
)
parser <- add_argument(
  parser,
  "--epsilon",
  type = "double",
  help = "tuning parameter for the online smoothing algorithm (smaller values increase computational cost but reduce the fixed-lag-smoothing approximation error)",
  default = 0.01
)
parser <- add_argument(
  parser,
  "--normalise_gradients",
  type = "boolean",
  help = "should the gradients be normalised to have unit $L_1$-norm?",
  default = FALSE
)
parser <- add_argument(
  parser,
  "--use_adam",
  type = "boolean",
  help = "should the ADAM optimiser be used (otherwise, it is plain stochastic gradient ascent)?",
  default = TRUE
)
parser <- add_argument(
  parser,
  "--n_steps_without_parameter_update",
  type = "integer",
  help = "The number of SMC steps between parameter updates. It seems to be necessary to set this to a value much larger than 1 due to the slow forgetting of the change-point model",
  default = 200
)
parser <- add_argument(
  parser,
  "--learning_rate_exponent",
  type = "double",
  help = "factor governing the learning-rate decay",
  default = 0.1
)
parser <- add_argument(
  parser,
  "--learning_rate_factor",
  type = "double",
  help = "exponent governing the learning-rate decay",
  default = 0.01
)
# ---------------------------------------------------------------------------- #
# That is, write
# a = learning_rate_factor
# b = learning_rate_exponent
# c = n_steps_without_parameter_update
# then the model parameters are only updated at the $t$th SMC step if
# $t \mod c = 0$
# and in this case, the step-size sequence/learing rate used by the
# gradient-ascent/ADAM update is
# $a / i^b$, where $i = t/c$.
# ---------------------------------------------------------------------------- #

parser <- add_argument(
  parser,
  "--root_dir",
  default = './src/r',
  help = "root directory for the src folder containing the R scripts"
)

parser <- add_argument(
  parser,
  "--randomise_rng_seed",
  type = "boolean",
  help = "should the seed for the random number generator be randomised?",
  default = TRUE
)
parser <- add_argument(
  parser,
  "--rng_seed",
  type = "integer",
  help = "seed for the random number generator (only used if randomise_rng_seed is FALSE)",
  default = -73
)

cmd_args <- parse_args(parser)

if (!cmd_args$randomise_rng_seed) {
  set.seed(cmd_args$rng_seed)
}

source(file = file.path(cmd_args$root_dir, "build.R"), chdir = TRUE)
source(file = file.path(cmd_args$root_dir, "input_output_functions.R"), chdir = TRUE)
source(file = file.path(cmd_args$root_dir, "model_functions.R"), chdir = TRUE)


# LOAD MODEL PARAMETERS FROM CSV FILES
# ============================================================================ #

mu    <- read_mu_from_csv_file(cmd_args$mu_csv_file)
sigma <- read_sigma_from_csv_file(cmd_args$sigma_csv_file)
omega <- read_omega_from_csv_file(cmd_args$omega_csv_file)
kappa <- read_kappa_from_csv_file(cmd_args$kappa_csv_file)
u     <- read_u_from_csv_file(cmd_args$u_csv_file)
p     <- read_p_from_csv_file(cmd_args$p_csv_file)


# RUN THE ESTIMATION ALGORITHM
# ============================================================================ #

aux <- get_known_parameters(
  u = u,
  mu = mu,
  sigma = sigma,
  is_kappa_fixed = cmd_args$is_kappa_fixed,
  kappa = kappa
)

vartheta    <- aux$vartheta # some other known parameters
dim_theta   <- aux$dim_theta # number of (transformed) unknown parameters

if (cmd_args$estimate_parameters) { # i.e., if we need to estimate the model parameters from the data
  theta_init <- sampleFromParameterPriorCpp(
    vartheta,
    cmd_args$randomise_rng_seed, # should the seed for the random number generator be randomised?
    cmd_args$rng_seed # the seed of the random number generator (if not randomised)
  )
} else { # i.e., if we already know the model parameters
  theta_init <- convert_model_parameters_to_theta(
    vartheta = vartheta,
    p = p,
    omega = omega,
    kappa = kappa
  )
}

# n_cpg_sites <- read_genomic_positions_from_csv_file(csv_file_genomic_positions) %>% length()
n_methylation_regimes <- length(mu) # number of methylation regimes

csv_file_genomic_positions  <- file.path(cmd_args$genomic_positions_csv_file)
csv_file_n_total_reads      <- file.path(cmd_args$n_total_reads_csv_file)
csv_file_n_methylated_reads <- file.path(cmd_args$n_methylated_reads_csv_file)

# Runs the algorithm. Here, n_samples is the number of biological samples and equal to
# dim(n_total_reads)[1]; n_cpg_sites is the total number of CpG sites and equal
# dim(n_total_reads)[2].
aux <- runOnlineCombinedInferenceCpp(
  vartheta, # hyperparameters and other auxiliary model parameters
  theta_init, # value for the model-parameter vector theta
  read_genomic_positions_from_csv_file(csv_file_genomic_positions), # (1, n_cpg_sites)-matrix of the genomic positions of the CpG sites
  read_n_total_reads_from_csv_file(csv_file_n_total_reads), # (n_samples, n_cpg_sites)-matrix of the total number of reads for each sample at each CpG site
  read_n_methylated_reads_from_csv_file(csv_file_n_methylated_reads), # (n_samples, n_cpg_sites)-matrix of the number of methylated reads for each sample at each CpG site
  cmd_args$n_particles, # number particles used by the SMC algorithm
  1, # the proposal kernel. 0: prior (i.e. bootstrap particle filter). 1: discrete particle filter proposal from Fearnhead (1998), Fearnhead et al. (2003)
  2, # the resampling scheme. 0: multinomial. 1: systematic. 2: optimal finite-state resampling
  cmd_args$estimate_regime_probabilities, # determines if we should estimate the expectations of the regimes under the marginal smoothing distributions
  cmd_args$epsilon, # tuning parameter for the online-smoothing algorithm
  cmd_args$estimate_parameters, # determines if the model parameters should be updated throughout the algorithm (otherwise they are fixed at their initial values)
  cmd_args$normalise_gradients, # should the gradients be normalised according to their $L_1$ norm?
  cmd_args$use_adam, # should the ADAM optimiser be used (otherwise it is plain stochastic gradient ascent)
  cmd_args$n_steps_without_parameter_update, # number of SMC steps without parameter updates
  cmd_args$learning_rate_exponent, # the exponent governing the learning-rate decay
  cmd_args$learning_rate_factor, # the factor governing the learning-rate decay
  cmd_args$randomise_rng_seed, # should the seed for the random number generator be randomised?
  cmd_args$rng_seed # the seed of the random number generator (if not randomised)
)


# store regime estimates
if (cmd_args$estimate_regime_probabilities) {

  aux$regimeProbabilityEstimates %>%
    unlist() %>%
    # matrix(n_methylation_regimes, n_cpg_sites) %>%
    matrix(., n_methylation_regimes, length(.) / n_methylation_regimes) %>% # automatically calculates n_cpg_sites from object size
    t() %>%
    as_tibble(.name_repair = ~ paste0("regime_", 1:n_methylation_regimes)) %>%
    readr::write_csv(
      file = file.path(cmd_args$regime_probabilities_csv_file)
    )
}

# store parameter estimates
if (cmd_args$estimate_parameters) {

  aux$thetaEstimates %>%
    unlist() %>%
    matrix(., dim_theta, length(.) / dim_theta) %>% # automatically calculates n_cpg_sites from object size (used instead of matrix(dim_theta, n_cpg_sites))
    t() %>%
    as_tibble(.name_repair = ~ paste0("theta_", 1:dim_theta)) %>%
    write_csv(file = file.path(cmd_args$theta_trace_csv_file))

  aux$thetaEstimates %>%
    unlist() %>%
    matrix(., dim_theta, length(.) / dim_theta) %>% # automatically calculates n_cpg_sites from object size (used instead of matrix(dim_theta, n_cpg_sites))
    t() %>%
    .[nrow(.), ] %>% # extracts the last row
    convert_theta_to_model_parameters(vartheta, .) -> final_parameter_estimates

  final_parameter_estimates$p %>%
    write_p_to_csv_file(file = cmd_args$p_csv_file)

  final_parameter_estimates$omega %>%
    write_omega_to_csv_file(file = cmd_args$omega_csv_file)

  if (!get_is_kappa_fixed(vartheta)) {
    final_parameter_estimates$kappa %>%
      write_kappa_to_csv_file(file = cmd_args$kappa_csv_file)
  }

}