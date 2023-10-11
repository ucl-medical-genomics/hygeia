################################################################################
#
# R script for simulating data (including the total number of reads)
# from the single-group model
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
  "--regimes_csv_file",
  help = "path to the regimes CSV file"
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
  "--number_of_samples",
  type = "integer",
  help = "the number of (biological) samples"
)
parser <- add_argument(
  parser,
  "--number_of_cpg_sites",
  type = "integer",
  help = "number of CpG sites"
)

parser <- add_argument(
  parser,
  "--lambda",
  type = "double",
  help = "expected number of reads per CpG site and sample",
  default = 100
)
parser <- add_argument(
  parser,
  "--root_dir",
  default = './src/r',
  help = "root directory for the src folder containing the R scripts"
)

cmd_args <- parse_args(parser)

# Load auxiliary functions
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


# SIMULATE DATA
# ============================================================================ #

# Simulates the data and writes the output to CSV files.
out <- simulate_data(
  regimes_csv_file = cmd_args$regimes_csv_file,
  n_methylated_reads_csv_file = cmd_args$n_methylated_reads_csv_file,
  genomic_positions_csv_file = cmd_args$genomic_positions_csv_file,
  n_total_reads_csv_file = cmd_args$n_total_reads_csv_file,
  mu = mu,
  sigma = sigma,
  p = p,
  omega = omega,
  kappa = kappa,
  u = u,
  n_samples = cmd_args$number_of_samples, # the number of (biological) samples
  n_cpg_sites = cmd_args$number_of_cpg_sites, # number of observations/time steps/CpG sites
  lambda = cmd_args$lambda # expected number of reads per CpG site and sample
)
