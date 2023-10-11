################################################################################
#
# R script for specifying some of the vector- and matrix-valued 
# model parameters in R which are then automatically written  
# into the right format of to CSV files for use within the data-simulation or
# inference algorithms. This avoids having to pass vector- or matrix valued
# argument via the command line.
#
################################################################################

library("argparser")
parser <- arg_parser(description = "---")

parser <- add_argument(
  parser, 
  "rootdir", 
  type = "character", 
  help = "root directory, i.e. the directory which holds this script file"
)
parser <- add_argument(
  parser, 
  "inputdir", 
  type = "character", 
  help = "folder from which some of the known/fixed parameter values will be read"
)
args <- parse_args(parser)
setwd(args$rootdir)
source(file = file.path(getwd(), "r", "input_output_functions.R"))

mu    <- c(0.99, 0.01, 0.80, 0.20, 0.50, 0.50)

sigma <- c(0.05, 0.05, 0.20, 0.20, 0.20, 1 / sqrt(12))

u <- 2

kappa  <- rep(2, times = 6) # suitable default

# transition matrix for the regimes:
p <- matrix(
  c(rep(c(0, rep(1 / 5, times = length(mu))), times = length(mu) - 1), 0),
  nrow = length(mu),
  ncol = length(mu)
) 

# Success-probability parameters of the regime-specific negative-binomial 
# distributions governing the function h() which determines the 
# change-point probabilities:
omega <- c(0.995, 0.975, 0.950, 0.925, 0.900, 0.900) 

mu %>% write_mu_to_csv_file(dir = args$inputdir)
sigma %>% write_sigma_to_csv_file(dir = args$inputdir)
kappa %>% write_kappa_to_csv_file(dir = args$inputdir)
omega %>% write_omega_to_csv_file(dir = args$inputdir)
u %>% write_u_to_csv_file(dir = args$inputdir)
p %>% write_p_to_csv_file(dir = args$inputdir)
