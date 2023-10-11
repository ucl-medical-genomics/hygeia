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
parser <- add_argument(
  parser, 
  "outputdir", 
  type = "character", 
  help = "folder in which the simulated data will be stored"
)
parser <- add_argument(
  parser, 
  "n_samples", 
  type = "integer", 
  help = "the number of (biological) samples"
)
parser <- add_argument(
  parser, 
  "n_cpg_sites", 
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

args <- parse_args(parser)
setwd(args$rootdir)


source(file = file.path(args$rootdir, "r", "build.R"))


# LOAD MODEL PARAMETERS FROM CSV FILES
# ============================================================================ #

print(paste0("inputdir", args$inputdir))

mu    <- read_mu_from_csv_file(args$inputdir)
sigma <- read_sigma_from_csv_file(args$inputdir)
omega <- read_omega_from_csv_file(args$inputdir)
kappa <- read_kappa_from_csv_file(args$inputdir)
u     <- read_u_from_csv_file(args$inputdir)
p     <- read_p_from_csv_file(args$inputdir)


# SIMULATE DATA
# ============================================================================ #

# Simulates the data and writes the output to CSV files.
out <- simulate_data(
  outputdir = args$outputdir,
  mu = mu,
  sigma = sigma,
  p = p,
  omega = omega,
  kappa = kappa,
  u = u,
  n_samples = args$n_samples, # the number of (biological) samples
  n_cpg_sites = args$n_cpg_sites, # number of observations/time steps/CpG sites
  lambda = args$lambda # expected number of reads per CpG site and sample
)