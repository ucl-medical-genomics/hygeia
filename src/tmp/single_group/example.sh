#!/bin/bash

##########################################################
#
# Runs the single-group algorithm on a data set simulated
# from the model or on a user-provided (e.g., real) data
# set
#
##########################################################


##########################################################
# SETUP
##########################################################

# Before you proceed, perform the following steps:

# 1. Edit this path so that it is the directory in which the present file is located:

rootdir="${PWD}/src"

# 2. If desired, change the following paths which give the directories
# from which the algorithm reads its input (such as the data) and to which
# the algorithm writes its output. You can just leave these lines unaltered.
# If the folders below do not exist, they will be automatically created later.

inputdir="$rootdir/../example_data/input"
outputdir="$rootdir/../example_data/output"

# 3. By default, this script will simulate data from the model. If you want to
# use your own data:
#
# a. Comment out all the lines in the "---" delineated box below as these
# simulate a suitable data set from the model.
#
# b. Ensure that your own data are located in inputdir and are stored in the
# following three CSV files (the names also need to match exactly):
# (1) genomic_positions.csv: holds the genomic positions of the CpG sites (in ascending order) in the first column, e.g., cell (2, 1) holds the genomic position of the second CpG site
# (2) n_total_reads.csv: holds the total number of reads for each CpG site (rows) and each sample (columns), e.g. the number in cell (2, 3) is the total number of reads at CpG Site 2 in Sample 3
# (3) n_methlyated_reads.csv: holds the number of methylated reads for each CpG site (rows) and each sample (columns), e.g. the number in cell (2, 3) is the number of methylated reads at CpG Site 2 in Sample 3

# 4. Set the following argument, which determines whether some model parameters
# (the regime-transition matrix P and the sojourn-time parameters omega) will
# be estimated from the data.

estimate_parameters=true # should the model parameters be estimated? Otherwise, all parameters will fixed to those specified by the user (default: true)

# 5. Edit the file specify_parameters.R which specifies the model parameters
# used by the algorithm. In particular, adjust the vectors mu and sigma to your
# needs as these determine the regime configurations. If you want a different
# number of regimes, n, then all the vectors that are currently of length 6,
# must be of length n (and p must be an (n, n) row-stochastic matrix with
# zeros on its diagonal).
#
# NOTE: If you have set "estimate_parameters=true" above, then the
# algorithm will not make use of the parameters omega and p specified in
# specify_parameters.R as these parameters will then instead be estimated from
# the data. If you have set "estimate_parameters=false" above, then
# all parameter values will be taken from specify_parameters.R.

# 6. Specify the following optional arguments (there are other optional
# arguments but only those listed here will be relevant to most users):

n_particles=300 # the number of particles; larger values increase computational cost but improve accuracy (default: 250)
estimate_regime_probabilities=true # should the regime probabilities be estimated (default: true)
randomise_rng_seed=false # shoud the random-number-generator seed be randomised?
rng_seed=73 # the random-number-generator seed (default: 73); only used if randomise_rng_seed is false


##########################################################
# RUN ESTIMATION ALGORITHM ON REAL OR SIMULATED DATA
##########################################################

cd $rootdir

# Make sure input and output directories exist:
mkdir -p $inputdir
mkdir -p $outputdir

# Specify the model parameters:
Rscript specify_parameters.R $rootdir $inputdir

# ---------------------------------------------------------------------------- #
# Comment out everything between the dashed lines if you want to provide your
# own data set.
#
n_samples=2 # number of (biological) samples
n_cpg_sites=250 # number of CpG sites
lambda=100 # optional parameter governing the expected number of reads per CpG site and sample (default: 100)

# Specify the model parameters for simulating data:
Rscript specify_parameters.R $rootdir $inputdir

# Simulate data using the parameters specified by the previous script:
Rscript simulate_data.R $rootdir $inputdir $outputdir $n_samples $n_cpg_sites --lambda $lambda --randomise_rng_seed $randomise_rng_seed --rng_seed $rng_seed

cd $outputdir # folder containing the simulated data
mv genomic_positions.csv n_total_reads.csv n_methylated_reads.csv $inputdir # move simulated data to inputdir
cd $rootdir # back to root directory
# ---------------------------------------------------------------------------- #

# Runs the single-group estimation algorithm:
Rscript estimate_parameters_and_regimes.R $rootdir $inputdir $outputdir --n_particles $n_particles --estimate_regime_probabilities $estimate_regime_probabilities --estimate_parameters $estimate_parameters --randomise_rng_seed $randomise_rng_seed --rng_seed $rng_seed

