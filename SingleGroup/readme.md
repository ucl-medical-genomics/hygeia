# Single Group Analysis

This is code for when you are comparing....

**TODO**
> read this paper for more info


# Installation.

> Ismail will do this


# Run With examplar data

Prerequisites: R and the R packages "tidyverse", "argparser", "Rcpp" and 
"RcppArmadillo" need to be installed. That is, after installing R, simply 
copy-and-paste the following code into the R terminal:
```
install.packages(c("tidyverse", "argparser", "Rcpp", "RcppArmadillo"))
```
For the latter two packages, the system also needs to have, e.g., the GCC compiler. 


Hereafter, we will call the directory containing this readme file the "root" directory, i.e., the root directory should contain the files 

- specify_parameters.R,
- simulate_data.R,
- estimate_parameters_and_regimes.R,
- readme.md,
- example.sh,

and also the folders

- cpp,
- r.

The file example.sh demonstrates how to use this programme. Simply edit the 
following path at the top of example.sh so that it contains the path to the root
directory.

```
rootdir="</path to root directory goes here>"
```
For example, if your root directory is "/home/john_smith/single-group_model" then change this line to
```
rootdir="/home/john_smith/single-group_model"
```
You may then run the script example.sh (e.g., by simply copying-and-pasting all of its contents into the command line).

The script will run the single-group algorithm to estimate regime probabilities
and potentially the unknown model parameters presuming that you have a data set
stored in three CSV files as follows (note that the names also need to match):

1. genomic_positions.csv: holds the genomic positions of the CpG sites (in ascending order) in the first column, e.g., cell (2, 1) holds the genomic position of the second CpG site.

2. n_total_reads.csv: holds the total number of reads for each CpG site (rows) and each sample (columns), e.g. the number in cell (2, 3) is the total number of reads at CpG Site 2 in Sample 3.

3. n_methylated_reads.csv: holds the number of methylated reads for each CpG site (rows) and each sample (columns), e.g. the number in cell (2, 3) is the number of methylated reads at CpG Site 2 in Sample 3

The script also provides the option of simulating such data from the model
in case you want to test the algorithm without providing your own data set.

See example.sh for further explanation.


