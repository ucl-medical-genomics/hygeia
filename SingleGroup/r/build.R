library("Rcpp")
library("RcppArmadillo")

# C++ COMPILATION
# ============================================================================ #

do_rebuild <- TRUE # only needs to be set to "TRUE" on first use (can be set to "FALSE" later to avoid re-compilation)

Sys.setenv(
  "PKG_CXXFLAGS" = paste0(
    "-Wall -std=c++11 -I\"", 
    getwd(), 
    "\" -I/usr/include -O3 -ffast-math -fno-finite-math-only -march=native"
  )
)
sourceCpp(
  file.path(getwd(), "cpp", "singleGroup.cpp"), 
  rebuild = do_rebuild, 
  verbose = FALSE
)

# LOAD OTHER AUXILIARY R FUNCTIONS
# ============================================================================ #

source(file = file.path(getwd(), "r", "input_output_functions.R"))
source(file = file.path(getwd(), "r", "model_functions.R"))
