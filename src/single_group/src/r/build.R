library("Rcpp")
library("RcppArmadillo")

if (exists("cmd_args")) {
  # If called from source function with chdir = TRUE
  # Note, that the src_root is relative to current file.
  cmd_args$skip_build <- TRUE
  cmd_args$src_root_dir <- '../'
} else {
  library("argparser")
  parser <- arg_parser(description = "")
  parser <- add_argument(
    parser,
    "--skip_build",
    flag = TRUE,
    help = "only needs to be built on first use (can be set to avoid re-compilation)"
  )

  parser <- add_argument(
    parser,
    "--src_root_dir",
    default = './src',
    help = "root directory for the src folder containing the R scripts"
  )
  cmd_args <- parse_args(parser)
}

# C++ COMPILATION
# ============================================================================ #

Sys.setenv(
  "PKG_CXXFLAGS" = paste0(
    "-Wall -std=c++11 -I\"",
    cmd_args$src_root_dir,
    "\" -I/usr/include -O3 -ffast-math -fno-finite-math-only -march=native"
  )
)
sourceCpp(
  file.path(cmd_args$src_root_dir, "cpp", "singleGroup.cpp"),
  rebuild = !cmd_args$skip_build,
  verbose = FALSE
)
