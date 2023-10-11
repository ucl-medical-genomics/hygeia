mkdir -p ./example2/params;
Rscript specify_parameters.R \
  --mu_csv_file ./example2/params/mu.csv \
  --sigma_csv_file ./example2/params/sigma.csv \
  --omega_csv_file ./example2/params/omega.csv \
  --kappa_csv_file ./example2/params/kappa.csv \
  --u_csv_file ./example2/params/u.csv \
  --p_csv_file ./example2/params/p.csv

Rscript specify_parameters.R $rootdir $inputdir

mkdir ./example2/simulated_data;
Rscript simulate_data.R \
  --mu_csv_file ./example2/params/mu.csv \
  --sigma_csv_file ./example2/params/sigma.csv \
  --omega_csv_file ./example2/params/omega.csv \
  --kappa_csv_file ./example2/params/kappa.csv \
  --u_csv_file ./example2/params/u.csv \
  --p_csv_file ./example2/params/p.csv \
  --regimes_csv_file ./example2/simulated_data/regimes.csv \
  --n_methylated_reads_csv_file ./example2/simulated_data/n_methylated_reads.csv \
  --genomic_positions_csv_file ./example2/simulated_data/genomic_positions.csv \
  --n_total_reads_csv_file ./example2/simulated_data/n_total_reads.csv \
  --number_of_samples 2 \
  --number_of_cpg_sites 250 \
  --randomise_rng_seed true \
  --rng_seed 73

mkdir ./example2/final_output
Rscript estimate_parameters_and_regimes.R \
  --mu_csv_file ./example2/params/mu.csv \
  --sigma_csv_file ./example2/params/sigma.csv \
  --omega_csv_file ./example2/params/omega.csv \
  --kappa_csv_file ./example2/params/kappa.csv \
  --u_csv_file ./example2/params/u.csv \
  --p_csv_file ./example2/params/p.csv \
  --n_methylated_reads_csv_file ./example2/simulated_data/n_methylated_reads.csv \
  --genomic_positions_csv_file ./example2/simulated_data/genomic_positions.csv \
  --n_total_reads_csv_file ./example2/simulated_data/n_total_reads.csv \
  --regime_probabilities_csv_file ./example2/final_output/regimes.csv \
  --theta_trace_csv_file ./example2/final_output/theta_trace.csv \
  --estimate_regime_probabilities \
  --estimate_parameters \
  --randomise_rng_seed true \
  --rng_seed 73
