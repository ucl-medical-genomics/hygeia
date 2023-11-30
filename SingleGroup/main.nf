process specifyParameters {
  // container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:0.0.1'
  container 'hygiea/single_group'
  if (params.output.specify_param_dir !== null) {
    publishDir "${params.output.specify_param_dir}", mode: 'copy'
  }

  output:
  path(params.mu_csv_filename), emit: mu_csv
  path(params.sigma_csv_filename), emit: sigma_csv
  path(params.omega_csv_filename), emit: omega_csv
  path(params.kappa_csv_filename), emit: kappa_csv
  path(params.u_csv_filename), emit: u_csv
  path(params.p_csv_filename), emit: p_csv

  script:
  """
  hygeia specify_parameters \\
    --mu_csv_file ${params.mu_csv_filename} \\
    --sigma_csv_file ${params.sigma_csv_filename} \\
    --omega_csv_file ${params.omega_csv_filename} \\
    --kappa_csv_file ${params.kappa_csv_filename} \\
    --u_csv_file ${params.u_csv_filename} \\
    --p_csv_file ${params.p_csv_filename} ${params.specify_parameters.extra_args.join(" ")}
  """
}

process simulateData {
  // container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:0.0.1'
  container 'hygiea/single_group'
  if (params.output.simulated_data_dir !== null) {
    publishDir "${params.output.simulated_data_dir}", mode: 'copy'
  }

  input:
  path mu_csv
  path sigma_csv
  path omega_csv
  path kappa_csv
  path u_csv
  path p_csv

  output:
  path(params.regimes_csv), emit: regimes_csv
  path(params.n_methylated_reads_csv), emit: n_methylated_reads_csv
  path(params.genomic_positions_csv), emit: genomic_positions_csv
  path(params.n_total_reads_csv), emit: n_total_reads_csv

  script:
  """
  hygeia simulate_data \\
    --mu_csv_file ${mu_csv} \\
    --sigma_csv_file ${sigma_csv} \\
    --omega_csv_file ${omega_csv} \\
    --kappa_csv_file ${kappa_csv} \\
    --u_csv_file ${u_csv} \\
    --p_csv_file ${p_csv} \\
    --regimes_csv_file ${params.regimes_csv} \\
    --n_methylated_reads_csv_file ${params.n_methylated_reads_csv} \\
    --genomic_positions_csv_file ${params.genomic_positions_csv} \\
    --n_total_reads_csv_file ${params.n_total_reads_csv} \\
    ${params.simulate_data.extra_args.join(" ")}
  """
}

process estimateParametersAndRegimes {
  // container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:0.0.1'
  container 'hygiea/single_group'
  if (params.output.parameters_and_regimes_dir !== null) {
    publishDir "${params.output.parameters_and_regimes_dir}", mode: 'copy'
  }

  input:
  path mu_csv
  path sigma_csv
  path omega_csv
  path kappa_csv
  path u_csv
  path p_csv
  path regimes_csv
  path n_methylated_reads_csv
  path genomic_positions_csv
  path n_total_reads_csv

  output:
  path(params.regime_probabilities_csv), emit: regime_probabilities_csv
  path(params.theta_trace_csv), emit: theta_trace_csv

  script:
  """
  hygeia estimate_parameters_and_regimes \\
    --mu_csv_file ${mu_csv} \\
    --sigma_csv_file ${sigma_csv} \\
    --omega_csv_file ${omega_csv} \\
    --kappa_csv_file ${kappa_csv} \\
    --u_csv_file ${u_csv} \\
    --p_csv_file ${p_csv} \\
    --n_methylated_reads_csv_file ${n_methylated_reads_csv} \\
    --genomic_positions_csv_file ${genomic_positions_csv} \\
    --n_total_reads_csv_file ${n_total_reads_csv} \\
    --regime_probabilities_csv_file ${params.regime_probabilities_csv} \\
    --theta_trace_csv_file ${params.theta_trace_csv} \\
    ${params.parameters_and_regimes.extra_args.join(" ")}
  """
}

workflow {
  specifyParameters()
  simulateData(
    specifyParameters.out.mu_csv,
    specifyParameters.out.sigma_csv,
    specifyParameters.out.omega_csv,
    specifyParameters.out.kappa_csv,
    specifyParameters.out.u_csv,
    specifyParameters.out.p_csv
  )
  estimateParametersAndRegimes(
    specifyParameters.out.mu_csv,
    specifyParameters.out.sigma_csv,
    specifyParameters.out.omega_csv,
    specifyParameters.out.kappa_csv,
    specifyParameters.out.u_csv,
    specifyParameters.out.p_csv,
    simulateData.out.regimes_csv,
    simulateData.out.n_methylated_reads_csv,
    simulateData.out.genomic_positions_csv,
    simulateData.out.n_total_reads_csv
  )
}
