process simulateData {
  // container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:0.0.1'
  container 'ucl-medical-genomics/hygeia_single_group'
  publishDir "${params.output_dir}/simulatedData", mode: 'copy'

  output:
  path(params.regimes_csv), emit: regimes_csv
  path(params.n_methylated_reads_csv), emit: n_methylated_reads_csv
  path(params.genomic_positions_csv), emit: genomic_positions_csv
  path(params.n_total_reads_csv), emit: n_total_reads_csv

  script:
  """
  hygeia simulate_data \\
    --mu 0.99,0.01,0.80,0.20,0.50,0.50 \\
    --sigma 0.05,0.05,0.20,0.20,0.20,0.2886751 \\
    --u 3 \\
    --kappa 2,2,2,2,2,2 \\
    --omega 0.995,0.975,0.950,0.925,0.900,0.900 \\
    --n_methylated_reads_csv_file ${params.n_methylated_reads_csv} \\
    --genomic_positions_csv_file ${params.genomic_positions_csv} \\
    --n_total_reads_csv_file ${params.n_total_reads_csv} \\
    --regimes_csv_file ${params.regimes_csv} \\
    ${params.simulate_data_extra_args.join(" ")}
  """
}

process estimateParametersAndRegimes {
  // container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:0.0.1'
  container 'ucl-medical-genomics/hygeia_single_group'
  publishDir "${params.output_dir}/estimatedParamatersAndRegimes", mode: 'copy'

  input:
  path n_methylated_reads_csv
  path genomic_positions_csv
  path n_total_reads_csv

  output:
  path(params.regime_probabilities_csv), emit: regime_probabilities_csv
  path(params.theta_trace_csv), emit: theta_trace_csv

  script:
  """
  hygeia estimate_parameters_and_regimes \\
    --mu 0.99,0.01,0.80,0.20,0.50,0.50 \\
    --sigma 0.05,0.05,0.20,0.20,0.20,0.2886751 \\
    --u 3 \\
    --n_methylated_reads_csv_file ${n_methylated_reads_csv} \\
    --genomic_positions_csv_file ${genomic_positions_csv} \\
    --n_total_reads_csv_file ${n_total_reads_csv} \\
    --regime_probabilities_csv_file ${params.regime_probabilities_csv} \\
    --theta_trace_csv_file ${params.theta_trace_csv} \\
    ${params.parameters_and_regimes.extra_args.join(" ")}
  """
}

workflow {
  simulateData()
  estimateParametersAndRegimes(
    simulateData.out.n_methylated_reads_csv,
    simulateData.out.genomic_positions_csv,
    simulateData.out.n_total_reads_csv
  )
}
