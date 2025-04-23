#!/usr/bin/env nextflow

process ESTIMATE_PARAMETERS_AND_REGIMES {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:v0.1.18'
    publishDir "${params.output_dir}/2_ESTIMATE_PARAMETERS_AND_REGIMES",
        mode: 'copy', pattern: 'nextflow_output/*',
        saveAs: { fn -> fn.replace("nextflow_output", "./") }

    memory { 8.GB + (4.GB * task.attempt) }

    input:
    tuple val(chrom),
          path(positions_chr, stageAs: "preprocessed_data/*"),
          path(n_total_reads_case_chr, stageAs: "preprocessed_data/*"),
          path(n_total_reads_control_chr, stageAs: "preprocessed_data/*"),
          path(n_methylated_reads_case_chr, stageAs: "preprocessed_data/*"),
          path(n_methylated_reads_control_chr, stageAs: "preprocessed_data/*"),
          path(cpg_sites_merged_chr, stageAs: "preprocessed_data/*")

    output:
    tuple val(chrom),
          path(positions_chr),
          path(n_total_reads_case_chr),
          path(n_total_reads_control_chr),
          path(n_methylated_reads_case_chr),
          path(n_methylated_reads_control_chr),
          path(cpg_sites_merged_chr),
          path("nextflow_output/regimes_${chrom}.csv.gz"),
          path("nextflow_output/theta_trace_${chrom}.csv.gz"),
          path("nextflow_output/p_${chrom}.csv.gz"),
          path("nextflow_output/kappa_${chrom}.csv.gz"),
          path("nextflow_output/omega_${chrom}.csv.gz"),
          path("nextflow_output/theta_${chrom}.csv.gz"), emit: single_group_estimation
    path "nextflow_output/versions.yml"

    script:
    """
    hygeia estimate_parameters_and_regimes \
        --mu ${params.meteor_mu} \
        --sigma ${params.meteor_sigma} \
        --u ${params.min_cpg_sites_between_change_points} \
        --n_methylated_reads_csv_file ${n_methylated_reads_control_chr} \
        --genomic_positions_csv_file ${positions_chr} \
        --n_total_reads_csv_file ${n_total_reads_control_chr} \
        --regime_probabilities_csv_file nextflow_output/regimes_${chrom}.csv.gz \
        --theta_trace_csv_file nextflow_output/theta_trace_${chrom}.csv.gz \
        --p_csv_file nextflow_output/p_${chrom}.csv.gz \
        --kappa_csv_file nextflow_output/kappa_${chrom}.csv.gz \
        --omega_csv_file nextflow_output/omega_${chrom}.csv.gz \
        --theta_file nextflow_output/theta_${chrom}.csv.gz \
        --estimate_regime_probabilities --estimate_parameters
    
    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p nextflow_output
    touch nextflow_output/regimes_${chrom}.csv.gz
    touch nextflow_output/theta_trace_${chrom}.csv.gz
    touch nextflow_output/p_${chrom}.csv.gz
    touch nextflow_output/kappa_${chrom}.csv.gz
    touch nextflow_output/omega_${chrom}.csv.gz
    touch nextflow_output/theta_${chrom}.csv.gz
    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}
