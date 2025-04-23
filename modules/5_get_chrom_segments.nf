#!/usr/bin/env nextflow

process GET_CHROM_SEGMENTS {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.19'
    publishDir "${params.output_dir}/3_GET_CHROM_SEGMENTS", mode: 'copy',
        pattern: 'nextflow_output/*',
        saveAs: { fn -> fn.replace("nextflow_output", "./") }

    input:
    tuple val(chrom),
          path(positions_chr, stageAs: 'preprocessed_data/*'),
          path(n_total_reads_case_chr, stageAs: 'preprocessed_data/*'),
          path(n_total_reads_control_chr, stageAs: 'preprocessed_data/*'),
          path(n_methylated_reads_case_chr, stageAs: 'preprocessed_data/*'),
          path(n_methylated_reads_control_chr, stageAs: 'preprocessed_data/*'),
          path(cpg_sites_merged_chr, stageAs: 'preprocessed_data/*'),
          path(regime_probabilities_csv, stageAs: "single_group_estimation/*"),
          path(theta_trace_csv, stageAs: "single_group_estimation/*"),
          path(p_csv, stageAs: "single_group_estimation/*"),
          path(kappa_csv, stageAs: "single_group_estimation/*"),
          path(omega_csv, stageAs: "single_group_estimation/*"),
          path(theta_csv, stageAs: "single_group_estimation/*")
    val(batch_size)

    output:
    tuple val(chrom),
          path(positions_chr),
          path(n_total_reads_case_chr),
          path(n_total_reads_control_chr),
          path(n_methylated_reads_case_chr),
          path(n_methylated_reads_control_chr),
          path(cpg_sites_merged_chr),
          path(regime_probabilities_csv),
          path(theta_trace_csv),
          path(p_csv),
          path(kappa_csv),
          path(omega_csv),
          path(theta_csv),
          path("nextflow_output/chrom_segments_${chrom}.csv"), emit: segments_files
    path "nextflow_output/versions.yml"

    script:
    """
    mkdir nextflow_output
    hygeia get_chrom_segments --input_file ${positions_chr} --chrom ${chrom} \
        --output_csv nextflow_output/chrom_segments_${chrom}.csv \
        --segment_size ${batch_size}

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p nextflow_output
    touch nextflow_output/chrom_segments_${chrom}.csv
    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}