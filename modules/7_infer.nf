#!/usr/bin/env nextflow

process INFER {
    tag "${chrom}-${batch_index}-${inference_seed}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group::v0.1.22'
    publishDir "${params.output_dir}/4_INFER/${chrom}_${batch_index}_${inference_seed}",
        mode: 'copy', pattern: 'nextflow_output/*',
        saveAs: { fn -> fn.replace("nextflow_output", "./") }

    memory { 16.GB + (4.GB * task.attempt) }

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
          path(theta_csv, stageAs: "single_group_estimation/*"),
          val(chrom_alt_from_batch_indexes),
          val(batch_index)
    each inference_seed

    output:
    tuple val(chrom),
          path(regime_probabilities_csv),
          path(theta_trace_csv),
          path(p_csv),
          path(kappa_csv),
          path(omega_csv),
          path(theta_csv),
          path("chrom_${chrom}_${batch_index}_${inference_seed}"), // two_group_results
          val(inference_seed), emit: infer_out
    path "nextflow_output/*"

    script:
    """
    hygeia infer --mu ${params.meteor_mu} --sigma ${params.meteor_sigma} \
        --chrom ${chrom} --single_group_dir ./single_group_estimation \
        --data_dir ./preprocessed_data \
        --results_dir chrom_${chrom}_${batch_index}_${inference_seed} \
        --seed ${inference_seed} --batch ${batch_index} \
        --segment_size ${params.batch_size}

    # prepare the output directory
    cp -r chrom_${chrom}_${batch_index}_${inference_seed} nextflow_output

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p chrom_${chrom}_${batch_index}_${inference_seed}
    touch chrom_${chrom}_${batch_index}_${inference_seed}/dummy_file
    echo "{chrom} {inference_seed} {batch_index}" > chrom_${chrom}_${batch_index}_${inference_seed}/dummy_file
    cp -r chrom_${chrom}_${batch_index}_${inference_seed} nextflow_output

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}
