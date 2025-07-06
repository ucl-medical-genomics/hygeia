#!/usr/bin/env nextflow

process AGGREGATE_RESULTS {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.23'
    publishDir "${params.output_dir}/5_AGGREGATE_RESULTS", mode: 'copy',
        pattern: 'nextflow_output/*',
        saveAs: { fn -> fn.replace("nextflow_output", "./") }

    cpus 4
    memory { 24.GB + (4.GB * task.attempt) }

    input:
    tuple val(chrom),
          path(regime_probabilities_csv, stageAs: "single_group_estimation/*"),
          path(theta_trace_csv, stageAs: "single_group_estimation/*"),
          path(p_csv, stageAs: "single_group_estimation/*"),
          path(kappa_csv, stageAs: "single_group_estimation/*"),
          path(omega_csv, stageAs: "single_group_estimation/*"),
          path(theta_csv, stageAs: "single_group_estimation/*"),
          path("two_group_results_${chrom}/*"), // two_group_results
          val(inference_seed)
    val(number_of_batches)

    output:
    tuple val(chrom),
          path("aggregated_out_${chrom}/*"), emit: aggregated_out
    path 'nextflow_output/*', emit: published_output

    script:
    """
    # structure is:
    # two_group_results_{chrom}/
    #      chrom_{chrom}_{batch_number}_{inference_seed}/
    #          chrom_{chrom}_{batch_number}/
    #              files...
    for file in two_group_results_${chrom}/chrom_${chrom}_*/*/*; do
        if [ -f "\$file" ]; then
            BASE_NAME=\$(basename "\$file")
            DEST_DIR=\$(basename "\$(dirname "\${file}")")
            DEST_FILE=merged_out_${chrom}/"\${DEST_DIR}/\${BASE_NAME}"

            mkdir -p merged_out_${chrom}/"\${DEST_DIR}"
            # Only create link if it doesn't already exist
            if [ ! -e "\${DEST_FILE}" ]; then
                ln -s "../../\$file" "\${DEST_FILE}"
            fi
        fi
    done

    hygeia aggregate --results_dir merged_out_${chrom} --chrom ${chrom} \
        --seeds ${params.num_of_inference_seeds} --num_particles 2400 \
        --output_dir aggregated_out_${chrom} --num_batches ${number_of_batches}

    # prepare the output directory
    cp -r aggregated_out_${chrom} nextflow_output

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p aggregated_out_${chrom}
    touch aggregated_out_${chrom}/dummy_file
    cp -r aggregated_out_${chrom} nextflow_output
    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}
