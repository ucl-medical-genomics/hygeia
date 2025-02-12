#!/usr/bin/env nextflow

process PREPROCESS {
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.10'
    publishDir "${params.output_dir}", mode: 'copy'

    memory 16.GB

    input:
    tuple val(case_group),
          val(case_ids),
          path(case_files),
          val(control_group),
          val(control_ids),
          path(control_files)
    path cpg_file_path
    val chrom

    output:
    tuple val(chrom),
          path("preprocessed_data/positions_${chrom}.txt"),
          path("preprocessed_data/n_total_reads_case_${chrom}.txt"),
          path("preprocessed_data/n_total_reads_control_${chrom}.txt"),
          path("preprocessed_data/n_methylated_reads_case_${chrom}.txt"),
          path("preprocessed_data/n_methylated_reads_control_${chrom}.txt"),
          path("preprocessed_data/cpg_sites_merged_${chrom}.txt")

    script:
    def caseIdArgs = case_ids.collect { "--case_id_names '$it'" }.join(" ")
    def caseFileArgs = case_files.collect { "--case_data_path '$it'" }.join(" ")
    def controlIdArgs = control_ids.collect { "--control_id_names '$it'" }.join(" ")
    def controlFileArgs = control_files.collect { "--control_data_path '$it'" }.join(" ")

    """
    hygeia preprocess ${caseIdArgs} ${caseFileArgs} ${controlIdArgs} \
        ${controlFileArgs} --cpg_file_path ${cpg_file_path} \
        --output_path preprocessed_data --chrom ${chrom}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p preprocessed_data
    touch preprocessed_data/positions_${chrom}.txt
    touch preprocessed_data/n_total_reads_case_${chrom}.txt
    touch preprocessed_data/n_total_reads_control_${chrom}.txt
    touch preprocessed_data/n_methylated_reads_case_${chrom}.txt
    touch preprocessed_data/n_methylated_reads_control_${chrom}.txt
    touch preprocessed_data/cpg_sites_merged_${chrom}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}

process ESTIMATE_PARAMETERS_AND_REGIMES {
    container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:v0.1.10'
    publishDir "${params.output_dir}", mode: 'copy', pattern: 'single_group_estimation/*'

    memory 4.GB

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
          path("single_group_estimation/regimes_${chrom}.csv"),
          path("single_group_estimation/theta_trace_${chrom}.csv"),
          path("single_group_estimation/p_${chrom}.csv"),
          path("single_group_estimation/kappa_${chrom}.csv"),
          path("single_group_estimation/omega_${chrom}.csv"),
          path("single_group_estimation/theta_${chrom}.csv")

    script:
    """
    hygeia estimate_parameters_and_regimes \
        --mu ${params.meteor_mu} \
        --sigma ${params.meteor_sigma} \
        --u ${params.min_cpg_sites_between_change_points} \
        --n_methylated_reads_csv_file ${n_methylated_reads_control_chr} \
        --genomic_positions_csv_file ${positions_chr} \
        --n_total_reads_csv_file ${n_total_reads_control_chr} \
        --regime_probabilities_csv_file single_group_estimation/regimes_${chrom}.csv \
        --theta_trace_csv_file single_group_estimation/theta_trace_${chrom}.csv \
        --p_csv_file single_group_estimation/p_${chrom}.csv \
        --kappa_csv_file single_group_estimation/kappa_${chrom}.csv \
        --omega_csv_file single_group_estimation/omega_${chrom}.csv \
        --theta_file single_group_estimation/theta_${chrom}.csv \
        --estimate_regime_probabilities --estimate_parameters

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p single_group_estimation
    touch single_group_estimation/regimes_${chrom}.csv
    touch single_group_estimation/theta_trace_${chrom}.csv
    touch single_group_estimation/p_${chrom}.csv
    touch single_group_estimation/kappa_${chrom}.csv
    touch single_group_estimation/omega_${chrom}.csv
    touch single_group_estimation/theta_${chrom}.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}

process GET_CHROM_SEGMENTS {
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.10'
    publishDir "${params.output_dir}/chrom_segments", mode: 'copy'

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
          path("chrom_segments_${chrom}.csv")

    script:
    """
    hygeia get_chrom_segments --input_file ${positions_chr} --chrom ${chrom} \
        --output_csv chrom_segments_${chrom}.csv --segment_size ${batch_size}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    touch chrom_segments_${chrom}.csv
    """
}

process INFER {
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.10'
    publishDir "${params.output_dir}/two_group_output", mode: 'copy',
        pattern: "infer_out_${chrom}_${inference_seed}/*"

    // memory 16.GB
    memory 250.GB

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
          val(inference_seed)

    script:
    """
    hygeia infer --mu ${params.meteor_mu} --sigma ${params.meteor_sigma} \
        --chrom ${chrom} --single_group_dir ./single_group_estimation \
        --data_dir ./preprocessed_data \
        --results_dir chrom_${chrom}_${batch_index}_${inference_seed} \
        --seed ${inference_seed} --batch ${batch_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p chrom_${chrom}_${batch_number}_${inference_seed}
    touch chrom_${chrom}_${batch_number}_${inference_seed}/dummy_file
    echo "{chrom} {inference_seed} {batch_number}" > chrom_${chrom}_${batch_number}_${inference_seed}/dummy_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}

process AGGREGATE_RESULTS {
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.10'
    publishDir "${params.output_dir}/aggregated", mode: 'copy'

    cpus 4
    memory 24.GB

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

    output:
    tuple val(chrom),
          path("aggregated_out_${chrom}/*")

    script:
    """
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
        --output_dir aggregated_out_${chrom} --num_batches ${params.batches - 1}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p aggregated_out_${chrom}
    touch aggregated_out_${chrom}/dummy_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}

process GET_DMPS {
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.10'
    publishDir "${params.output_dir}/dmps/", mode: 'copy'

    cpus 4
    memory 24.GB

    input:
    tuple val(chrom),
          path(aggregated_out_chr, stageAs: 'aggregated_data/*')

    output:
    path "dmps_${chrom}"

    script:
    """
    hygeia get_dmps --results_dir aggregated_data --output_dir dmps_${chrom} \
        --chrom ${chrom}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dmps_${chrom}
    touch dmps_${chrom}/dummy_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}

workflow {
    ch_chroms = Channel.of(params.chroms.split(','))
    ch_inference_seeds = Channel.of(0..params.num_of_inference_seeds - 1)
    ch_samples = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> tuple(row.group, row.id, row.file) }
        .groupTuple()
        .collect()

    PREPROCESS(ch_samples, params.cpg_file_path, ch_chroms)
    ESTIMATE_PARAMETERS_AND_REGIMES(PREPROCESS.out)

    GET_CHROM_SEGMENTS(ESTIMATE_PARAMETERS_AND_REGIMES.out, params.batch_size)

    // Modify the ch_batch_indexes creation to properly handle file paths
    ch_batch_indexes = GET_CHROM_SEGMENTS.out
        .flatMap { chrom, positions_chr, n_total_reads_case_chr, n_total_reads_control_chr,
                   n_methylated_reads_case_chr, n_methylated_reads_control_chr,
                   cpg_sites_merged_chr, regime_probabilities_csv,theta_trace_csv,
                   p_csv, kappa_csv, omega_csv, theta_csv, segments_file ->
            file(segments_file)
                .splitCsv(header: true)
                .collect { row ->
                    tuple(chrom, positions_chr, n_total_reads_case_chr, n_total_reads_control_chr,
                            n_methylated_reads_case_chr, n_methylated_reads_control_chr,
                            cpg_sites_merged_chr, regime_probabilities_csv,theta_trace_csv,
                            p_csv, kappa_csv, omega_csv, theta_csv, row.chrom, row.segment_index.toInteger())
                }
        }

    INFER(
        ch_batch_indexes,
        ch_inference_seeds
    )

    INFER.out
        .groupTuple()
        .map { r -> tuple(
            r[0],  // chrom
            r[1].first(),  // regime_probabilities_csv
            r[2].first(),  // theta_trace_csv
            r[3].first(),  // p_csv
            r[4].first(),  // kappa_csv
            r[5].first(),  // omega_csv
            r[6].first(),  // theta_csv
            r[7],  // infer_out
            r[8],  // inference_seed
        )}
        .set { merged_infer_outputs }

    AGGREGATE_RESULTS(merged_infer_outputs)
    GET_DMPS(AGGREGATE_RESULTS.out)
}
