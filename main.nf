#!/usr/bin/env nextflow

process PREPROCESS {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.14'
    publishDir "${params.output_dir}/1_PREPROCESS", mode: 'copy',
        pattern: 'nextflow_output/*',
        saveAs: { fn -> fn.replace("nextflow_output", "./") }

    memory { 16.GB + (4.GB * task.attempt) }

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
          path("nextflow_output/positions_${chrom}.txt.gz", arity: '1'),
          path("nextflow_output/n_total_reads_case_${chrom}.txt.gz", arity: '1'),
          path("nextflow_output/n_total_reads_control_${chrom}.txt.gz", arity: '1'),
          path("nextflow_output/n_methylated_reads_case_${chrom}.txt.gz", arity: '1'),
          path("nextflow_output/n_methylated_reads_control_${chrom}.txt.gz", arity: '1'),
          path("nextflow_output/cpg_sites_merged_${chrom}.txt.gz", arity: '1'), emit: preprocessed_data
    path "nextflow_output/versions.yml"

    script:
    def caseIdArgs = case_ids.collect { "--case_id_names '$it'" }.join(" ")
    def caseFileArgs = case_files.collect { "--case_data_path '$it'" }.join(" ")
    def controlIdArgs = control_ids.collect { "--control_id_names '$it'" }.join(" ")
    def controlFileArgs = control_files.collect { "--control_data_path '$it'" }.join(" ")

    """
    hygeia preprocess ${caseIdArgs} ${caseFileArgs} ${controlIdArgs} \
        ${controlFileArgs} --cpg_file_path ${cpg_file_path} \
        --output_path nextflow_output --chrom ${chrom}

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p nextflow_output
    touch nextflow_output/positions_${chrom}.txt.gz
    touch nextflow_output/n_total_reads_case_${chrom}.txt.gz
    touch nextflow_output/n_total_reads_control_${chrom}.txt.gz
    touch nextflow_output/n_methylated_reads_case_${chrom}.txt.gz
    touch nextflow_output/n_methylated_reads_control_${chrom}.txt.gz
    touch nextflow_output/cpg_sites_merged_${chrom}.txt.gz

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}

process ESTIMATE_PARAMETERS {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:v0.1.14'
    publishDir "${params.output_dir}/2_ESTIMATE_PARAMETERS", mode: 'copy',
        pattern: 'nextflow_output/*',
        saveAs: { fn -> fn.replace("nextflow_output", "./") }

    memory { 4.GB * task.attempt }

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
          path("nextflow_output/theta_trace_${chrom}.csv.gz", arity: '1'),
          path("nextflow_output/p_${chrom}.csv.gz", arity: '1'),
          path("nextflow_output/kappa_${chrom}.csv.gz", arity: '1'),
          path("nextflow_output/omega_${chrom}.csv.gz", arity: '1'),
          path("nextflow_output/theta_${chrom}.csv.gz", arity: '1'), emit: single_group_estimation
    path "nextflow_output/versions.yml"

    script:
    """
    mkdir nextflow_output
    hygeia estimate_parameters_and_regimes \
        --mu ${params.meteor_mu} \
        --sigma ${params.meteor_sigma} \
        --u ${params.min_cpg_sites_between_change_points} \
        --n_methylated_reads_csv_file ${n_methylated_reads_control_chr} \
        --genomic_positions_csv_file ${positions_chr} \
        --n_total_reads_csv_file ${n_total_reads_control_chr} \
        --theta_trace_csv_file nextflow_output/theta_trace_${chrom}.csv.gz \
        --p_csv_file nextflow_output/p_${chrom}.csv.gz \
        --kappa_csv_file nextflow_output/kappa_${chrom}.csv.gz \
        --omega_csv_file nextflow_output/omega_${chrom}.csv.gz \
        --theta_file nextflow_output/theta_${chrom}.csv.gz \
        --estimate_parameters

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

process ESTIMATE_REGIMES {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:v0.1.14'
    publishDir "${params.output_dir}/3_ESTIMATE_REGIMES", mode: 'copy',
        pattern: 'nextflow_output/*',
        saveAs: { fn -> fn.replace("nextflow_output", "./") }

    memory { 4.GB * task.attempt }

    input:
    tuple val(chrom),
          path(positions_chr, stageAs: 'preprocessed_data/*'),
          path(n_total_reads_case_chr, stageAs: 'preprocessed_data/*'),
          path(n_total_reads_control_chr, stageAs: 'preprocessed_data/*'),
          path(n_methylated_reads_case_chr, stageAs: 'preprocessed_data/*'),
          path(n_methylated_reads_control_chr, stageAs: 'preprocessed_data/*'),
          path(cpg_sites_merged_chr, stageAs: 'preprocessed_data/*'),
          path(theta_trace_csv, stageAs: "single_group_estimation/*"),
          path(p_csv, stageAs: "single_group_estimation/*"),
          path(kappa_csv, stageAs: "single_group_estimation/*"),
          path(omega_csv, stageAs: "single_group_estimation/*"),
          path(theta_csv, stageAs: "single_group_estimation/*")

    output:
    path "nextflow_output/*"

    script:
    """
    hygeia estimate_parameters_and_regimes \
        --mu ${params.meteor_mu} \
        --sigma ${params.meteor_sigma} \
        --u ${params.min_cpg_sites_between_change_points} \
        --p_input_csv_file ${p_csv} \
        --kappa_input_csv_file ${kappa_csv} \
        --omega_input_csv_file ${omega_csv} \
        --n_methylated_reads_csv_file ${n_methylated_reads_control_chr} \
        --genomic_positions_csv_file ${positions_chr} \
        --n_total_reads_csv_file ${n_total_reads_control_chr} \
        --regime_probabilities_csv_file nextflow_output/regimes_${chrom}.csv.gz \
        --estimate_regime_probabilities

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p nextflow_output
    touch nextflow_output/regimes_${chrom}.csv.gz

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}

process ESTIMATE_PARAMETERS_AND_REGIMES {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:v0.1.14'
    publishDir "${params.output_dir}/2_ESTIMATE_PARAMETERS_AND_REGIMES", mode: 'copy', pattern: 'nextflow_output/*', saveAs: { fn -> fn.replace("nextflow_output", "./") }

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
          path("single_group_estimation/regimes_${chrom}.csv"),
          path("single_group_estimation/theta_trace_${chrom}.csv"),
          path("single_group_estimation/p_${chrom}.csv"),
          path("single_group_estimation/kappa_${chrom}.csv"),
          path("single_group_estimation/omega_${chrom}.csv"),
          path("single_group_estimation/theta_${chrom}.csv"), emit: single_group_estimation
    path "nextflow_output/*", emit: published_output

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
    
    # prepare the output directory
    cp -r single_group_estimation nextflow_output
    gzip nextflow_output/regimes_${chrom}.csv
    gzip nextflow_output/theta_trace_${chrom}.csv

    cat <<-END_VERSIONS > nextflow_output/versions.yml
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
    cp -r single_group_estimation nextflow_output
    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}

process GET_CHROM_SEGMENTS {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.14'
    publishDir "${params.output_dir}/3_GET_CHROM_SEGMENTS", mode: 'copy', pattern: 'nextflow_output/*', saveAs: { fn -> fn.replace("nextflow_output", "./") }

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
          path("chrom_segments_${chrom}.csv"), emit: segments_files
    path "nextflow_output/*", emit: published_output

    script:
    """
    hygeia get_chrom_segments --input_file ${positions_chr} --chrom ${chrom} \
        --output_csv chrom_segments_${chrom}.csv --segment_size ${batch_size}

    # prepare the output directory
    mkdir nextflow_output
    cp chrom_segments_${chrom}.csv nextflow_output

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    touch chrom_segments_${chrom}.csv
    mkdir -p nextflow_output
    cp chrom_segments_${chrom}.csv nextflow_output/
    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}

process INFER {
    tag "${chrom}-${batch_index}-${inference_seed}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.14'
    publishDir "${params.output_dir}/4_INFER/${chrom}_${batch_index}_${inference_seed}", mode: 'copy', pattern: 'nextflow_output/*', saveAs: { fn -> fn.replace("nextflow_output", "./") }

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
    path "nextflow_output/*", emit: published_output

    script:
    """
    hygeia infer --mu ${params.meteor_mu} --sigma ${params.meteor_sigma} \
        --chrom ${chrom} --single_group_dir ./single_group_estimation \
        --data_dir ./preprocessed_data \
        --results_dir chrom_${chrom}_${batch_index}_${inference_seed} \
        --seed ${inference_seed} \
        --batch ${batch_index} --segment_size ${params.batch_size}

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

process AGGREGATE_RESULTS {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.14'
    publishDir "${params.output_dir}/5_AGGREGATE_RESULTS", mode: 'copy', pattern: 'nextflow_output/*', saveAs: { fn -> fn.replace("nextflow_output", "./") }

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

process GET_DMPS {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.14'
    publishDir "${params.output_dir}/6_GET_DMPS/${chrom}", mode: 'copy', 
        pattern: 'nextflow_output/*', saveAs: { fn -> fn.replace("nextflow_output", "./") }

    cpus 4
    memory { 24.GB + (4.GB * task.attempt) }

    input:
    tuple val(chrom),
          path(aggregated_out_chr, stageAs: 'aggregated_data/*')

    output:
    path "nextflow_output/*", emit: published_output

    script:
    """
    hygeia get_dmps --results_dir aggregated_data --output_dir dmps_${chrom} \
        --chrom ${chrom}
    
    cp -r dmps_${chrom} nextflow_output

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dmps_${chrom}
    touch dmps_${chrom}/dummy_file

    cat <<-END_VERSIONS > nextflow_output/versions.yml
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

    if ( !params.two_group ) {
        ESTIMATE_PARAMETERS(PREPROCESS.out.preprocessed_data)
        ESTIMATE_REGIMES(ESTIMATE_PARAMETERS.out.single_group_estimation)
    } else {
        ESTIMATE_PARAMETERS_AND_REGIMES(PREPROCESS.out.preprocessed_data)

        GET_CHROM_SEGMENTS(ESTIMATE_PARAMETERS_AND_REGIMES.out.single_group_estimation, params.batch_size)

        // Modify the ch_batch_indexes creation to properly handle file paths
        ch_batch_indexes = GET_CHROM_SEGMENTS.out.segments_files
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

        // calculate the number of batches
        number_of_batches = GET_CHROM_SEGMENTS.out.segments_files
            .map { _chrom, _positions_chr, _n_total_reads_case_chr, _n_total_reads_control_chr,
                _n_methylated_reads_case_chr, _n_methylated_reads_control_chr,
                _cpg_sites_merged_chr, _regime_probabilities_csv, _theta_trace_csv,
                _p_csv, _kappa_csv, _omega_csv, _theta_csv, segments_file ->
                // Subtract 2 to exclude the header row and 0 indexing
                file(segments_file).readLines().size() - 2  
            }

        INFER(
            ch_batch_indexes,
            ch_inference_seeds
        )

        INFER.out.infer_out
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

        AGGREGATE_RESULTS(merged_infer_outputs, number_of_batches)
        GET_DMPS(AGGREGATE_RESULTS.out.aggregated_out)
    }
}
