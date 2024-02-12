#!/usr/bin/env nextflow

params.cpg_file_path = "/home/imoghul/d/hygeia/data/ref/cpg.tsv.gz"
params.sample_sheet = "/home/imoghul/d/hygeia/data/t1d/sample_sheet.csv"
params.output_dir = "results"
params.meteor_mu = "0.95,0.05,0.8,0.2,0.50,0.50"
params.meteor_sigma = "0.05,0.05,0.1,0.1,0.1,0.2886751"
params.min_cpg_sites_between_change_points = 3

Channel
    .of(1..22)
    .set { chroms }

Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true, sep: ',', strip: true)
    .map { row -> tuple(row.group, row.id, row.file) }
    .groupTuple()
    .collect()
    .set { samples }

process preprocess {
    container 'ucl-medical-genomics/hygeia_two_group'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(case_group), val(case_ids), path(case_files), val(control_group), val(control_ids), path(control_files)
    path cpg_file_path
    val chrom

    output:
    path("preprocessed_data/positions_${chrom}.txt"), emit: positions_chr
    path("preprocessed_data/n_total_reads_case_${chrom}.txt"), emit: n_total_reads_case_chr
    path("preprocessed_data/n_total_reads_control_${chrom}.txt"), emit: n_total_reads_control_chr
    path("preprocessed_data/n_methylated_reads_case_${chrom}.txt"), emit: n_methylated_reads_case_chr
    path("preprocessed_data/n_methylated_reads_control_${chrom}.txt"), emit: n_methylated_reads_control_chr
    path("preprocessed_data/cpg_sites_merged_${chrom}.txt"), emit: cpg_sites_merged_chr
    val chrom, emit: chrom

    script:
    def caseIdArgs = case_ids.collect { "--case_id_names '$it'" }.join(" ")
    def caseFileArgs = case_files.collect { "--case_data_path '$it'" }.join(" ")
    def controlIdArgs = control_ids.collect { "--control_id_names '$it'" }.join(" ")
    def controlFileArgs = control_files.collect { "--control_data_path '$it'" }.join(" ")

    """
    hygeia preprocess \
        ${caseIdArgs} \
        ${caseFileArgs} \
        ${controlIdArgs} \
        ${controlFileArgs} \
        --cpg_file_path ${cpg_file_path} \
        --output_path preprocessed_data \
        --chrom ${chrom}
    """
}

process estimateParametersAndRegimes {
    // container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:0.0.1'
    container 'ucl-medical-genomics/hygeia_single_group'
    publishDir "${params.output_dir}/estimatedParamatersAndRegimes", mode: 'copy'

    input:
    path n_methylated_reads_control_chr
    path positions_chr
    path n_total_reads_control_chr
    val chrom
    // just passing this forward...
    path n_methylated_reads_case_chr
    path n_total_reads_case_chr
    path cpg_sites_merged_chr


    output:
    path n_methylated_reads_control_chr, emit: n_methylated_reads_control_chr
    path positions_chr, emit: positions_chr
    path n_total_reads_control_chr, emit: n_total_reads_control_chr
    path n_methylated_reads_case_chr, emit: n_methylated_reads_case_chr
    path n_total_reads_case_chr, emit: n_total_reads_case_chr
    path cpg_sites_merged_chr, emit: cpg_sites_merged_chr
    path("single_group_estimation/regimes_${chrom}.csv"), emit: regime_probabilities_csv
    path("single_group_estimation/theta_trace_${chrom}.csv"), emit: theta_trace_csv
    path("single_group_estimation/p_${chrom}.csv"), emit: p_csv
    path("single_group_estimation/kappa_${chrom}.csv"), emit: kappa_csv
    path("single_group_estimation/omega_${chrom}.csv"), emit: omega_csv
    path("single_group_estimation/theta_${chrom}.csv"), emit: theta_csv
    val chrom, emit: chrom

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
        --estimate_regime_probabilities \
        --estimate_parameters
    """
}

process infer {
    container 'ucl-medical-genomics/hygeia_two_group'
    publishDir "${params.output_dir}/two_group", mode: 'copy'

    input:
    // preprocessed data output
    path positions_chr
    path n_total_reads_case_chr
    path n_total_reads_control_chr
    path n_methylated_reads_case_chr
    path n_methylated_reads_control_chr
    path cpg_sites_merged_chr
    // single group output
    path regime_probabilities_csv
    path theta_trace_csv
    path p_csv
    path kappa_csv
    path omega_csv
    path theta_csv
    val chrom

    output:
    path("two_group_results_${chrom}/*")
    val chrom, emit: chrom

    script:
    """
    hygeia infer \
        --mu ${params.meteor_mu} \
        --sigma ${params.meteor_sigma} \
        --chrom ${chrom} \
        --single_group_dir ./ \
        --data_dir ./ \
        --results_dir two_group_results_${chrom}
    """
}

workflow {
    preprocess(samples, params.cpg_file_path, chroms)
    estimateParametersAndRegimes(
        preprocess.out.n_methylated_reads_control_chr,
        preprocess.out.positions_chr,
        preprocess.out.n_total_reads_control_chr,
        preprocess.out.chrom,
        preprocess.out.n_methylated_reads_case_chr,
        preprocess.out.n_total_reads_case_chr
    )
    infer(
        estimateParametersAndRegimes.out.positions_chr,
        estimateParametersAndRegimes.out.n_total_reads_case_chr,
        estimateParametersAndRegimes.out.n_total_reads_control_chr,
        estimateParametersAndRegimes.out.n_methylated_reads_case_chr,
        estimateParametersAndRegimes.out.n_methylated_reads_control_chr,
        estimateParametersAndRegimes.out.cpg_sites_merged_chr,
        estimateParametersAndRegimes.out.regime_probabilities_csv,
        estimateParametersAndRegimes.out.theta_trace_csv,
        estimateParametersAndRegimes.out.p_csv,
        estimateParametersAndRegimes.out.kappa_csv,
        estimateParametersAndRegimes.out.omega_csv,
        estimateParametersAndRegimes.out.theta_csv,
        estimateParametersAndRegimes.out.chrom
    )
}
