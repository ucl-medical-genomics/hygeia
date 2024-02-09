#!/usr/bin/env nextflow

params.cpg_file_path = "/home/imoghul/d/hygeia/run_08_02_2024/data/ref/cpg.tsv.gz"
params.sample_sheet = "/home/imoghul/d/hygeia/run_08_02_2024/data/t1d/sample_sheet.csv"
params.output_dir = "results"
params.output_path = '/home/imoghul/d/hygeia/run_08_02_2024/data/t1d/preprocessed2'
params.chrom = 22

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
    path("preprocessed_data/total_cpg_sites_merged_${chrom}.txt"), emit: total_cpg_sites_merged_chr
    path("preprocessed_data/total_cpg_sites_merged.txt"), emit: total_cpg_sites_merged

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

    output:
    path("single_group_estimation/regimes_${chrom}.csv"), emit: regime_probabilities_csv
    path("single_group_estimation/theta_trace_${chrom}.csv"), emit: theta_trace_csv

    script:
    """
    hygeia estimate_parameters_and_regimes \
        --mu 0.95,0.05,0.8,0.2,0.50,0.50 \
        --sigma 0.05,0.05,0.1,0.1,0.1,0.2886751 \
        --u 3 \
        --n_methylated_reads_csv_file ${n_methylated_reads_control_chr} \
        --genomic_positions_csv_file ${positions_chr} \
        --n_total_reads_csv_file ${n_total_reads_control_chr} \
        --regime_probabilities_csv_file single_group_estimation/regimes_${chrom}.csv \
        --theta_trace_csv_file single_group_estimation/theta_trace_${chrom}.csv \
        --estimate_regime_probabilities \
        --estimate_parameters
    """
}

process infer {
    container 'ucl-medical-genomics/hygeia_two_group'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    // preprocessed data output
    path positions_chr
    path n_total_reads_case_chr
    path n_total_reads_control_chr
    path n_methylated_reads_case_chr
    path n_methylated_reads_control_chr
    path cpg_sites_merged_chr
    path total_cpg_sites_merged_chr
    path total_cpg_sites_merged
    // single group output
    path regime_probabilities_csv
    path theta_trace_csv
    val chrom

    output:

    script:
    """
    hygeia infer \
        --mu 0.95,0.05,0.8,0.2,0.50,0.50 \
        --sigma 0.05,0.05,0.1,0.1,0.1,0.2886751 \
        --chrom ${chrom} \
        --single_group_dir ./ \
        --data_dir ./ \
        --results_dir two_group_results
    """

}


workflow {
    preprocess(samples, params.cpg_file_path, params.chrom)
    estimateParametersAndRegimes(
        preprocess.out.n_methylated_reads_control_chr,
        preprocess.out.positions_chr,
        preprocess.out.n_total_reads_control_chr,
        params.chrom
    )
    infer(
        preprocess.out.positions_chr,
        preprocess.out.n_total_reads_case_chr,
        preprocess.out.n_total_reads_control_chr,
        preprocess.out.n_methylated_reads_case_chr,
        preprocess.out.n_methylated_reads_control_chr,
        preprocess.out.cpg_sites_merged_chr,
        preprocess.out.total_cpg_sites_merged_chr,
        preprocess.out.total_cpg_sites_merged,
        estimateParametersAndRegimes.out.regime_probabilities_csv,
        estimateParametersAndRegimes.out.theta_trace_csv,
        params.chrom
    )
}

