#!/usr/bin/env nextflow

include { PREPROCESS } from './modules/1_preprocess.nf'
include { ESTIMATE_PARAMETERS } from './modules/2_estimate_parameters.nf'
include { ESTIMATE_REGIMES } from './modules/3_estimate_regimes.nf'
include { GENERATE_SINGLE_GROUP_BED_FILES } from './modules/4_generate_single_group_bed_files.nf'

include { ESTIMATE_PARAMETERS_AND_REGIMES } from './modules/5_estimate_parameters_and_regimes.nf'
include { GET_CHROM_SEGMENTS } from './modules/6_get_chrom_segments.nf'
include { INFER } from './modules/7_infer.nf'
include { AGGREGATE_RESULTS } from './modules/8_aggregate_results.nf'
include { GET_DMPS } from './modules/9_get_dmps.nf'

workflow {
    ch_inference_seeds = Channel.of(0..params.num_of_inference_seeds - 1)
    ch_samples = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row -> tuple(row.chrom, row.sample_name, row.sample_file) }

    PREPROCESS(ch_samples)

    if ( !params.two_group ) {
        ESTIMATE_PARAMETERS(PREPROCESS.out.preprocessed_data)
        ESTIMATE_REGIMES(ESTIMATE_PARAMETERS.out.single_group_estimation)
        GENERATE_SINGLE_GROUP_BED_FILES(ESTIMATE_REGIMES.out.single_group_regimes)
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
