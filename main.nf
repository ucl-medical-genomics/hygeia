#!/usr/bin/env nextflow

include { SINGLE_GRP_PREPROCESS } from './modules/single_group/1_preprocess.nf'
include { ESTIMATE_PARAMETERS } from './modules/single_group/2_estimate_parameters.nf'
include { ESTIMATE_REGIMES } from './modules/single_group/3_estimate_regimes.nf'
include { GENERATE_SINGLE_GROUP_BED_FILES } from './modules/single_group/4_generate_single_group_bed_files.nf'

// TWO group
include { TWO_GRP_PREPROCESS } from './modules/two_group/1_preprocess.nf'
include { ESTIMATE_PARAMETERS_AND_REGIMES } from './modules/two_group/2_estimate_parameters_and_regimes.nf'
include { GET_CHROM_SEGMENTS } from './modules/two_group/3_get_chrom_segments.nf'
include { INFER } from './modules/two_group/4_infer.nf'
include { AGGREGATE_RESULTS } from './modules/two_group/5_aggregate_results.nf'
include { GET_DMPS } from './modules/two_group/6_get_dmps.nf'

workflow {
    ch_chroms = Channel.of(params.chroms.split(','))
    ch_inference_seeds = Channel.of(0..params.num_of_inference_seeds - 1)

    if ( !params.two_group ) {
        ch_samples_chroms = Channel
            .fromPath(params.sample_sheet)
            .splitCsv(header: true, sep: ',', strip: true)
            .map { row -> tuple(row.id, row.file) }
            .combine(ch_chroms) // Combine samples with chromosomes to run separately for each chrom + chr combination

        SINGLE_GRP_PREPROCESS(ch_samples_chroms, params.cpg_file_path)
        ESTIMATE_PARAMETERS(SINGLE_GRP_PREPROCESS.out.preprocessed_data)
        ESTIMATE_REGIMES(ESTIMATE_PARAMETERS.out.single_group_estimation)
        GENERATE_SINGLE_GROUP_BED_FILES(ESTIMATE_REGIMES.out.single_group_regimes)
    } else {
        // TWO GROUP ANALYSIS
        ch_samples = Channel
            .fromPath(params.sample_sheet)
            .splitCsv(header: true, sep: ',', strip: true)
            .map { row -> tuple(row.group, row.id, row.file) }
            .groupTuple()
            .collect()

        TWO_GRP_PREPROCESS(ch_samples, params.cpg_file_path, ch_chroms)

        ESTIMATE_PARAMETERS_AND_REGIMES(TWO_GRP_PREPROCESS.out.preprocessed_data)

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
