#!/usr/bin/env nextflow

process SINGLE_GRP_PREPROCESS {
    tag "${case_id}_${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.27'
    publishDir "${params.output_dir}/1_PREPROCESS", mode: 'copy',
        pattern: 'nextflow_output/*',
        saveAs: { fn -> fn.replace("nextflow_output", "./") }

    memory { 16.GB + (4.GB * task.attempt) }

    input:
    tuple val(case_id), path(case_file), val(chrom)
    path cpg_file_path

    output:
    tuple val(case_id),
          val(chrom),
          path("nextflow_output/${case_id}_positions_${chrom}.txt.gz", arity: 1),
          path("nextflow_output/${case_id}_n_total_reads_case_${chrom}.txt.gz", arity: 1),
          path("nextflow_output/${case_id}_n_methylated_reads_case_${chrom}.txt.gz", arity: 1),
          path("nextflow_output/${case_id}_cpg_sites_merged_${chrom}.txt.gz", arity: 1), emit: preprocessed_data
    path "nextflow_output/versions.yml"

    script:
    """
    hygeia preprocess --case_id_names '${case_id}' --case_data_path '${case_file}' \
        --cpg_file_path ${cpg_file_path} --output_path nextflow_output --chromosome ${chrom}

    mv nextflow_output/positions_${chrom}.txt.gz nextflow_output/${case_id}_positions_${chrom}.txt.gz
    mv nextflow_output/n_total_reads_case_${chrom}.txt.gz nextflow_output/${case_id}_n_total_reads_case_${chrom}.txt.gz
    mv nextflow_output/n_total_reads_control_${chrom}.txt.gz nextflow_output/${case_id}_n_total_reads_control_${chrom}.txt.gz
    mv nextflow_output/n_methylated_reads_case_${chrom}.txt.gz nextflow_output/${case_id}_n_methylated_reads_case_${chrom}.txt.gz
    mv nextflow_output/n_methylated_reads_control_${chrom}.txt.gz nextflow_output/${case_id}_n_methylated_reads_control_${chrom}.txt.gz
    mv nextflow_output/cpg_sites_merged_${chrom}.txt.gz nextflow_output/${case_id}_cpg_sites_merged_${chrom}.txt.gz

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p nextflow_output
    touch nextflow_output/${case_id}_positions_${chrom}.txt.gz
    touch nextflow_output/${case_id}_n_total_reads_case_${chrom}.txt.gz
    touch nextflow_output/${case_id}_n_total_reads_control_${chrom}.txt.gz
    touch nextflow_output/${case_id}_n_methylated_reads_case_${chrom}.txt.gz
    touch nextflow_output/${case_id}_n_methylated_reads_control_${chrom}.txt.gz
    touch nextflow_output/${case_id}_cpg_sites_merged_${chrom}.txt.gz

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}
