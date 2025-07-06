#!/usr/bin/env nextflow

process PREPROCESS {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.23'
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
          path("nextflow_output/positions_${chrom}.txt.gz", arity: 1),
          path("nextflow_output/n_total_reads_case_${chrom}.txt.gz", arity: 1),
          path("nextflow_output/n_total_reads_control_${chrom}.txt.gz", arity: 1),
          path("nextflow_output/n_methylated_reads_case_${chrom}.txt.gz", arity: 1),
          path("nextflow_output/n_methylated_reads_control_${chrom}.txt.gz", arity: 1),
          path("nextflow_output/cpg_sites_merged_${chrom}.txt.gz", arity: 1), emit: preprocessed_data
    path "nextflow_output/versions.yml"

    script:
    def caseIdArgs = case_ids.collect { "--case_id_names '$it'" }.join(" ")
    def caseFileArgs = case_files.collect { "--case_data_path '$it'" }.join(" ")
    def controlIdArgs = control_ids.collect { "--control_id_names '$it'" }.join(" ")
    def controlFileArgs = control_files.collect { "--control_data_path '$it'" }.join(" ")

    """
    hygeia preprocess ${caseIdArgs} ${caseFileArgs} ${controlIdArgs} \
        ${controlFileArgs} --cpg_file_path ${cpg_file_path} \
        --output_path nextflow_output --chromosome ${chrom}

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
