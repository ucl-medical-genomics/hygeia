#!/usr/bin/env nextflow

process PREPROCESS {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.27'
    publishDir "${params.output_dir}/1_PREPROCESS", mode: 'copy',
        pattern: 'nextflow_output/*',
        saveAs: { fn -> fn.replace("nextflow_output", "./") }

    memory { 16.GB + (4.GB * task.attempt) }

    input:
    tuple val(chrom),
          val(sample_name),
          path(sample_file)

    output:
    tuple val(chrom),
          path("nextflow_output/positions_${chrom}.txt.gz", arity: 1),
          path("nextflow_output/n_total_reads_case_${chrom}.txt.gz", arity: 1),
          path("nextflow_output/n_methylated_reads_case_${chrom}.txt.gz", arity: 1),
          path("nextflow_output/cpg_sites_merged_${chrom}.txt.gz", arity: 1), emit: preprocessed_data
    path "nextflow_output/versions.yml"

    script:
    """
    mkdir -p nextflow_output
    tail -n +2 ${sample_file} | cut -f 1 | gzip > nextflow_output/positions_${chrom}.txt.gz
    tail -n +2 ${sample_file} | cut -f 2 | gzip > nextflow_output/n_total_reads_case_${chrom}.txt.gz
    tail -n +2 ${sample_file} | cut -f 3 | gzip > nextflow_output/n_methylated_reads_case_${chrom}.txt.gz
    tail -n +2 ${sample_file} | wc -l | gzip > nextflow_output/cpg_sites_merged_${chrom}.txt.gz

    touch nextflow_output/versions.yml
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
