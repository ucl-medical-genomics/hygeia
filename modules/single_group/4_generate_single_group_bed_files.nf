#!/usr/bin/env nextflow

process GENERATE_SINGLE_GROUP_BED_FILES {
    tag "${case_id}_${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_single_group:v0.1.27'
    publishDir "${params.output_dir}/4_SINGLE_GROUP_OUTPUT", mode: 'copy',
        pattern: 'nextflow_output/*',
        saveAs: { fn -> fn.replace("nextflow_output", "./") }

    memory { 4.GB * task.attempt }

    time { 5.minutes + (10.minutes * task.attempt) }

    input:
    tuple val(case_id),
          val(chrom),
          path(regimes_csv)

    output:
    path "nextflow_output/*"

    script:
    """
    hygeia make_bed_file \
        --chr ${chrom} \
        --regimes_file ${regimes_csv} \
        --output_file nextflow_output/${case_id}_regimes_${chrom}.bed

    bgzip nextflow_output/${case_id}_regimes_${chrom}.bed
    tabix -p bed nextflow_output/${case_id}_regimes_${chrom}.bed.gz

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: \$(hygeia --version | sed 's/hygeia version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p nextflow_output
    touch nextflow_output/${case_id}_regimes_${chrom}.csv.gz

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}
