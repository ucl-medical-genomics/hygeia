#!/usr/bin/env nextflow

process GET_DMPS {
    tag "${chrom}"
    container 'ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.1.19'
    publishDir "${params.output_dir}/6_GET_DMPS/${chrom}", mode: 'copy', 
        pattern: 'nextflow_output/*',
        saveAs: { fn -> fn.replace("nextflow_output", "./") }

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
    cp -r dmps_${chrom} nextflow_output

    cat <<-END_VERSIONS > nextflow_output/versions.yml
    "${task.process}":
        hygeia: stub
    END_VERSIONS
    """
}
