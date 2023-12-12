#!/usr/bin/env nextflow

params.cpg_file_path = "./tmp/test_data/cpg.tsv.gz"
params.case_data_path = "./tmp/test_data/EGAZ00001016574_90_cpg.txt.gz"
params.case_id_names = "EGAZ00001016574_90"
params.control_data_path = "./tmp/test_data/EGAZ00001016575_90_cpg.txt.gz"
params.control_id_names = "EGAZ00001016575_90"
params.output_path_preprocessed = "./results/preprocessed"
params.results_dir = "./results/out"
params.data_dir = "./results/preprocessed"
params.single_group_dir = "../SingleGroup/result/params"
params.chrom = 22

process preprocess {
  container 'hygiea/two_group'
    input:
    path cpg_file_path
    path case_data_path
    path control_data_path

    output:
    path "${params.output_path_preprocessed}", emit: processed_data_dir

    script:
    """
    hygeia preprocess \
      -cpg_file_path ${cpg_file_path} \
      -case_data_path ${case_data_path} \
      -case_id_names ${params.case_id_names} \
      -control_data_path ${control_data_path} \
      -control_id_names ${params.control_id_names} \
      -output_path ${params.output_path_preprocessed}
    """
}

process infer {
  container 'hygiea/two_group'
    input:
    path data_dir

    output:
    path "${params.results_dir}"

    script:
    """
    hygeia infer \
      -results_dir $params.results_dir \
      -data_dir $data_dir \
      -single_group_dir $params.single_group_dir \
      -chrom $params.chrom
    """
}

workflow {
    preprocess(
        file(params.cpg_file_path),
        file(params.case_data_path),
        file(params.control_data_path)
    )
    infer(preprocess.out.processed_data_dir)
}
