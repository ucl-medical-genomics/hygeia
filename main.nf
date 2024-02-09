#!/usr/bin/env nextflow

params.cpg_file_path = "/home/imoghul/d/hygeia/run_08_02_2024/data/ref/cpg.tsv.gz"
params.sample_sheet = "/home/imoghul/d/hygeia/run_08_02_2024/data/t1d/sample_sheet.csv"
params.output_path = '/home/imoghul/d/hygeia/run_08_02_2024/data/t1d/preprocessed2'

Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true, sep: ',', strip: true)
    .map { row -> tuple(row.group, row.id, row.file) }
    .groupTuple()
    .collect()
    .set { samples }

process preprocess {
    container 'ucl-medical-genomics/hygeia_two_group'

    input:
    tuple val(case_group), val(case_ids), path(case_files), val(control_group), val(control_ids), path(control_files)
    path cpg_file_path
    path output_path

    script:
    def caseIdArgs = case_ids.collect { "--case_id_names '$it'" }.join(" ")
    def caseFileArgs = case_files.collect { "--case_data_path '$it'" }.join(" ")
    def controlIdArgs = control_ids.collect { "--control_id_names '$it'" }.join(" ")
    def controlFileArgs = control_files.collect { "--control_data_path '$it'" }.join(" ")

    """
    hygeia preprocess \\
        ${caseIdArgs} \\
        ${caseFileArgs} \\
        ${controlIdArgs} \\
        ${controlFileArgs} \\
        --cpg_file_path ${cpg_file_path} \\
        --output_path ${output_path} \\
        --chrom 22
    """
}

workflow {
    preprocess(samples, params.cpg_file_path, params.output_path)
}


// #!/usr/bin/env nextflow

// params.cpg_file_path = "/home/imoghul/d/hygeia/run_08_02_2024/data/t1d/ref/cpg.tsv.gz"
// params.sample_sheet = "/home/imoghul/d/hygeia/run_08_02_2024/data/t1d/sample_sheet.csv"
// params.output_path = '/home/imoghul/d/hygeia/hygeia/result1'
// Channel
//     .fromPath(params.sample_sheet)
//     .splitCsv(header: true, sep: ',', strip: true)
//     .map { row -> tuple(row.group, row.id, row.file) }
//     .groupTuple()
//     .collect()
//     .set { samples }

// process preprocess {
//     container 'ucl-medical-genomics/hygeia_two_group'

//     input:
//     tuple (val(group), val(id), path(file))
//     path cpg_file_path
//     path output_path

//     script:


//     """
//     echo hygeia preprocess \\
//         --cpg_file_path ${cpg_file_path} \\
//         ${args} \\
//         --output_path ${output_path}
//         exit 1
//     """
// }

// workflow {
//     preprocess(samples, params.cpg_file_path, params.output_path)
// }



// #!/usr/bin/env nextflow
// params.cpg_file_path = "/home/imoghul/d/hygeia/run_08_02_2024/data/t1d/ref/cpg.tsv.gz"
// params.sample_sheet = "/home/imoghul/d/hygeia/run_08_02_2024/data/t1d/sample_sheet.csv"
// params.output_path = './results' // Default output path

// Channel
//     .fromPath(params.sample_sheet)
//     .splitCsv(header: true, sep: ',', strip: true)
//     .set { samples }

// process preprocess {
//     container 'ucl-medical-genomics/hygeia_two_group'
//     input:
//     path cpg_file_path
//     path output_path
//     val samples

//     script:
//     // Adjusted to dynamically construct the command with proper handling for each sample
//     def args = samples.collect { sample ->
//         return "--${sample.group}_id_names ${sample.id} --${sample.group}_data_path ${sample.file}"
//     }.join(' ')

//     """
//     echo hygeia preprocess \
//         -cpg_file_path ${cpg_file_path} \
//         ${args.join(' ')} \
//         -output_path ${params.output_path}
//       exit 1
//     """
// }

// workflow {
//     preprocess(
//         file(params.cpg_file_path),
//         file(params.output_path),
//         samples
//     )
// }




// // // params.cpg_file_path = "./tmp/test_data/cpg.tsv.gz"
// // // params.case_data_path = "./tmp/test_data/EGAZ00001016574_90_cpg.txt.gz"
// // // params.case_id_names = "EGAZ00001016574_90"
// // // params.control_data_path = "./tmp/test_data/EGAZ00001016575_90_cpg.txt.gz"
// // // params.control_id_names = "EGAZ00001016575_90"
// // // params.output_path_preprocessed = "./results/preprocessed"
// // // params.results_dir = "./results/out"
// // // params.data_dir = "./results/preprocessed"
// // // params.single_group_dir = "../SingleGroup/result/params"
// // params.chrom = 22

// // process preprocess {
// //   container 'ucl-medical-genomics/hygeia_two_group'
// //     input:
// //     path cpg_file_path
// //     path case_data_path
// //     path control_data_path

// //     output:
// //     path "${params.output_path_preprocessed}", emit: processed_data_dir

// //     script:
// //     """
// //     hygeia preprocess \
// //       -cpg_file_path ${cpg_file_path} \
// //       -case_data_path ${case_data_path} \
// //       -case_id_names ${params.case_id_names} \
// //       -control_data_path ${control_data_path} \
// //       -control_id_names ${params.control_id_names} \
// //       -output_path ${params.output_path_preprocessed}
// //     """
// // }

// // process infer {
// //   container 'ucl-medical-genomics/hygeia_two_group'
// //     input:
// //     path data_dir

// //     output:
// //     path "${params.results_dir}"

// //     script:
// //     """
// //     hygeia infer \
// //       -results_dir $params.results_dir \
// //       -data_dir $data_dir \
// //       -single_group_dir $params.single_group_dir \
// //       -chrom $params.chrom
// //     """
// // }

// // workflow {
// //     preprocess(
// //         file(params.cpg_file_path),
// //         file(params.case_data_path),
// //         file(params.control_data_path)
// //     )
// //     infer(preprocess.out.processed_data_dir)
// // }
