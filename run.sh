nextflow run ucl-medical-genomics/hygeia \
  --cpg_file_path "s3://hygeia-analysis/data/ref/cpg.tsv.gz" \
  --sample_sheet "s3://hygeia-analysis/data/aging/sample_sheet.csv" \
  --output_dir "s3://hygeia-analysis/results/1" \
  -bucket-dir "s3://hygeia-analysis/tmp/work/1" \
  -c aws.config -resume -with-tower