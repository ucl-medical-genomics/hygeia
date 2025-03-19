nextflow run ucl-medical-genomics/hygeia \
  -r main \
  --cpg_file_path "s3://hygeia-analysis/data/ref/cpg.tsv.gz" \
  --sample_sheet "s3://hygeia-analysis/data/aging/sample_sheet.csv" \
  --output_dir "s3://hygeia-analysis/results/1" \
  -bucket-dir "s3://hygeia-analysis/tmp/work/1" \
  -c aws.config -resume -with-tower

# Beck Server
nextflow run main.nf \
  --cpg_file_path "/scratch/imoghul/hygeia_data/ref/cpg.tsv.gz" \
  --sample_sheet "/scratch/imoghul/hygeia_data/aging/old/sample_sheet.csv" \
  -resume -with-tower -c nextflow.config --output_dir "out"
