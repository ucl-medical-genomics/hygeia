# Hygeia

> [!note]
>
> Hygeia is still in Alpha. However, feel free to use the app on your data and
> reach out to us if you need any help.
>

## Requirements

* Nextflow [see here](https://www.nextflow.io/docs/latest/getstarted.html)
* [Docker](https://www.docker.com/) / [Singularity](https://docs.sylabs.io/guides/3.7/user-guide/)

## Two Group Analysis

To run the pipeline, after installing nextflow, you can run the following command. This will automatically download the latest version of the pipeline from Github.

### Run the pipeline with a config file (recommended)

```bash
nextflow run ucl-medical-genomics/hygeia -c run.config
```

> Note you can use -r <git_branch_name> to use a different github branch.

You would need to create a config file with the following paramaters defined. You can also overwrite any parameters in the default nextflow config file.

```bash
params.cpg_file_path = "/scratch/imoghul/hygeia_data/ref/cpg.tsv.gz"
params.sample_sheet = "/scratch/imoghul/hygeia_data/aging/sample_sheet.csv"
params.output_dir = "results"
params.meteor_mu = "0.95,0.05,0.8,0.2,0.50,0.50"
params.meteor_sigma = "0.05,0.05,0.1,0.1,0.1,0.2886751"
params.min_cpg_sites_between_change_points = 3
params.num_of_inference_seeds = 2
```

Run the pipeline without a config file
All the paramaters in the config file can be set via the CLI. This may be useful for scripts.

```bash
nextflow run ucl-medical-genomics/hygeia \
  --cpg_file_path "/scratch/imoghul/hygeia_data/ref/cpg.tsv.gz" \
  --sample_sheet "/scratch/imoghul/hygeia_data/aging/sample_sheet.csv" \
  --output_dir "results" \
  --meteor_mu "0.95,0.05,0.8,0.2,0.50,0.50" \
  --meteor_sigma "0.05,0.05,0.1,0.1,0.1,0.2886751" \
  --min_cpg_sites_between_change_points 3 \
  --num_of_inference_seeds 2 \
  -with-report report.html -with-dag flowchart.pdf
```

## Development

After running the above commands, the pipeline will be cloned into `~/.nextflow/assets/ucl-medical-genomics/hygeia`. If you make any changes here, you will need to commit them before running the above command again.

You may want to add `process.errorStrategy = 'terminate'` to a local nextflow config to override the default (which is to ignore errors).

### Release a new Docker build

The pipeline uses docker images on Github Docker Registry. If you make any changes to the underlying files including in the dockerfile, please push them to Dockerhub:

1. Login using a token with access to Github Packages. See [here](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry) for more info.

```bash
export CR_PAT=YOUR_TOKEN

echo $CR_PAT | docker login ghcr.io -u USERNAME --password-stdin
```

2. Build docker image and upload to Github Packages

```bash
VERSION=v0.1.28
IMAGES=(
  "hygeia_single_group:src/single_group"
  "hygeia_two_group:src/two_group"
)

for img in "${IMAGES[@]}"; do
    d_label=${img%%:*}
    d_src_path=${img#*:}
    docker build -t ${d_label} ${d_src_path}
    docker tag ${d_label} ghcr.io/ucl-medical-genomics/${d_label}:${VERSION}
    docker tag ${d_label} ghcr.io/ucl-medical-genomics/${d_label}:latest
    docker push ghcr.io/ucl-medical-genomics/${d_label}:${VERSION}
    docker push ghcr.io/ucl-medical-genomics/${d_label}:latest
done
```

# Tutorial - Run Single Group analysis with NA12878

> TODO: Complete Tutorial

1. Download data

```bash
curl -LO https://www.encodeproject.org/files/ENCFF608CXC/@@download/ENCFF608CXC.bigWig

https://www.encodeproject.org/files/ENCFF446HUA/@@download/ENCFF446HUA.bed.gz
```

2. Run Hygeia

```bash
nextflow run ....

```
