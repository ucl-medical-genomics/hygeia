# Hygeia

> TODO: Finish Readme

## Requirements

* Nextflow (https://www.nextflow.io/docs/latest/getstarted.html)
* Docker

## Single Group Analysis

To run the pipeline, after installing nextflow, you can run the following command. This will automatically download the latest version of the pipeline from Github.

### Run the pipeline with a config file (recommended)

```bash
nextflow run ucl-medical-genomics/hygeia -entry SingleGroup/main.nf -c nextflow.config
```


> Note you can use -r <git_branch_name> to use a different github branch.

You would need to create a config file with the following paramaters defined. You can also overwrite any parameters in the default nextflow config file.

```bash
params {
    base_output_dir = "./result"
    params_data_dir = "params"
    simulated_data_dir = "simulated_data"
    regimes_dir = "regimes"

    input_mu_csv_path = null
    input_sigma_csv_path = null
    input_omega_csv_path = null
    input_kappa_csv_path = null
    input_u_csv_path = null
    input_p_csv_path = null

    regimes_csv_path = null
    n_methylated_reads_csv = null
    genomic_positions_csv = null
    n_total_reads_csv = null
    regime_probabilities_csv = null
    theta_trace_csv = null
}
```

Run the pipeline without a config file
All the paramaters in the config file can be set via the CLI. This may be useful for scripts.

```bash
nextflow run ucl-medical-genomics/hygeia \
  --base_output_dir /path/to/output/ \
  --params_data_dir params \
  --simulated_data_dir simulated_data \
  --regimes_dir regimes / \
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
docker build -t hygeia/single_group src/single_group
docker build -t hygeia/two_group src/two_group

docker tag hygeia/single_group ghcr.io/ucl-medical-genomics/hygeia_single_group:v0.0.1
docker tag hygeia/single_group ghcr.io/ucl-medical-genomics/hygeia_single_group:latest

docker tag hygeia/two_group ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.0.1
docker tag hygeia/two_group ghcr.io/ucl-medical-genomics/hygeia_two_group:latest

docker push ghcr.io/ucl-medical-genomics/hygeia_single_group:latest
docker push ghcr.io/ucl-medical-genomics/hygeia_single_group:v0.0.1
docker push ghcr.io/ucl-medical-genomics/hygeia_two_group:latest
docker push ghcr.io/ucl-medical-genomics/hygeia_two_group:v0.0.1
```

# Tutorial - Run Single Group analysis with NA12878

1. Download data

```bash
curl -LO https://www.encodeproject.org/files/ENCFF608CXC/@@download/ENCFF608CXC.bigWig

https://www.encodeproject.org/files/ENCFF446HUA/@@download/ENCFF446HUA.bed.gz
```

2. Run Hygeia

```bash
nextflow run ....

```

