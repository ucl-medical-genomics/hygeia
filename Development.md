# Local Development

## Requirements

The following dependencies are required

* Docker

```bash
# Single Group
docker build -t ucl-medical-genomics/hygeia_single_group src/single_group

# Two Group
docker build -t ucl-medical-genomics/hygeia_two_group src/two_group
```

## Run locally

```bash
nextflow run single_group.nf -
```
