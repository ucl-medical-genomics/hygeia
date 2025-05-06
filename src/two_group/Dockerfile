FROM python:3.8.5-slim-buster

# procps  Needed by Nextflow
RUN apt-get update && apt-get install -y procps  build-essential && rm -rf /var/lib/apt/lists/*

WORKDIR /hygeia/two_group
COPY requirements.txt ./
RUN pip install -r requirements.txt

COPY preprocess_data.py ./
COPY run_inference_two_groups.py ./
COPY aggregate_results.py ./
COPY get_dmps.py ./
COPY multiple_testing.py ./
COPY get_chrom_segments.py ./

COPY ./bin ./bin/
COPY ./hygeia ./hygeia/
COPY ./hygeia.docker ./bin/hygeia

ENV PYTHONPATH=/hygeia/two_group/
ENV HYGEIA_TWO_GRP_BASE_PATH="/hygeia/two_group/"
ENV PATH="/hygeia/two_group/bin:${PATH}"
ENV HYGEIA_VERSION=0.1.14
