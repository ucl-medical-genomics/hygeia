FROM python:3.8.5-slim-buster

# procps  Needed by Nextflow
RUN apt-get update && apt-get install -y procps  build-essential && rm -rf /var/lib/apt/lists/*

WORKDIR /hygeia/TwoGroup
COPY requirements.txt ./
RUN pip install -r requirements.txt

COPY preprocess_data.py ./
COPY run_inference_two_groups.py ./

COPY ./bin ./bin/
ENV PATH="/hygeia/TwoGroup/bin:${PATH}"
