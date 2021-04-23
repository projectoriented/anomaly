FROM continuumio/miniconda3

LABEL description="anomaly- an RNA-seq pipeline to detect anomalies in the transcriptome"

ENV PATH="/opt/conda/bin/:${PATH}"

COPY environment.yaml /
RUN conda env update --name root --file /environment.yaml && conda clean --all --yes