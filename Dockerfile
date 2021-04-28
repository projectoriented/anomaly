FROM continuumio/miniconda3

LABEL description="anomaly- an RNA-seq pipeline to detect anomalies in the transcriptome"

ENV PATH="/opt/conda/bin/:${PATH}"

COPY environment.yaml /
RUN conda install --name base -c conda-forge mamba && \
    mamba env update --name base --file /environment.yaml && \
    conda clean --all --yes
