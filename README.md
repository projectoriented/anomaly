# anomaly
A bioinformatics pipeline (in development) to detect aberrant RNA events in rare disease for clinical diagnostics

anomaly is a pipeline that uses Snakemake as its workflow manager for modularity and scalability for end-to-end RNA-seq analyses. The goal is to increase clinical diagnostics in Rare Diseases with the additional layer of RNA-seq to uncover insights that WGS and WES couldn't alone. The workflow will undergo the standard bioinformatic RNA-seq procedures (e.g. QC, trim, align, mark duplicate, produce normalized counts) and finally utilize [DROP](https://github.com/gagneurlab/drop) to detect aberrant activity in the transcripts.

## Overview
![pipeline vector](https://github.com/projectoriented/anomaly/blob/main/images/dag.svg)

## Getting started
1. Clone the repo
2. Install Snakemake version >= 6.0.0
3. Build the container using Docker `docker build --tag local/ Dockerfile` or if Singularity is preferred: `singularity pull --name anomaly.sif docker://lettucerap/anomaly:latest` (default for Snakefile is singularity)
4. Have at least 10 samples per tissue (requirement from DROP) and use a public RNA-seq dataset with as close library preparation and tissue. This is for DROP to have enough data to perform statistical analyses and infer true biological variance. (need to write more about how the samples.tsv is written out)
5. Start the analysis!
    * Begin with a dry-run
    ```
    snakemake -np
    ```
    * If dry-run looks good, proceed with:
    ```
    snakemake --profile slurm --jobs [put # here] --use-singularity --singularity-args "--bind [paths that lives outside the container] --fakeroot"
    ```
     A slurm profile located in `~/.config/snakemake/slurm/config.yaml` is required for the `--profile` option. I'd recommend doing this on a separate screen using `screen -S rna_analysis` but otherwise `--immediate-submit --notemp --scheduler greedy` options can be added and use the slurm bash script.
