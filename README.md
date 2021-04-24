# anomaly
A pipeline to detect aberrant RNA events in rare disease for clinical diagnostics

anomaly is a pipeline that uses Snakemake as its workflow manager for comprehensibility and scalability for end-to-end RNA-seq analyses. The goal is to increase clinical diagnostics in Rare Diseases with the additional layer of RNA-seq to uncover insights that WGS and WES couldn't alone. The workflow will undergo the standard bioinformatic RNA-seq procedures (e.g. QC, trim, align, mark duplicate, produce normalized counts) and finally utilize [DROP](https://github.com/gagneurlab/drop) to detect aberrant activity in the transcripts.

## Overview
![pipeline vector](https://github.com/projectoriented/anomaly/blob/main/images/dag.svg)

## Getting started
#### Prerequisites
