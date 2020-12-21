import pandas as pd

# -------- Config file and sample sheets --------#
configfile: "config.yaml"

samples = pd.read_table(config["samples"], dtype=str).set_index(["sample", "lane", "sample_number"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index


# -------- Wildcard constraints --------#
wildcard_constraints:
    sample="|".join(samples["sample"]),
    lane="|".join(samples["lane"]),
    sample_number="|".join(samples["sample_number"])

def get_fastq(wildcards):
    """Get fastq files of given wildcards."""
    return samples.loc[(wildcards.sample, wildcards.lane, wildcards.sample_number), ["fq1", "fq2"]].dropna()