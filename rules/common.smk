import pandas as pd


# -------- Config file and sample sheets --------#
configfile: "config.yaml"

samples = pd.read_table(config["samples"], dtype=str).set_index(["sample", "lane", "sample_number"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index


# -------- Global variables --------#
out = config["common_out"]
tmp = config["tmp"]
script_dir = config["script_dir"]


# -------- Wildcard constraints --------#
wildcard_constraints:
    sample="|".join(samples["sample"]),
    lane="|".join(samples["lane"]),
    sample_number="|".join(samples["sample_number"]),


# -------- Functions --------#
def get_fastq(wildcards):
    """Get fastq files with given wildcards."""
    return samples.loc[(wildcards.sample, wildcards.lane, wildcards.sample_number),
                       ["fq1", "fq2"]].dropna()

def get_trimmed_reads(wildcards):
    """Get PE trimmed reads with given wildcards."""
    return expand(out + "/{sample}/trim-galore/{lane}_{sample}_{sample_number}_val_{group}.fq.gz",
                  group=[1, 2], **wildcards)

def get_read_group(wildcards):
    """Denote ID, LB, PL, PU, SM in read group."""
    sample=wildcards.sample
    platform=samples.loc[(wildcards.sample, wildcards.lane, wildcards.sample_number),"platform"]
    identifier=samples.loc[(wildcards.sample, wildcards.lane, wildcards.sample_number),"id"]
    pu=samples.loc[(wildcards.sample, wildcards.lane, wildcards.sample_number),"pu"]
    return f"ID:{identifier} LB:{sample} PL:{platform} PU:{pu} SM:{sample}"