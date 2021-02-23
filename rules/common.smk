import pandas as pd


# -------- Config files -------- #
configfile: "config/config.yaml"
configfile: "config/cluster.yaml"


# -------- Global variables -------- #
out = config["common_out"]
tmp = config["tmp"]
script_dir = config["script_dir"]
proj_dir = config["proj_dir"]


# --------  Load sample sheet -------- #
samples = pd.read_table(config["samples"], dtype=str).set_index(["sample", "lane", "sample_number"], drop=False)
# samples = pd.read_table("/home/mei.wu/rna-seq/samples.tsv", dtype=str).set_index(["sample", "lane", "sample_number"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index
samples = samples.sort_index() # increase downstream performance


# -------- Wildcard constraints -------- #
wildcard_constraints:
    sample="|".join(samples["sample"]),
    lane="|".join(samples["lane"]),
    sample_number="|".join(samples["sample_number"]),


# -------- Helper functions -------- #
def get_fastq(wildcards):
    """Get fastq files with given wildcards."""
    return samples.loc[(wildcards.sample, wildcards.lane, wildcards.sample_number),
                       ["fq1", "fq2"]].dropna()

def get_trimmed_read1(wildcards):
    """Get trimmed read 1s with given wildcards."""
    lanes = samples.loc[wildcards.sample]["lane"].unique()
    sn = samples.loc[wildcards.sample]["sample_number"].unique()
    return expand([out + "/{sample}/trim-galore/{l}_{sample}_{sample_number}_val_1.fq.gz"], l=lanes, sample_number=sn, **wildcards)

def get_trimmed_read2(wildcards):
    """Get trimmed read 2s with given wildcards."""
    lanes = samples.loc[wildcards.sample]["lane"].unique()
    sn = samples.loc[wildcards.sample]["sample_number"].unique()
    return expand([out + "/{sample}/trim-galore/{l}_{sample}_{sample_number}_val_2.fq.gz"], l=lanes, sample_number=sn, **wildcards)

def get_read_group(wildcards):
    """Denote ID, LB, PL, PU, SM in read group."""
    sample = wildcards.sample
    all = samples.loc[sample, ["id", "sample", "platform", "pu", "sample"]].values
    all = ["ID:{} LB:{} PL:{} PU:{} SM:{}".format(*x) for x in all]
    return [all[0]] + [", " + s for s in all[1:]]