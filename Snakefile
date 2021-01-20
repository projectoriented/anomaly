include: "rules/common.smk"


# -------- Target Rules -------- #

rule all:
    params:
        job_name="localrule-target_all",
    input:
        out + "/multiqc_report.html",
        out + "/star_index",
        out + "/rsem_index",
        proj_dir + "/sample_annotation.tsv",


# -------- Setup Singularity -------- #
# singularity:


# -------- Load Rules -------- #

include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/mapping.smk"
include: "rules/counting.smk"
include: "rules/utils.smk"