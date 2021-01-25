include: "rules/common.smk"
localrules: done,

# -------- Target Rules -------- #
rule done:
    params:
        tmp=tmp
    input:
        out + "/multiqc_report.html",
        proj_dir + "/sample_annotation.tsv",
    shell:
        "rm -rf {params.tmp}/tmp.*"


# -------- Setup Singularity -------- #
# singularity:


# -------- Load Rules -------- #

include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/mapping.smk"
include: "rules/counting.smk"
include: "rules/utils.smk"