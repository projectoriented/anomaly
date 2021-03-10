include: "rules/common.smk"
localrules: done,

# -------- Target Rules -------- #
rule done:
    input:
        out + "/multiqc_report.html",
        proj_dir + "/sample_annotation.tsv",
    params:
        tmp=tmp,
    shell:
        """        
        rm -rf ~/slurm/logs/*;         
        find {params.tmp} -mindepth 1 -name tmp* -exec rm -r {{}} +;        
        """


# -------- Setup Singularity -------- #
# singularity:


# -------- Load Rules -------- #

include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/mapping.smk"
include: "rules/counting.smk"
include: "rules/utils.smk"