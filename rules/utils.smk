# rule generate_sample_list:
#     output: proj_dir + "/samples.tsv"
#     params:
#         script=script_dir + "/generate_sample_list.py",
#         arg=config["common_fastq_dir"]
#     shell:
#         "{params.script} {params.arg}"

rule generate_drop_sample_annot:
    input:
        sample_tsv=proj_dir + "/samples.tsv",
        multiqc_report=out + "/multiqc_report.html",
    output: proj_dir + "/sample_annotation.tsv"
    params:
        script=script_dir + "/generate_drop_sample_annot.py"
    script:
        "../scripts/generate_drop_sample_annot.py"