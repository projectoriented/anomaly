rule generate_drop_sample_annot:
    input:
        sample_tsv=proj_dir + "/samples.tsv",
        rna_meta_tsv=proj_dir + "/rna_meta.tsv",
        multiqc_report=out + "/multiqc_report.html",
    output:
        proj_dir + "/sample_annotation.tsv"
    params:
        job_name="snakemake_post_processing",
    resources:
        cores=config["cores"]["default"],
        mem=config["mem"]["default"],
        time_min=config["time_min"]["default"],
    script:
        "../scripts/generate_drop_sample_annot.py"