localrules: generate_drop_sample_annot,

rule generate_drop_sample_annot:
    input:
        sample_tsv=proj_dir + "/samples.tsv",
        rna_meta_tsv=proj_dir + "/rna_meta.tsv",
        multiqc_report=out + "/multiqc_report.html",
    output:
        proj_dir + "/sample_annotation.tsv"
    script:
        "../scripts/generate_drop_sample_annot.py"