localrules: generate_drop_sample_annot,

rule remove_alt_chr_for_drop:
    input: out + "/{sample}/picard_markdupe/{lane}_{sample}_{sample_number}_trim_star_marked.bam",
    output: out + "/{sample}/picard_markdupe/{lane}_{sample}_{sample_number}_trim_star_marked_alt_dropped.bam",
    params:
        tmp_dir = tmp,
        job_name = "alt4removal_{lane}_{sample}_{sample_number}",
    resources:
        cores = config["cores"]["default"],
        mem = config["mem"]["default"],
        time_min = config["time_min"]["default"],
    shell:
        """
        chr_list=$(samtools view {input} | cut -f3 | uniq | grep -E '^.{{1,5}}$') &&
        samtools view -h -b ${{chr_list}} > {output} &&
        samtools index {output}
        """


rule generate_drop_sample_annot:
    input:
        sample_tsv=proj_dir + "/samples.tsv",
        rna_meta_tsv=proj_dir + "/rna_meta.tsv",
        multiqc_report=out + "/multiqc_report.html",
    output:
        proj_dir + "/sample_annotation.tsv"
    script:
        "../scripts/generate_drop_sample_annot.py"