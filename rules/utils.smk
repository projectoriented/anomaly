localrules: generate_drop_sample_files,

rule remove_alt_chr_for_drop:
    input:
        out + "/{sample}/picard_markdupe/{sample}_trim_star_marked.bam",
    output:
        proj_dir + "/drop_files/{sample}/{sample}_trim_star_marked_no_alt.bam",
    params:
        job_name="{sample}_chr_alt4removal",
    resources:
        cores=config["cores"]["trimming"],
        time_min=config["time_min"]["default"],
    shell:
        """
        chr_list=$(samtools view {input} | cut -f3 | uniq | grep -E \"^.{{1,5}}$\") && 
        samtools view -@ {resources.cores} -b -h {input} $chr_list -o {output} && 
        samtools index -@ {resources.cores} -b {output}
        """

rule generate_drop_sample_files:
    input:
        rna_meta_tsv=proj_dir + "/merged_rna_meta.txt2",
        fexist=expand([proj_dir + "/drop_files/{s.sample}/{s.sample}_trim_star_marked_no_alt.bam"], s=samples.itertuples()),
    output:
        proj_dir + "/sample_annotation.tsv",
    script:
        "../scripts/generate_drop_sample_annot.py"

rule generate_drop_config:
    input:
        proj_dir + "/sample_annotation.tsv",
    output:
        proj_dir + "/drop/config.yaml"
    script:
        "../scripts/modify_drop_config.py"
