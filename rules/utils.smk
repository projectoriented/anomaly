localrules: generate_drop_sample_files,


rule remove_alt_chr_for_drop:
    input:
        out + "/{sample}/picard_markdupe/{lane}_{sample}_{sample_number}_trim_star_marked.bam",
    output:
        proj_dir + "/drop_files/{sample}/{lane}_{sample}_{sample_number}_trim_star_marked_no_alt.bam",
    params:
        job_name="{lane}_{sample}_{sample_number}_chr_alt4removal",
    resources:
        cores=config["cores"]["default"],
        time_min=25,
    shell:
        """
        chr_list=$(samtools view {input} | cut -f3 | uniq | grep -E \"^.{{1,5}}$\") && 
        samtools view -b -h {input} $chr_list -o {output} && 
        samtools index -b {output}
        """

rule generate_drop_sample_files:
    input:
        sample_tsv=proj_dir + "/samples.tsv",
        rna_meta_tsv=proj_dir + "/rna_meta.tsv",
        fexist=expand([proj_dir + "/drop_files/{s.sample}/{s.lane}_{s.sample}_{s.sample_number}_trim_star_marked_no_alt.bam"], s=samples.itertuples()),
    output:
        proj_dir + "/sample_annotation.tsv",
    script:
        "../scripts/generate_drop_sample_annot.py"
