container: "docker://quay.io/biocontainers/drop:1.0.5--pyhdfd78af_0"

rule remove_alt_chr_for_drop:
    input:
        out + "/{sample}/picard_markdupe/{sample}_trim_star_marked.bam",
    output:
        proj_dir + "/drop_files/{sample}/{sample}_trim_star_marked_no_alt.bam",
    group: "group_drop"
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
        rna_meta_tsv=config['rna_meta'],
        fexist=expand([proj_dir + "/drop_files/{s.sample}/{s.sample}_trim_star_marked_no_alt.bam"], s=samples.itertuples()),
    output:
        drop_proj + "/sample_annotation.tsv",
    resources:
        cores=config["cores"]["trimming"],
        time_min=config["time_min"]["default"],
    group: "group_drop"
    script:
        "../scripts/generate_drop_sample_annot.py"

rule generate_drop_config:
    input:
        drop_proj + "/sample_annotation.tsv",
    output:
        drop_proj + "/config.yaml",
    resources:
        cores=config["cores"]["trimming"],
        time_min=config["time_min"]["default"],
    group: "group_drop"
    script:
        "../scripts/modify_drop_config.py"

rule drop_init:
    input:
        sa=drop_proj + "/sample_annotation.tsv",
        config=drop_proj + "/config.yaml",
    params:
        dir=drop_proj,
    log:
        temp(drop_proj + "/done.txt")
    resources:
        cores=config["cores"]["trimming"],
        time_min=config["time_min"]["default"],
    group: "group_drop"
    shell:
        """        
        cd {params.dir} &&
        drop init > {log}
        """

rule run_drop_AE:
    input:
        initialized=drop_proj + "/done.txt",
    output:
        html=drop_proj + "/output/html/AberrantExpression/Outrider/v34/Summary_blood.html",
    params:
        dir=drop_proj,
    resources:
        cores=config["cores"]["drop"],
        mem=config["mem"]["drop"],
        time_min=config["time_min"]["drop"],
    log:
        drop_proj + "/outAE.log"
    shell:
        """
        cd {params.dir} &&         
        snakemake aberrantSplicing --cores {resources.cores} >& {log}
        """

rule run_drop_AS:
    input:
        initialized=drop_proj + "/done.txt",
    output:
        html=drop_proj + "/output/html/AberrantSplicing/blood--v34_summary.html",
    params:
        dir=drop_proj,
    resources:
        cores=config["cores"]["drop"],
        mem=config["mem"]["drop"],
        time_min=config["time_min"]["drop"],
    log:
        drop_proj + "/outAS.log"
    shell:
        """
        cd {params.dir} &&         
        snakemake aberrantSplicing --cores {resources.cores} >& {log}
        """