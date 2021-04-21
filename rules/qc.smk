rule fastqc:
    input:
        get_fastq
    output:
        html=out + "/{sample}/fastqc/{sample}_{pu}_1_fastqc.html",
        zip=out + "/{sample}/fastqc/{sample}_{pu}_1_fastqc.zip",
    params:
        outpath=out + "/{sample}/fastqc",
        reg="{sample}.{pu}",
        tmp_dir=tmp,
        job_name="fastqc_{sample}_{pu}",
        script_dir=script_dir,
    resources:
        cores=config["cores"]["default"],
        time_min=config["time_min"]["default"],
    shell:
        """
        tmpdir=$(mktemp --directory {params.tmp_dir}/tmp.XXXXX) && 
        fastqc --noextract --dir ${{tmpdir}} --outdir {params.outpath} {input} && 
        {params.script_dir}/rename_files.sh {params.outpath} fastqc {params.reg}
        """


rule multiqc:
    input:
        expand([out + "/{s.sample}/fastqc/{s.sample}_{s.pu}_1_fastqc.zip"], s=samples.itertuples()),
        expand([out + "/{s.sample}/trim-galore/{s.sample}_{s.pu}_1.fastq.gz_trimming_report.txt"], s=samples.itertuples()),
        expand([out + "/{s.sample}/star_aln/{s.sample}_trim_star.Log.final.out"],s=samples.itertuples()),
        expand([out + "/{s.sample}/picard_markdupe/{s.sample}_trim_star_marked.metrics.txt"], s=samples.itertuples()),
        expand([out + "/{s.sample}/rsem/{s.sample}.stat/{s.sample}.cnt"], s=samples.itertuples()),
    output:
        report(
            out + "/multiqc_report.html",
            category="Quality control",
        ),
    params:
        common_dir=out,
        multiqc_config=proj_dir + "/config/multiqc_config.yaml",
        job_name="multiqc",
    resources:
        cores=config["cores"]["default"],
        time_min=config["time_min"]["default"],
    shell:
         "multiqc -f {params.common_dir} --filename {output} --config {params.multiqc_config}"