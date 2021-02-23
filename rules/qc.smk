rule fastqc:
    input:
        get_fastq
    output:
        html=out + "/{sample}/fastqc/{lane}_{sample}_{sample_number}_1_fastqc.html",
        zip=out + "/{sample}/fastqc/{lane}_{sample}_{sample_number}_1_fastqc.zip",
    params:
        tmp_dir=tmp,
        outpath=out + "/{sample}/fastqc",
        job_name="fastqc_{lane}_{sample}_{sample_number}",
        capture_group="{sample}_{sample_number}",
        script_dir=script_dir,
    resources:
        cores=config["cores"]["default"],
        time_min=config["time_min"]["default"],
    shell:
        """
        tmpdir=$(mktemp --directory {params.tmp_dir}/tmp.XXXXX) && 
        fastqc --noextract --dir ${{tmpdir}} --outdir {params.outpath} {input} && 
        {params.script_dir}/rename_files.sh {params.outpath} fastqc {wildcards.lane} {params.capture_group}
        """


rule multiqc:
    input:
        expand([out + "/{s.sample}/fastqc/{s.lane}_{s.sample}_{s.sample_number}_1_fastqc.zip"], s=samples.itertuples()),
        expand([out + "/{s.sample}/trim-galore/{s.lane}_{s.sample}_{s.sample_number}_1.fastq.gz_trimming_report.txt"], s=samples.itertuples()),
        expand([out + "/{s.sample}/star_aln/{s.sample}_trim_star.Log.final.out"],s=samples.itertuples()),
        expand([out + "/{s.sample}/picard_markdupe/{s.sample}_trim_star_marked.metrics.txt"], s=samples.itertuples()),
        expand([out + "/{s.sample}/rsem/{s.sample}.stat/{s.sample}.cnt"], s=samples.itertuples()),
    output:
        report(out + "/multiqc_report.html", category="Quality control")
    params:
        common_dir=out,
        job_name="multiqc",
    resources:
        cores=config["cores"]["default"],
        time_min=config["time_min"]["default"],
    shell:
         "multiqc -f --dirs {params.common_dir} --filename {output}"