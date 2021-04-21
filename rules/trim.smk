rule trim_galore_pe:
    input:
        get_fastq
    output:
        trimmed_1=temp(out + "/{sample}/trim-galore/{sample}_{pu}_val_1.fq.gz"),
        trimmed_2=temp(out + "/{sample}/trim-galore/{sample}_{pu}_val_2.fq.gz"),
        report=out + "/{sample}/trim-galore/{sample}_{pu}_1.fastq.gz_trimming_report.txt",
    params:
        outpath=out + "/{sample}/trim-galore",
        reg="{sample}.{pu}",
        job_name="trim-galore_{sample}_{pu}",
        prefix="{sample}_{pu}",
        script_dir=script_dir
    resources:
        cores=config["cores"]["trimming"],
        mem=config["mem"]["trimming"],
        time_min=config["time_min"]["trimming"],
    benchmark:
        "benchmarks/{sample}_{pu}.trim.benchmark.txt",
    shell:
         "trim_galore "
         "--cores {resources.cores} "
         "--basename {params.prefix} "
         "--fastqc "
         "--gzip "
         "--output_dir {params.outpath} "
         "--paired {input} && "
         "{params.script_dir}/rename_files.sh {params.outpath} trimming_report {params.reg}"