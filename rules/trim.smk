rule trim_galore_pe:
    input:
        get_fastq
    output:
        trimmed_1=temp(out + "/{sample}/trim-galore/{lane}_{sample}_{sample_number}_val_1.fq.gz"),
        trimmed_2=temp(out + "/{sample}/trim-galore/{lane}_{sample}_{sample_number}_val_2.fq.gz"),
        report=out + "/{sample}/trim-galore/{lane}_{sample}_{sample_number}_1.fastq.gz_trimming_report.txt",
    params:
        outpath=out + "/{sample}/trim-galore",
        job_name="trim-galore_{lane}_{sample}_{sample_number}",
        prefix="{lane}_{sample}_{sample_number}",
        capture_group="{sample}_{sample_number}",
        script_dir=script_dir
    resources:
        cores=9,
    shell:
         "trim_galore "
         "--cores {resources.cores} "
         "--basename {params.prefix} "
         "--fastqc "
         "--gzip "
         "--output_dir {params.outpath} "
         "--paired {input} && "
         "{params.script_dir}/rename_files.sh {params.outpath} trimming_report {wildcards.lane} {params.capture_group}"

