## merge_fastq: Demultiplexed samples will come from more than 1 lane and these samples need to be merged/concat.
##              If sample isn't multiplexed, then this is skipped and copied to target directory.
rule merge_fastq:
    input:
        fq1=get_read1,
        fq2=get_read2,
    output:
        merged_1=temp(out + "/{sample}/fastq/{sample}_1.fastq.gz"),
        merged_2=temp(out + "/{sample}/fastq/{sample}_2.fastq.gz")
    resources:
        cores=config["cores"]["default"],
        mem=config["mem"]["default"],
        time_min=config["time_min"]["default"],
    run:
        if len(input.fq1) == 1:
            shell("cp {input.fq1} {output.merged_1}; \
                    cp {input.fq2} {output.merged_2};")
        else:
            shell("cat {input.fq1} > {output.merged_1}; \
                    cat {input.fq2} > {output.merged_2};")
