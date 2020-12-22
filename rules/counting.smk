rule rsem_index:
    input:
         fa = config["genome"],
         gtf = config["gtf"]
    output:
         directory(out + "/rsem_index")
    resources:
        cores=9,
    params:
        job_name = "rsem_index_hg38",
        prefix = out + "/rsem_index/hg38"
    shell:
        'rsem-prepare-reference --num-threads {resources.cores} '
        '--gtf {input.gtf} '
        '{input.fa} '
        '{params.prefix}'

rule rsem:
    input:
        star_aln=out + "/{sample}/star_aln/{lane}_{sample}_{sample_number}_trim_star.Aligned.toTranscriptome.out.bam",
        ref_dir=out + "/rsem_index"
    output:
        out + "/{sample}/rsem/{lane}_{sample}_{sample_number}.stat/{lane}_{sample}_{sample_number}.cnt",
    params:
        prefix=out + "/{sample}/rsem/{lane}_{sample}_{sample_number}",
        job_name="rsem_{lane}_{sample}_{sample_number}",
        ref_prefix=out + "/rsem_index/hg38"
    resources:
        cores=18,
        mem_mb=90000,
        time_min=480,
    shell:
        "rsem-calculate-expression "
        "--paired-end "
        "--bam "
        "--num-threads {resources.cores} "        
        "{input.star_aln} "
        "{params.ref_prefix} "
        "{params.prefix}"