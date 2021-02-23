rule rsem_index:
    input:
         fa = config["genome"],
         gtf = config["gtf"]
    output:
         directory(out + "/rsem_index")
    resources:
        cores=config["cores"]["genome_index"],
        mem=config["mem"]["genome_index"],
        time_min=config["time_min"]["genome_index"],
    params:
        job_name="rsem_index_hg38",
        prefix=out + "/rsem_index/hg38"
    shell:
        "mkdir {output} && "
        "rsem-prepare-reference --num-threads {resources.cores} " 
        "--gtf {input.gtf} "
        "{input.fa} "
        "{params.prefix}"

rule rsem:
    input:
        index=out + "/rsem_index",
        star_aln=out + "/{sample}/star_aln/{sample}_trim_star.Aligned.toTranscriptome.out.bam",
    output:
        report(out + "/{sample}/rsem/{sample}.stat/{sample}.cnt")
    params:
        prefix=out + "/{sample}/rsem/{sample}",
        job_name="rsem_{sample}",
        ref_prefix=out + "/rsem_index/hg38"
    benchmark:
        "benchmarks/{sample}.rsem.benchmark.txt",
    resources:
        cores=config["cores"]["counting"],
        mem=config["mem"]["counting"],
        time_min=config["time_min"]["counting"],
    shell:
        "rsem-calculate-expression "
        "--paired-end "
        "--bam "
        "--no-bam-output "
        "--num-threads {resources.cores} "        
        "{input.star_aln} "
        "{params.ref_prefix} "
        "{params.prefix}"