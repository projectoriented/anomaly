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
        star_aln=out + "/{sample}/star_aln/{lane}_{sample}_{sample_number}_trim_star.Aligned.toTranscriptome.out.bam",
    output:
        report(out + "/{sample}/rsem/{lane}_{sample}_{sample_number}.stat/{lane}_{sample}_{sample_number}.cnt")
    params:
        prefix=out + "/{sample}/rsem/{lane}_{sample}_{sample_number}",
        job_name="rsem_{lane}_{sample}_{sample_number}",
        ref_prefix=out + "/rsem_index/hg38"
    benchmark:
        "benchmarks/{lane}_{sample}_{sample_number}.rsem.benchmark.txt",
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