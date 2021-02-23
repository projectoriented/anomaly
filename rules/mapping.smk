rule star_index:
    input:
         fa = config["genome"],
         gtf = config["gtf"]
    output:
         directory(out + "/star_index")
    resources:
        cores=config["cores"]["genome_index"],
        mem=config["mem"]["genome_index"],
        time_min=config["time_min"]["genome_index"],
    params:
        job_name="star_index",
    benchmark:
        "benchmarks/star_index.benchmark.txt",
    shell:
         "mkdir {output} && "
         "STAR --runThreadN {resources.cores} " 
         "--runMode genomeGenerate "
         "--genomeDir {output} "
         "--genomeFastaFiles {input.fa} "
         "--sjdbGTFfile {input.gtf} "
         "--sjdbOverhang 99"

rule star_align:
    input:
        read1=get_trimmed_read1,
        read2=get_trimmed_read2,
        ref_dir=out + "/star_index"
    output:
        genomic=out + "/{sample}/star_aln/{sample}_trim_star.Aligned.sortedByCoord.out.bam",
        transcipts=out + "/{sample}/star_aln/{sample}_trim_star.Aligned.toTranscriptome.out.bam",
        log=out + "/{sample}/star_aln/{sample}_trim_star.Log.final.out"
    params:
        prefix=out + "/{sample}/star_aln/{sample}_trim_star.",
        job_name="starNindex_{sample}",
        rg=get_read_group,
    benchmark:
        "benchmarks/{sample}.starNindex.benchmark.txt",
    resources:
        cores=config["cores"]["mapping"],
        mem=config["mem"]["mapping"],
        time_min=config["time_min"]["mapping"],
    shell:
        "format1=$(echo {input.read1} | tr ' ' ',');"
        "format2=$(echo {input.read2} | tr ' ' ',');"
        "STAR --genomeDir {input.ref_dir} "
        "--readFilesIn ${{format1}} ${{format2}} "
        "--readFilesCommand zcat "        
        "--twopassMode Basic "
        "--quantMode TranscriptomeSAM " # for RSEM
        "--peOverlapNbasesMin 10 "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMattrRGline {params.rg} "
        "--runThreadN {resources.cores} && "
        "samtools index -b -@ {resources.cores} {output.genomic}; "

rule mark_dupes:
    input:
        out + "/{sample}/star_aln/{sample}_trim_star.Aligned.sortedByCoord.out.bam"
    output:
        bam=out + "/{sample}/picard_markdupe/{sample}_trim_star_marked.bam",
        metrics=out + "/{sample}/picard_markdupe/{sample}_trim_star_marked.metrics.txt"
    params:
        tmp_dir=tmp,
        job_name="mdupeNindex_{sample}",
    benchmark:
        "benchmarks/{sample}.mark_dupes.benchmark.txt",
    resources:
        cores=config["cores"]["mapping"],
        mem=config["mem"]["mapping"],
        time_min=config["time_min"]["mapping"],
    shell:
        "tmpdir=$(mktemp --directory {params.tmp_dir}/tmp.XXXXX) &&"        
        "gatk MarkDuplicates "
        "-I {input} "
        "-O {output.bam} "
        "-M {output.metrics} "
        "-REMOVE_DUPLICATES true "
        "-TMP_DIR ${{tmpdir}} && "
        "samtools index -b -@ {resources.cores} {output.bam}; "



