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
        reads=get_trimmed_reads,
        ref_dir=out + "/star_index"
    output:
        genomic=temp(out + "/{sample}/star_aln/{lane}_{sample}_{sample_number}_trim_star.Aligned.sortedByCoord.out.bam"),
        transcipts=temp(out + "/{sample}/star_aln/{lane}_{sample}_{sample_number}_trim_star.Aligned.toTranscriptome.out.bam"),
        log=out + "/{sample}/star_aln/{lane}_{sample}_{sample_number}_trim_star.Log.final.out"
    params:
        prefix=out + "/{sample}/star_aln/{lane}_{sample}_{sample_number}_trim_star.",
        job_name="star_{lane}_{sample}_{sample_number}",
        rg=get_read_group,
    benchmark:
        "benchmarks/{lane}_{sample}_{sample_number}.starNindex.benchmark.txt",
    resources:
        cores=config["cores"]["mapping"],
        mem=config["mem"]["mapping"],
        time_min=config["time_min"]["mapping"],
    shell:
        "STAR --genomeDir {input.ref_dir} "
        "--readFilesIn {input.reads} "
        "--readFilesCommand zcat "        
        "--twopassMode Basic "
        "--quantMode TranscriptomeSAM " # for RSEM
        "--peOverlapNbasesMin 10 "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMattrRGline {params.rg} "
        "--runThreadN {resources.cores} "

rule samtools_merge:
    input:
        genomic_bams=get_sample_bams_genomic,
        transcripts_bams=get_sample_bams_transcript,
    output:
        genomic=out + "/{sample}/star_aln/{sample}_{sample_number}_trim_star.Aligned.sortedByCoord.out.bam",
        transcripts=out + "/{sample}/star_aln/{sample}_{sample_number}_trim_star.Aligned.toTranscriptome.out.bam",
    params:
        job_name="mergeNindex_{sample}_{sample_number}",
    benchmark:
        "benchmarks/{sample}_{sample_number}.mergeNindex.benchmark.txt",
    resources:
        cores=config["cores"]["genome_index"],
        mem=config["mem"]["genome_index"],
        time_min=config["time_min"]["genome_index"],
    shell:
        """
        samtools merge -@ {resources.cores} -cp {output.genomic} {input.genomic_bams} &&
        samtools merge -@ {resources.cores} -cp {output.transcripts} {input.transcripts_bams} &&
        samtools index -b -@ {resources.cores} {output.genomic} &&
        samtools index -b -@ {resources.cores} {output.transcripts};
        """

rule mark_dupes:
    input:
        out + "/{sample}/star_aln/{sample}_{sample_number}_trim_star.Aligned.sortedByCoord.out.bam"
    output:
        bam=out + "/{sample}/picard_markdupe/{sample}_{sample_number}_trim_star_marked.bam",
        metrics=out + "/{sample}/picard_markdupe/{sample}_{sample_number}_trim_star_marked.metrics.txt"
    params:
        tmp_dir=tmp,
        job_name="mdupeNindex_{sample}_{sample_number}",
    benchmark:
        "benchmarks/{sample}_{sample_number}.mark_dupes.benchmark.txt",
    resources:
        cores=config["cores"]["mapping"],
        mem=config["mem"]["mapping"],
        time_min=config["time_min"]["mapping"],
    shell:
        "tmpdir=$(mktemp --directory {params.tmp_dir}/tmp.XXXXX) &&"
        "java -jar /home/proj/bin/conda/envs/D_rna-seq_MW/share/picard-2.22.1-0/picard.jar MarkDuplicates "
        "I={input} "
        "O={output.bam} "
        "M={output.metrics} "
        "REMOVE_DUPLICATES=true "
        "TMP_DIR=${{tmpdir}} && "
        "samtools index -b -@ {resources.cores} {output.bam}; "



