rule star_index:
    input:
         fa = config["genome"],
         gtf = config["gtf"]
    output:
         directory(out + "/star_index")
    resources:
        cores=18,
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
        genomic=out + "/{sample}/star_aln/{lane}_{sample}_{sample_number}_trim_star.Aligned.sortedByCoord.out.bam",
        transcipts=out + "/{sample}/star_aln/{lane}_{sample}_{sample_number}_trim_star.Aligned.toTranscriptome.out.bam",
        log=out + "/{sample}/star_aln/{lane}_{sample}_{sample_number}_trim_star.Log.final.out"
    params:
        prefix=out + "/{sample}/star_aln/{lane}_{sample}_{sample_number}_trim_star.",
        job_name="starNindex_{lane}_{sample}_{sample_number}",
        rg=get_read_group,
    resources:
        cores=18,
        mem_mb=90000,
        time_min=480,
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
        "--runThreadN {resources.cores} && "
        "samtools index -b -@ {resources.cores} {output.genomic}; "
        "samtools index -b -@ {resources.cores} {output.transcipts};"

rule mark_dupes:
    input:
        out + "/{sample}/star_aln/{lane}_{sample}_{sample_number}_trim_star.Aligned.sortedByCoord.out.bam"
    output:
        bam=out + "/{sample}/picard_markdupe/{lane}_{sample}_{sample_number}_trim_star_marked.bam",
        metrics=out + "/{sample}/picard_markdupe/{lane}_{sample}_{sample_number}_trim_star_marked.metrics.txt"
    params:
        tmp_dir=tmp,
        job_name="mdupeNindex_{lane}_{sample}_{sample_number}",
    resources:
        cores=18,
        mem_mb=90000,
        time_min=180,
    shell:
        "tmpdir=$(mktemp --directory {params.tmp_dir}/tmp.XXXXX) &&"
        "java -jar /home/proj/bin/conda/envs/D_rna-seq_MW/share/picard-2.22.1-0/picard.jar MarkDuplicates "
        "I={input} "
        "O={output.bam} "
        "M={output.metrics} "
        "REMOVE_DUPLICATES=true "
        "TMP_DIR=${{tmpdir}} && "
        "samtools index -b -@ {resources.cores} {output.bam}; "