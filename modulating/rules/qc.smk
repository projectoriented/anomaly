out = config["common_out"]
tmp = config["tmp"]

rule fastqc:
    input:
        get_fastq
    output:
        html=out + "/{sample}/fastqc/{lane}_{sample}_{sample_number}.html",
        zip=out + "/{sample}/fastqc/{lane}_{sample}_{sample_number}.zip",
    params:
        tmp_dir=tmp,
        outpath=out + "/{sample}/fastqc/",
        job_name="fastqc_{lane}_{sample}_{sample_number}"
    shell:
        "tmpdir=$(mktemp --directory {params.tmp_dir}/tmp.XXXXX) && "
        "fastqc --noextract --dir ${{tmpdir}} --outdir {params.outpath} {input};"

rule multiqc:
    input:
        expand([out + "/{s.sample}/fastqc/{s.lane}_{s.sample}_{s.sample_number}.zip"], s=samples.itertuples())
    output:
        out + "/multiqc_report.html"
    params:
        common_dir=out,
        job_name="multiqc",
    shell:
         "multiqc {params.common_dir} --filename {output}"