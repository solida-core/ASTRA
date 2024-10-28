rule merge:
    input:
        lambda wildcards: get_bams_by_sample(wildcards),
    output: 
        bam=resolve_results_filepath('reads', "{sample}.cram"),
    params:
        cmd="samtools",
        genome=config.get("resources").get("reference"),
        output_fmt="CRAM",
    conda:
        "../envs/samtools.yaml"
    log:
        resolve_logs_filepath('merge', "{sample}.log"),
    benchmark:
        resolve_benchmarks_filepath('merge', "{sample}.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    script:
        "workflow/scripts/merge.py"

rule index:
    input:
        bam=rules.merge.output
    output:
        bai=resolve_results_filepath('reads', "{sample}.cram.crai")
    conda:
        "../envs/samtools.yaml"
    log:
        resolve_logs_filepath('index',"{sample}.log"),
    benchmark:
        resolve_benchmarks_filepath('index', "{sample}.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "samtoools index "
        "--threads {threads} "
        "--output {output.bai} "
        "{input.bam} "
        ">& {log} "




