
rule merge:
    input:
        lambda wildcards: get_bams_by_sample(wildcards),
    output: 
        temp(
            resolve_results_filepath('reads', "{sample}.cram")
        ),
    params:
        cmd="samtools",
        genome=config.get("resources").get("reference"),
        output_fmt="CRAM",
    conda:
        resolve_envs_filepath("samtools.yaml")
    log:
        resolve_logs_filepath('merge', "{sample}.log"),
    benchmark:
        resolve_benchmarks_filepath('merge', "{sample}.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=temp_path(),
    script:
        "../scripts/merge.py"

rule index:
    input:
        bam=rules.merge.output
    output:
        bai=resolve_results_filepath('reads', "{sample}.cram.crai")
    conda:
        resolve_envs_filepath("samtools.yaml")
    log:
        resolve_logs_filepath('index',"{sample}.log"),
    benchmark:
        resolve_benchmarks_filepath('index', "{sample}.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=temp_path(),

    shell:
        "samtoools index "
        "--threads {threads} "
        "--output {output.bai} "
        "{input.bam} "
        ">& {log} "




