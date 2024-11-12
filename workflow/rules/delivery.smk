rule delivery_bam:
    input:
        bam=lambda wildcards: get_sample_by_client(wildcards, folder="reads", pattern="samplename.bam"),
    output:
        bam=resolve_results_filepath('delivery',"{client}/{client}.bam"),
    conda:
        resolve_envs_filepath("samtools.yaml")
    log:
        resolve_logs_filepath('reheader', "{client}.log"),
    benchmark:
        resolve_benchmarks_filepath('reheader', "{client}.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=temp_path(),
    shell:
        "samtools view -H {input.bam}  | "
        "sed \"s/SM:[^\t]*/SM:{wildcards.client}/g\" | "
        "samtools reheader - {input.bam} > {output.bam} "
        "&& "
        "samtools index "
        "{output.bam} "
