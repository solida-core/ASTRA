rule delivery_bam:
    input:
        bam=lambda wildcards: get_sample_by_client(wildcards, folder="reads", filename="{sample}.bam"),
    output:
        bam=resolve_results_filepath('delivery',"{Client}/{Client}.bam"),
    params:
        client_id=get_client_id_by_sample(sample)
    conda:
        resolve_envs_filepath("samtools.yaml")
    log:
        resolve_logs_filepath('reheader', "{Client}.log"),
    benchmark:
        resolve_benchmarks_filepath('reheader', "{Client}.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=temp_path(),
    shell:
        "samtools view -H {input.bam}  | "
        "sed \"s/SM:[^\t]*/SM:TEST_SAMPLE_NAME/g\" | "
        "samtools reheader - {input.bam} > {output.bam} "
        "&& "
        "samtools index "
        "{output.bam} "
