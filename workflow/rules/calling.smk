rule call_variants:
    input:
        bam=rules.preprocessing.output.bam,
        bai=rules.preprocessing.output.bai
    output:
        vcf=resolve_results_filepath('calling', "{sample}/{sample}.vcf.gz"),
    params:
        genome=config.get("resources").get("reference"),
        regions=config.get("resources").get("regions"),
        type=config.get("params").get("deepVariant").get("model_type"),
        version=config.get("params").get("deepVariant").get("version"),
        outdir=resolve_results_filepath('calling', "{sample}"),
    singularity:
        "docker://google/deepvariant:1.6.1"
    log:
        resolve_logs_filepath('calling', "{sample}.log"),
    benchmark:
        resolve_benchmarks_filepath('calling', "{sample}.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=temp_path(),
    shell:
        "/opt/deepvariant/bin/run_deepvariant "
        "--model_type={params.type} "
        "--ref={params.genome} "
        "--reads={input.cram} "
        "--regions {params.regions} "
        "--output_vcf={output} "
        "--intermediate_results_dir {params.outdir} "
        "--num_shards=1 "
