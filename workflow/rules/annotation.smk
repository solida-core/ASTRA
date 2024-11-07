rule vep2vcf:
    input:
        vcf=rules.call_variants.output,
    output:
        vcf=resolve_results_filepath('annotation',"{sample}/{sample}.annotated.vcf"),
    params:
        genome_version=get_vep_genome_version(),
        resources=config.get("params").get('vep').get("resources"),
        cache_version=config.get("params").get('vep').get("cache_version"),
    conda:
        resolve_envs_filepath("vep.yaml")
    log:
        resolve_logs_filepath('annotation',"{sample}.vcf.log"),
    benchmark:
        resolve_benchmarks_filepath('annotation',"{sample}.vcf.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=temp_path(),
    shell:
        "vep "
        "--input_file {input.vcf} "
        "--output_file {output.vcf} "
        "--dir {params.resources} "
        "--assembly {params.genome_version} "
        "--cache --cache_version {params.cache_version} "
        "--offline --everything "
        "--fork {threads} "
        "--vcf "

rule vep2tsv:
    input:
        vcf=rules.call_variants.output,
    output:
        vcf=resolve_results_filepath('annotation',"{sample}/{sample}.annotated.tsv"),
    params:
        genome_version=get_vep_genome_version(),
        resources=config.get("params").get('vep').get("resources"),
        cache_version=config.get("params").get('vep').get("cache_version"),
    conda:
        resolve_envs_filepath("vep.yaml")
    log:
        resolve_logs_filepath('annotation',"{sample}.tsv.log"),
    benchmark:
        resolve_benchmarks_filepath('annotation',"{sample}.tsv.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=temp_path(),
    shell:
        "vep "
        "--input_file {input.vcf} "
        "--output_file {output.vcf} "
        "--dir {params.resources} "
        "--assembly {params.genome_version} "
        "--cache --cache_version {params.cache_version} "
        "--offline --everything "
        "--fork {threads} "
        "--tsv "
