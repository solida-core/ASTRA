rule delivery_bam:
    input:
        bam=lambda wildcards: get_filepath_by_client_id(wildcards, folder="reads", pattern="samplename.bam"),
    output:
        bam=resolve_results_filepath('delivery',"{client}/{client}.bam"),
    conda:
        resolve_envs_filepath("samtools.yaml")
    log:
        resolve_logs_filepath('delivery_bam', "{client}.log"),
    benchmark:
        resolve_benchmarks_filepath('delivery_bam', "{client}.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=temp_path(),
    shell:
        "samtools view -H {input.bam}  | "
        "sed \"s/SM:[^\t]*/SM:{wildcards.client}/g\" | "
        "samtools reheader - {input.bam} > {output.bam} ; "
        "samtools index {output.bam} "

rule delivery_calling:
    input:
        vcf=lambda wildcards: get_filepath_by_client_id(wildcards, folder="calling", pattern="samplename/samplename.vcf.gz"),
    output:
        vcf=resolve_results_filepath('delivery',"{client}/{client}.vcf.gz"),
    params:
        header=resolve_results_filepath('delivery',"{client}/{client}.header.calling.txt"),
        sample_id=lambda wildcards: get_sample_id(wildcards.client),   
    conda:
        resolve_envs_filepath("samtools.yaml")
    log:
        resolve_logs_filepath('delivery_calling',"{client}.log"),
    benchmark:
        resolve_benchmarks_filepath('delivery_calling',"{client}.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=temp_path(),
    shell:
        "bcftools view -h {input.vcf} | sed \"s/{params.sample_id}/{wildcards.client}/g\" > {params.header} ; "
        "bcftools reheader -h {params.header} -o {output.vcf} {input.vcf} ; "
        "rm -f {params.header} ; "
        "tabix -p vcf {output.vcf} "

rule delivery_annotation:
    input:
        vcf=lambda wildcards: get_filepath_by_client_id(wildcards, folder="annotation", pattern="samplename/samplename.annotated.vcf"),
        tsv=lambda wildcards: get_filepath_by_client_id(wildcards,folder="annotation",pattern="samplename/samplename.annotated.tsv"),
    output:
        vcf=temp(resolve_results_filepath('delivery',"{client}/{client}.annotated.vcf")),
        tsv=resolve_results_filepath('delivery',"{client}/{client}.annotated.tsv"),
        gz=resolve_results_filepath('delivery',"{client}/{client}.annotated.vcf.gz"),
    params:
        header=resolve_results_filepath('delivery',"{client}/{client}.header.annotated.txt"),
        sample_id=lambda wildcards: get_sample_id(wildcards.client),   
    conda:
        resolve_envs_filepath("samtools.yaml")
    log:
        resolve_logs_filepath('delivery_annotated',"{client}.log"),
    benchmark:
        resolve_benchmarks_filepath('delivery_annotated',"{client}.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=temp_path(),
    shell:
        "bcftools view -h {input.vcf} | sed \"s/{params.sample_id}/{wildcards.client}/g\" > {params.header} ; "
        "bcftools reheader -h {params.header} -o {output.vcf} {input.vcf} ; "
        "rm -f {params.header} ; "
        "bgzip -c {output.vcf} > {output.gz} ; "
        "tabix -p vcf {output.gz} ; "
        "sed 's/{params.sample_id}/{wildcards.client}/g' {input.tsv} > {output.tsv} "

rule delivery_metrics:
    input:
        dat=lambda wildcards: get_filepath_by_client_id(wildcards,folder="metrics",pattern="samplename.hsmetrics.dat"),
    output:
        dat=resolve_results_filepath('delivery',"{client}/{client}.hsmetrics.dat"),
    conda:
        resolve_envs_filepath("rsync.yaml")
    log:
        resolve_logs_filepath('delivery_metrics',"{client}.log"),
    benchmark:
        resolve_benchmarks_filepath('delivery_metrics',"{client}.txt"),
    threads: conservative_cpu_count(reserve_cores=2,max_cores=99)
    resources:
        tmpdir=temp_path(),
    shell:
        "rsync -az {input.dat} {output.dat} "