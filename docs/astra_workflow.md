# ASTRA WORKFLOW
**ASTRA** (Advanced Snakemake Pipeline for Thorough Variant Analysis) 
is a pipeline designed for comprehensive variant calling and annotation from Next-Generation Sequencing (NGS) data.

______________________________



**ASTRA** processes **BAM files** to perform variant calling and annotation. The workflow is designed for efficient and accurate analysis of Next-Generation Sequencing (NGS) data, and it includes the following key steps:

1. **BAM Merging**  
   - If a sample has multiple BAM files (e.g., from different sequencing runs or lanes), these are merged into a single BAM file per sample. This step ensures that downstream analyses are performed on consolidated data.

2. **Hybrid-Selection (HS) Metrics Collection**  
   - This step computes metrics specific to hybrid-selection datasets, which are commonly used in targeted sequencing experiments such as exome sequencing.  
   - The tool used is **GATK CollectHsMetrics**, which provides detailed statistics on target coverage, hybrid-selection efficiency, and potential biases in the sequencing data.

3. **Variant Calling**  
   - The merged BAM files are processed for variant calling using **DeepVariant**, a deep learning-based variant caller.  
   - **DeepVariant** is highly accurate and works across multiple sequencing platforms, detecting SNPs and indels with high sensitivity and specificity.  
   - The workflow supports various DeepVariant model types, such as WGS, WES, PACBIO, ONT, and hybrid sequencing platforms, as specified in the configuration file.

4. **Variant Annotation**  
   - Detected variants are annotated using **VEP (Variant Effect Predictor)**.  
   - **VEP** integrates information from multiple genomic databases to predict the functional effects of variants, such as coding impact, regulatory region changes, and associated diseases.  
   - The annotations include both functional insights and population frequency data, aiding downstream interpretation.

5. **Results Generation**  
   - For each sample, the workflow produces:  
     - A **raw VCF file** containing the detected variants.  
     - An **annotated VCF file** with detailed variant annotations.  
     - A **tabular annotation file** for easier integration with other tools or reports.  
     - A **metrics file** summarizing sequencing and variant quality metrics.

This pipeline ensures comprehensive and streamlined processing of NGS data, making it suitable for both research and clinical applications.

A complete view of the analysis workflow is provided by the pipeline's [graph](images/astra.png).



