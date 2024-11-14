# Configuration File: `config/config.yaml`

The `config/config.yaml` file is the main configuration file for the Astra pipeline. 
It contains several sections that define paths, resources, and parameters required for pipeline execution. 
You can modify the values according to your project and computational environment.

### 1. Samples, Units, and Reheader Files

```yaml
samples: config/samples.tsv
units: config/units.tsv
reheader: config/reheader.tsv
```

- `samples`: Path to the `samples.tsv` file, which defines the list of samples to be processed.
- `units`: Path to the `units.tsv` file, which specifies the technical units associated with the samples.
- `reheader`: Optional path to a `reheader.tsv` file, used to modify sample identifiers during processing.

### 2. Paths

```yaml
paths:
    workdir: "/path/to/workdir"
    results_dir: "/path/to/results_dir"
```

- `workdir`: Directory where the pipeline will store intermediate files and logs. Ensure that there is sufficient disk space.
- `results_dir`: Directory where final results will be saved after pipeline execution.

### 3. Resources

```yaml
resources:
    reference: "/path/to/reference/reference_genome.fasta"
    regions: "/path/to/reference/regions.bed"
```

- `reference`: Path to the reference genome file in FASTA format. This will be used for alignment and variant calling.
- `regions`: Path to a BED file specifying genomic regions of interest (optional). This file helps restrict analyses to specific regions of the genome.

### 4. Parameters

```yaml
params:
  deepVariant:
    model_type: "WES"  # options [WGS, WES, PACBIO, ONT_R104, HYBRID_PACBIO_ILLUMINA]
  vep:
    resources: "/path/to/vep_resources"
    reference_version: "hg19"  # options [hg19, hg38] or [GRCh37, GRCh38]
    cache_version: "106"
```

- `deepVariant.model_type`: Defines the model type for the DeepVariant tool. Choose from the following options:
  - `WGS`: Whole genome sequencing
  - `WES`: Whole exome sequencing
  - `PACBIO`: PacBio sequencing
  - `ONT_R104`: Oxford Nanopore R104 model
  - `HYBRID_PACBIO_ILLUMINA`: Hybrid PacBio and Illumina sequencing
- `vep.resources`: Path to the resources required for the Variant Effect Predictor (VEP).
- `vep.reference_version`: Specifies the reference genome version for VEP. Options include `hg19`, `hg38`, `GRCh37`, and `GRCh38`.
- `vep.cache_version`: Defines the VEP cache version. Make sure it matches the version of the resources.

# Sample Files

The Astra pipeline requires two essential files for defining sample-related information: `config/samples.tsv` and `config/units.tsv`.

## Units File: `config/units.tsv`

In the `units.tsv` file, each row represents a **unit** for a given sample. A **unit** typically corresponds to a specific sequencing run or file associated with the sample, and it contains paths to the relevant BAM, BAI, and MD5 checksum files. This file is essential for linking the raw sequencing data with the samples defined in `samples.tsv`.

The file contains 3 tab-separated columns:

- `sample`: The generic sample name. This name will also be listed in the `samples.tsv` file. A sample can appear multiple times in the `units.tsv` file if it has multiple units associated with it (e.g., different sequencing runs or lanes).
- `unit`: A unique identifier for the unit. The `unit` is composed of three parts:
  - `flowcell_id`: An identifier for the flowcell or sequencing instrument.
  - `lane`: The lane ID for the specific sequencing run. If the BAM file comes from multiple lanes or if the lane is unknown, the lane parameter can be set to `L000`.
  - `sample_id`: The same sample ID as in the `sample` column.
- `bam`: The absolute path to the BAM file for this unit. This is the aligned sequence data.

### Example of `units.tsv`:
```tsv
sample        unit                           bam_path
SampleA       Flowcell1.L001.SampleA       /abs_path/to/data/Flowcell1.L001.SampleA.bam
SampleB       Flowcell2.L002.SampleB       /abs_path/to/data/Flowcell2.L002.SampleB.bam
SampleA       Flowcell3.L000.SampleA       /abs_path/to/data/Flowcell3.L000.SampleA.bam  # Unknown or multiple lanes
SampleC       Flowcell4.L001.SampleC       /abs_path/to/data/Flowcell4.L001.SampleC.bam
```

- **SampleA** has two units: one from lane `L001` and another from lane `L000`, where `L000` indicates that the lane is unknown or that the BAM file is from multiple lanes.

## Samples File: `config/samples.tsv`


In the `config/samples.tsv` file, each row contains information for a single **sample**, including a list of all associated units. 

The file has 2 tab-separated columns:

* **`sample`**: The generic sample name. This should match the name listed in the `units.tsv` file.
* **`units`**: A comma-separated list of units (the unit names reported in the `units.tsv` file) associated with the given sample.

### Example of `samples.tsv`:
```tsv
sample        units
SampleA       Flowcell1.L001.SampleA,Flowcell3.L000.SampleA
SampleB       Flowcell2.L002.SampleB
SampleC       Flowcell4.L001.SampleC
```

In the example:
- **SampleA** has two units associated with it, `Flowcell1.L001.SampleA` and `Flowcell3.L000.SampleA`.
- **SampleB** has one unit, `Flowcell2.L002.SampleB`.
- **SampleC** has one unit, `Flowcell4.L001.SampleC`.

This structure ensures that each sample is linked to its respective units, which can come from different sequencing runs, lanes, or flowcells.

## Hash File (Optional): `config/reheader.tsv` 

The `reheader.tsv` file contains a mapping between the LIMS (Laboratory Information Management System) identifiers and the client identifiers for each sample. This mapping is used to update or reheader the sample names in the workflow.

The file has 2 tab-separated columns:

* `LIMS`: The identifier used by the LIMS system for a sample. This identifier is often used for internal tracking and management purposes.
* `Client`: The identifier used by the client for the same sample. This name will replace the LIMS identifier at the end of the workflow, for delivering results to the client.

### Example of `config/reheader.tsv`:

```tsv
LIMS              Client
SampleA           C1234
SampleB           C5678
SampleC           C91011
```

In the example:
- The sample with LIMS identifier **`SampleA`** is associated with client identifier **`C1234`**.
- The sample with LIMS identifier **`SampleB`** is associated with client identifier **`C5678`**.
- The sample with LIMS identifier **`SampleC`** is associated with client identifier **`C91011`**.
