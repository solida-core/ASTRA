## User files
**ASTRA** (Advanced Snakemake Pipeline for Thorough Variant Analysis)
processes deduplicated and recalibrated BAM files to generate high-quality outputs, including raw VCF files and fully annotated VCF files.
The standardization provided by [solida-core](https://github.com/solida-core) requires perhaps accessory user-defined files to have a given organization.

### Required Input Files
**ASTRA** pipeline requires three different mandatory user inputs:
* [samples](../config/README.md#samples-file-configsamplestsv)
* [units](../config/README.md#units-file-configunitstsv)
* [reheader](../config/README.md#hash-file-optional-configreheadertsv-)

These files contain **tab-separated** information about samples to be analyzed and must be declared into the `config.yaml` file:
```
samples: "path_to_input_files/samples.tsv"
units: "path_to_input_files/units.tsv"
reheader: "path_to_input_files/reheader.tsv"
```

For more details on configuration files and required user files, refer to the [detailed documentation](../config/README.md).