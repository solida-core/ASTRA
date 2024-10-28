#######################
import errno
import pandas as pd
import os
import multiprocessing
import psutil
from snakemake.utils import validate

report: "../report/workflow.rst"

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config.get("samples"), index_col="sample")
units = pd.read_table(config.get("units"), index_col=["unit"], dtype=str)
reheader = pd.read_table(config.get("reheader"), index_col=["LIMS"], dtype=str)

def resolve_single_filepath(basepath, filename):
    return os.path.join(basepath, filename)

def resolve_results_filepath(dirname, filename):
    path = os.path.join(config.get('paths').get('results_dir'), dirname)
    return resolve_single_filepath(path, filename)

def resolve_logs_filepath(dirname, filename):
    path = os.path.join(config.get('paths').get('results_dir'), 'logs', dirname)
    return resolve_single_filepath(path, filename)

def resolve_benchmarks_filepath(dirname, filename):
    path = os.path.join(config.get('paths').get('results_dir'), 'benchmarks', dirname)
    return resolve_single_filepath(path, filename)

def expand_filepath(filepath):
    filepath = os.path.expandvars(os.path.expanduser(filepath))
    if not os.path.isabs(filepath):
        raise FileNotFoundError(
            errno.ENOENT,
            os.strerror(errno.ENOENT) + " (path must be absolute)",
            filepath,
        )
    return filepath

def  get_units_by_sample(sample, label="units",):
    return [ unit for unit in samples.loc[sample, [label]][0].split(",")]

def get_bams_by_sample(wildcards, samples, label='bam'):
    sample = wildcards.sample
    _units = get_units_by_sample(sample)
    return units.loc[_units, label].tolist()

def cpu_count():
    return multiprocessing.cpu_count()

def conservative_cpu_count(reserve_cores=1, max_cores=8):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return max(cores - reserve_cores, 1)

