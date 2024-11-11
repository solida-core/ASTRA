#######################
import errno
import pandas as pd
import os
import multiprocessing
import psutil
from snakemake.utils import validate

report: "../report/workflow.rst"

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config.get("samples"), sep='\t')
units = pd.read_csv(config.get("units"), sep='\t')
reheader = pd.read_csv(config.get("reheader"), sep='\t')

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

def resolve_envs_filepath(filename):
    path = os.path.join(config.get('paths').get('workdir'), 'workflow', 'envs')
    return resolve_single_filepath(path, filename)

def resolve_scripts_filepath(filename):
    path = os.path.join(config.get('paths').get('workdir'), 'workflow', 'scripts')
    return resolve_single_filepath(path, filename)

def temp_path(path=""):
    default_path = os.path.join(config.get('paths').get('results_dir'), 'tmp')
    if path:
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                return default_path
        return path
    return default_path

def expand_filepath(filepath):
    filepath = os.path.expandvars(os.path.expanduser(filepath))
    if not os.path.isabs(filepath):
        raise FileNotFoundError(
            errno.ENOENT,
            os.strerror(errno.ENOENT) + " (path must be absolute)",
            filepath,
        )
    return filepath

def get_bams_by_sample(wildcards):
    return units.loc[units['sample'] == wildcards.sample, 'bam'].tolist()

def get_vep_genome_version(version=None):
    version = version if version else config.get("params").get("vep").get("reference_version")
    if version in ['hg19', 'hg38']:
        return 'GRCh37' if version in 'hg19' else 'GRCh38'
    return version

def cpu_count():
    return multiprocessing.cpu_count()

def conservative_cpu_count(reserve_cores=1, max_cores=8):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return max(cores - reserve_cores, 1)

def get_client_id_by_sample(sample_id):
    return reheader.loc[reheader['LIMS'] == sample_id, 'Client']

def resolve_bam_delivery_filepath(sample_id):
    client_id = get_client_id_by_sample(sample_id)
    return resolve_results_filepath("delivery",f"{client_id}/{client_id}.bam")