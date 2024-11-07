#!/usr/bin bash
#
# ASTRA launcher
#

usage="$(basename "$0") [-h] [-p \"parameters\"]

where:
    -h  show this help text
    -p  snakemake parameters as \"--rerun-incomplete --dryrun --keep-going --profile label_of_cluster_profile\"
"

while getopts ':hc:p:w:s:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    p) SM_PARAMETERS=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    *) echo "$usage"
       exit
  esac
done

shift $((OPTIND - 1))


WORKDIR=$(pwd)
CFGFILE="${WORKDIR}/config/config.yml"

REFFILE=$(yq e '.resources.reference' "${CFGFILE}")
BEDFILE=$(yq e '.resources.regions' "${CFGFILE}")
OUTDIR=$(yq e '.paths.results_dir' "${CFGFILE}")

REFDIR=$(dirname "$REFFILE")
BEDDIR=$(dirname "$BEDFILE")

snakemake \
  --use-conda \
  --use-apptainer \
  --apptainer-args "--bind ${OUTDIR} --bind ${REFDIR} --bind ${BEDDIR} " \
  --printshellcmds \
  --restart-times 3 \
  --latency-wait 120 \
  --jobname "astra.{rulename}.{jobid}" \
  ${SM_PARAMETERS}