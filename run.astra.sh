#!/usr/bin bash
#
# ASTRA launcher
#

usage="$(basename "$0") [-h] [-p \"parameters\"]

where:
    -h  show this help text
    -p  snakemake parameters as \"--rerun-incomplete --dryrun --keep-going --profile label_of_cluster_profile\"
"

while getopts ':h:p:' option; do
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
CFGFILE="${WORKDIR}/config/config.yaml"

REFFILE=$(yq -Y '.resources.reference' "${CFGFILE}")
BEDFILE=$(yq -Y '.resources.regions' "${CFGFILE}")
OUTDIR=$(yq -Y '.paths.results_dir' "${CFGFILE}")

REFDIR=$(dirname "$REFFILE")
BEDDIR=$(dirname "$BEDFILE")
OUTDIR=$(dirname "$OUTDIR")

snakemake \
  --use-conda \
  --use-apptainer \
  --apptainer-args " --bind ${OUTDIR} --bind ${REFDIR} --bind ${BEDDIR} " \
  --printshellcmds \
  --restart-times 1 \
  --latency-wait 120 \
  --jobname "astra.{rulename}.{jobid}" \
  ${SM_PARAMETERS}