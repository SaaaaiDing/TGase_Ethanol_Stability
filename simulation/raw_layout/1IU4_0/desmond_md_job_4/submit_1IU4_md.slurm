#!/bin/bash
#SBATCH --job-name=1IU4_MD_A100
#SBATCH --account=yjz5549_cr_default
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:a100:1
#SBATCH --time=120:00:00
#SBATCH --export=NONE
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

module purge
module load schrodinger/2023.4

WORKDIR=/storage/home/sqd5856/Desktop/default/HanLiu/prep_files/yaml/1IU4/boltz_results_1IU4/predictions/1IU4/1IU4_0/desmond_md_job_4
INPUT=desmond_md_job_4
JOBNAME=1IU4_md_batch

cd "$WORKDIR"

for ext in msj cfg cms; do
    if [[ ! -s "${INPUT}.${ext}" ]]; then
        echo "Missing or empty input: ${WORKDIR}/${INPUT}.${ext}" >&2
        exit 1
    fi
done

# Single-node localhost execution inside this batch allocation.
unset SLURM_STEP_ID SLURM_STEPID SLURM_PROCID SLURM_LOCALID SLURM_GTIDS
export SLURM_EXPORT_ENV=ALL

export SCHRODINGER_TMPDIR="/scratch/${USER}/${JOBNAME}_${SLURM_JOB_ID}"
export TMPDIR="$SCHRODINGER_TMPDIR"
mkdir -p "$SCHRODINGER_TMPDIR"

cleanup() {
    status=$?
    if [[ "$status" -eq 0 ]]; then
        rm -rf -- "$SCHRODINGER_TMPDIR"
    else
        echo "Job failed (exit ${status}); temporary files retained: ${SCHRODINGER_TMPDIR}" >&2
    fi
    exit "$status"
}
trap cleanup EXIT

ulimit -s unlimited

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "SCHRODINGER: ${SCHRODINGER}"
echo "Work directory: ${WORKDIR}"
echo "Temporary directory: ${SCHRODINGER_TMPDIR}"

nvidia-smi

"$SCHRODINGER/utilities/multisim" \
    -JOBNAME "$JOBNAME" \
    -HOST localhost \
    -maxjob 1 \
    -cpu 1 \
    -m "${INPUT}.msj" \
    -c "${INPUT}.cfg" \
    -description "Molecular Dynamics" \
    "${INPUT}.cms" \
    -mode umbrella \
    -o "${JOBNAME}-out.cms" \
    -lic DESMOND_GPGPU:16 \
    -WAIT

echo "End time: $(date)"
