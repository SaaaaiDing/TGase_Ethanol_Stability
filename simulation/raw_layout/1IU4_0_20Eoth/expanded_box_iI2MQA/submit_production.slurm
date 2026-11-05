#!/bin/bash
#SBATCH --job-name=1IU4_EtOH20_200ns
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

SOURCE_DIR="$SLURM_SUBMIT_DIR"
JOBNAME="1IU4_EtOH20_md_${SLURM_JOB_ID}"
RUN_DIR="${SOURCE_DIR}/production_${SLURM_JOB_ID}"

cd "$SOURCE_DIR"

for file in equilibrated.cms production.msj production.cfg settings.json; do
    [[ -s "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
done

grep -q '^Multisim completed\.' 1IU4_EtOH20_large_equil_multisim.log || {
    echo "The equilibration log does not confirm successful completion." >&2
    exit 1
}

"$SCHRODINGER/run" python3 - <<'PY'
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd().parent))
import expand_1IU4_EtOH20 as workflow

st = workflow.read_full("equilibrated.cms")
groups = workflow.audit(st, workflow.NEW_COUNTS)
workflow.forcefield_audit("equilibrated.cms")

gap, diameter = workflow.solute_geometry(st, groups)
box = workflow.box(st)
density = st.total_weight * 1.66053906660 / box.prod()
settings = json.loads(Path("settings.json").read_text())
minimum_gap = max(float(settings["cutoff_A"]), 15.0)

print(f"Final box (A): {box.tolist()}")
print(f"Final density (g/cm3): {density:.4f}")
print(f"Solute image distance (A): {gap:.2f}")
print(f"Production duration (ns): {settings['production_ns']}")

if not 0.85 < density < 1.30:
    raise SystemExit("Final density needs review before production MD.")
if gap <= minimum_gap:
    raise SystemExit("The solute is too close to its periodic image.")
PY

mkdir -p "$RUN_DIR"
cp equilibrated.cms production.msj production.cfg settings.json "$RUN_DIR/"
cd "$RUN_DIR"

unset SLURM_STEP_ID SLURM_STEPID SLURM_PROCID SLURM_LOCALID SLURM_GTIDS
unset OMPI_MCA_plm_slurm_args PMI_RANK PMI_SIZE PMI_FD PMIX_RANK PMIX_NAMESPACE
export SLURM_EXPORT_ENV=ALL

export SCHRODINGER_TMPDIR="/scratch/${USER}/${JOBNAME}"
export TMPDIR="$SCHRODINGER_TMPDIR"
mkdir -p "$SCHRODINGER_TMPDIR"

ulimit -s unlimited

echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Run directory: $RUN_DIR"
echo "Start time: $(date)"
nvidia-smi

"$SCHRODINGER/utilities/multisim" \
    -JOBNAME "$JOBNAME" \
    -HOST localhost \
    -maxjob 1 \
    -cpu 1 \
    -m production.msj \
    -c production.cfg \
    -description "1IU4 20 percent ethanol production MD" \
    equilibrated.cms \
    -mode umbrella \
    -o "${JOBNAME}-out.cms" \
    -lic DESMOND_GPGPU:16 \
    -WAIT

[[ -s "${JOBNAME}-out.cms" ]] || {
    echo "Production output CMS is missing." >&2
    exit 1
}

echo "End time: $(date)"
echo "Output directory: $RUN_DIR"
