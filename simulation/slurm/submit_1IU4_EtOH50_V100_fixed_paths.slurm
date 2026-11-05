#!/bin/bash
#SBATCH --job-name=1IU4_EtOH50_200ns
#SBATCH --account=yjz5549_cr_default
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --gres=gpu:v100:1
#SBATCH --time=120:00:00
#SBATCH --export=NONE
#SBATCH --chdir=/storage/group/yjz5549/default/HanLiu/prep_files/yaml/1IU4/boltz_results_1IU4/predictions/1IU4/1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

module purge
module load schrodinger/2023.4

SYSTEM_ROOT=/storage/group/yjz5549/default/HanLiu/prep_files/yaml/1IU4/boltz_results_1IU4/predictions/1IU4/1IU4_0_20Eoth/1IU4_EtOH50
CFG_DIR="${SYSTEM_ROOT}/neutralization_50"
SOURCE_DIR="${CFG_DIR}/npt_500ps"
JOBNAME="1IU4_EtOH50_md_${SLURM_JOB_ID}"
RUN_DIR="${SOURCE_DIR}/production_${SLURM_JOB_ID}"

cd "$SOURCE_DIR"
for file in equilibrated.cms npt_500ps.msj 1IU4_EtOH50_npt500_multisim.log; do
    [[ -s "$SOURCE_DIR/$file" ]] || {
        echo "Missing NPT input: $SOURCE_DIR/$file" >&2
        exit 1
    }
done
[[ -d "$SOURCE_DIR/npt_trj" ]] || {
    echo "Missing NPT trajectory: $SOURCE_DIR/npt_trj" >&2
    exit 1
}
for file in production.msj production.cfg settings.json; do
    [[ -s "$CFG_DIR/$file" ]] || {
        echo "Missing production input: $CFG_DIR/$file" >&2
        exit 1
    }
done
[[ -s "$SYSTEM_ROOT/prepare_1IU4_EtOH50.py" ]] || {
    echo "Missing preparation helper: $SYSTEM_ROOT/prepare_1IU4_EtOH50.py" >&2
    exit 1
}
[[ -s "$SYSTEM_ROOT/../expand_1IU4_EtOH20.py" ]] || {
    echo "Missing validated structure helper: $SYSTEM_ROOT/../expand_1IU4_EtOH20.py" >&2
    exit 1
}
grep -q '^Multisim completed\.' 1IU4_EtOH50_npt500_multisim.log || {
    echo "The 500 ps NPT continuation did not complete." >&2
    exit 1
}

export SYSTEM_ROOT CFG_DIR
"$SCHRODINGER/run" python3 - <<'PY'
import json
import os
import sys
from pathlib import Path

import numpy as np
from schrodinger.application.desmond.packages import traj
from schrodinger.utils import sea

system_root = Path(os.environ["SYSTEM_ROOT"])
config_dir = Path(os.environ["CFG_DIR"])
sys.path.insert(0, str(system_root))
import prepare_1IU4_EtOH50 as prep

st = prep.helper.read_full("equilibrated.cms")
groups = prep.check(st, prep.NEUTRAL, 0, 113710)
prep.helper.forcefield_audit("equilibrated.cms")
gap, diameter = prep.helper.solute_geometry(st, groups)
box = prep.helper.box(st)

settings = json.loads((config_dir / "settings.json").read_text())
cfg = sea.Map((config_dir / "production.cfg").read_text())
if float(settings["production_ns"]) != 200.0:
    raise SystemExit("Unexpected production duration in settings.json.")
if abs(float(cfg["time"].val) - 200000.0) > 1e-3:
    raise SystemExit("Unexpected production duration in production.cfg.")

frames = traj.read_traj("npt_trj")
if len(frames) < 20:
    raise SystemExit("Missing NPT continuation trajectory.")
times = np.asarray([frame.time for frame in frames], dtype=float)
volumes = np.asarray([abs(np.linalg.det(frame.box)) for frame in frames])
recent = times >= times[-1] - 200.0
t, v = times[recent], volumes[recent]
if len(v) < 15 or not (np.isfinite(v).all() and np.all(v > 0)):
    raise SystemExit("Insufficient valid volumes in the final 200 ps.")
mid = len(v) // 2
blocks = abs(v[:mid].mean() - v[mid:].mean()) / v.mean()
drift = abs(np.polyfit(t - t[0], v, 1)[0]) * (t[-1] - t[0]) / v.mean()
rho = st.total_weight * 1.66053906660 / v.mean()
rho_final = st.total_weight * 1.66053906660 / np.prod(box)
cutoff = float(settings["cutoff_A"])

print(f"Final box (A): {box.tolist()}", flush=True)
print(f"Final density (g/cm3): {rho_final:.4f}", flush=True)
print(f"Mean density, final 200 ps (g/cm3): {rho:.4f}", flush=True)
print(f"Volume block change: {blocks:.3%}; drift: {drift:.3%}", flush=True)
print(f"Solute image gap (A): {gap:.2f}; diameter (A): {diameter:.2f}", flush=True)
print(f"Cutoff (A): {cutoff:.2f}", flush=True)

if not 0.85 < rho < 1.05 or abs(rho_final - rho) / rho >= 0.02:
    raise SystemExit("The mixture density needs further NPT equilibration.")
if blocks >= 0.02 or drift >= 0.02:
    raise SystemExit("The box volume is still drifting; extend NPT equilibration.")
if gap <= max(cutoff, 15.0):
    raise SystemExit("Insufficient separation from a periodic solute image.")
print("PRODUCTION PREFLIGHT PASSED", flush=True)
PY

mkdir -p "$RUN_DIR"
cp equilibrated.cms "$CFG_DIR/production.msj" \
    "$CFG_DIR/production.cfg" "$CFG_DIR/settings.json" "$RUN_DIR/"
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
    -description "1IU4 50 percent v/v ethanol production MD" \
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
echo "Production trajectory: ${RUN_DIR}/production_trj"
