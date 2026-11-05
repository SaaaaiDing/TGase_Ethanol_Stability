#!/bin/bash
#SBATCH --job-name=1IU4_EtOH50_pack
#SBATCH --account=yjz5549_cr_default
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --export=NONE
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

module purge
module load schrodinger/2023.4

SOURCE_DIR="$SLURM_SUBMIT_DIR"
RUN_DIR="${SOURCE_DIR}/packing_${SLURM_JOB_ID}"

[[ -s "$SOURCE_DIR/solute_1IU4.maegz" ]] || {
    echo "Missing solute_1IU4.maegz" >&2
    exit 1
}

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
cp "$SOURCE_DIR/solute_1IU4.maegz" .

unset SLURM_STEP_ID SLURM_STEPID SLURM_PROCID SLURM_LOCALID SLURM_GTIDS
unset OMPI_MCA_plm_slurm_args PMI_RANK PMI_SIZE PMI_FD PMIX_RANK PMIX_NAMESPACE
export SLURM_EXPORT_ENV=ALL
export SOURCE_DIR

echo "Run directory: $RUN_DIR"
echo "Start time: $(date)"

"$SCHRODINGER/run" python3 - <<'PY'
import os
import sys
from pathlib import Path

import numpy as np
from schrodinger import structure

source_dir = Path(os.environ["SOURCE_DIR"]).resolve()
sys.path.insert(0, str(source_dir.parent))
import expand_1IU4_EtOH20 as helper

source = source_dir.parent / "disordered_system_1" / "disordered_system_1_system-out.cms"
st = helper.read_full(source)
groups = helper.groups(st)

expected = {"water": 9296, "ethanol": 704, "chloride": 0, "solute": 2}
counts = {name: len(molecules) for name, molecules in groups.items()}
if counts != expected:
    raise SystemExit(f"Unexpected solvent template source: {counts}")

coordinates = helper.whole_centered(st, groups)

with structure.StructureWriter("solvent_templates_50.maegz") as writer:
    for kind in ("ethanol", "water"):
        ids = groups[kind][0]
        centered = coordinates.copy()
        indices = np.asarray(ids) - 1
        centered[indices] -= centered[indices].mean(axis=0)

        template = helper.clean_copy(st, ids, centered)
        template.title = kind
        writer.append(template)
        print(f"{kind} template atoms: {template.atom_total}")

print("Saved: solvent_templates_50.maegz")
PY

"$SCHRODINGER/run" \
    disordered_system_builder_gui_dir/disordered_system_builder_driver.py \
    -molecules 25000 \
    -composition 5584:19416 \
    -substrate solute_1IU4.maegz \
    -substrate_mode immersed \
    -pbc new \
    -pack \
    -density 0.500 \
    -scale_vdw 0.50 \
    -obey density \
    -tries_per_mol 50 \
    -tries_per_dc 20 \
    -seed 3456 \
    -preserve_res_info \
    -water_fftype SPC \
    -split_components \
    -forcefield OPLS_2005 \
    solvent_templates_50.maegz \
    mix50 > packing.log 2>&1

[[ -s mix50_system-out.cms ]] || {
    echo "Packed CMS is missing." >&2
    exit 1
}

"$SCHRODINGER/run" python3 - <<'PY'
import os
import sys
from pathlib import Path

source_dir = Path(os.environ["SOURCE_DIR"]).resolve()
sys.path.insert(0, str(source_dir.parent))
import expand_1IU4_EtOH20 as helper

st = helper.read_full("mix50_system-out.cms")
groups = helper.groups(st)
counts = {name: len(molecules) for name, molecules in groups.items()}

expected = {"water": 19416, "ethanol": 5584, "chloride": 0, "solute": 2}
print("Packed composition:", counts)
print("Atoms:", st.atom_total)
print("Formal charge:", st.formal_charge)
print("Initial box (A):", helper.box(st).tolist())

if counts != expected or st.formal_charge != 1 or st.atom_total != 113712:
    raise SystemExit("Packed system validation failed.")

print("PACKING PASSED")
PY

echo "End time: $(date)"
echo "Packed system: $RUN_DIR/mix50_system-out.cms"
