# Reproduction guide

## Published tables and analysis figures

From the repository root:

```bash
python -m pip install -r requirements-analysis.txt
python scripts/validate_publication.py
python scripts/reproduce_from_tables.py --output reproduced
```

Validation checks the imported file hashes, compressed-file recovery, Python syntax, CSV frame counts/times, dense near-state counts, sampled frame mapping, pair-energy additivity and reproduction of the two main summary tables. It does not rerun MD or certify physical convergence.

The reproduction command writes a fresh output directory, rebuilds the comparison tables, and regenerates analysis plots for RMSD/proximity, reactive geometry, state/solvation and pair energies. These are the analysis-stage layouts. The six authoritative final SVG layouts are preserved unchanged in `figures/`; some labels, colors and formatting were edited between the analysis-stage figures and the final website versions. The deposited RMSF CSV and final RMSF figure are included; trajectory-based recalculation of RMSF requires the raw data described below.

## Full trajectory-level analysis

Use the Schrödinger 2023-4 Python environment for the vendor `topo`, `traj`, `analysis` and `energygroup` APIs. NumPy/Matplotlib are also used by the scripts. These APIs and the MD engine are not included or installed by the repository's pip requirements.

First restore the deposited CMS/configuration layout to a new directory:

```bash
python scripts/restore_simulation_files.py --output-root /path/to/restored_1IU4
export TGASE_DATA_ROOT=/path/to/restored_1IU4
export ANALYSIS_DIR=/path/to/new_analysis
mkdir -p "$ANALYSIS_DIR"
```

Obtain the DTR directories listed in `provenance/external_trajectory_inventory.json` from the project maintainers and place them at the recorded paths relative to `TGASE_DATA_ROOT`. They are not embedded in this Git repository. Keep all DTR component files together. The deposited output CMS files correspond to these trajectories; do not substitute unrelated topologies.

Run the main observable stages in order, using a fresh analysis directory:

```bash
"$SCHRODINGER/run" python3 analysis/step01_rmsd.py
"$SCHRODINGER/run" python3 analysis/step02_active_site_proximity.py
"$SCHRODINGER/run" python3 analysis/step04_reactive_geometry.py
"$SCHRODINGER/run" python3 analysis/step06_active_site_solvent.py
"$SCHRODINGER/run" python3 analysis/step07_state_conditioned_contacts.py
"$SCHRODINGER/run" python3 analysis/step08_carbonyl_backbone_hbonds.py
"$SCHRODINGER/run" python3 analysis/step09_rmsf_representatives.py
for system in water ethanol20 ethanol50; do
  "$SCHRODINGER/run" python3 analysis/step12_desmond_pair_energies.py --system "$system" --stride 5 --batch-size 25
done
python analysis/make_1IU4_md_comparison.py --input-dir "$ANALYSIS_DIR" --output-dir "$ANALYSIS_DIR/comparison"
"$SCHRODINGER/run" python3 analysis/step13_representative_atom_contacts.py --analysis-dir "$ANALYSIS_DIR"
```

Steps 10 and 11 are historical energy preflight/pilot tools. They are retained for provenance and are not required when reproducing the final stride-5 analysis. Portable versions of the analysis scripts accept the two environment variables above; their original versions and path-only patches are preserved. Historical tables may retain original absolute snapshot paths; `data/tables/representative_frames_portable.csv` supplies repository-relative paths.

## Figures and PyMOL scenes

```bash
cd visualization
pymol -cq render_1IU4_pymol.py
```

The script expects the three aligned PDBs next to it and writes `pymol_figures/`. It checks displayed distances before rendering. Pillow is optional for assembling its combined panel. Pre-rendered PNGs and gzip-compressed editable PSE sessions are deposited. Decompress a session with `gzip -dk sessions/water_contacts.pse.gz`, then open the resulting `.pse` in PyMOL.

The archive README reflects an earlier environment in which PyMOL was unavailable. Later rendered PNGs/PSE sessions are present in the supplied source files and are included here. Publication validation does not claim to have rerendered those sessions or reproduced the final Figure 6 layout.

## Simulation preparation and reruns

`simulation/preparation/` contains the 20% expansion/workflow helpers and the 50% neutralization/preparation helper. `simulation/raw_layout/` retains actual packing, relaxation and production configurations in their original relative directory layout, including the extra 500 ps NPT stage for 50% ethanol. Production input/output CMS files are gzip-compressed without changing their decompressed contents.

Historical Slurm files retain their original accounts, partitions, GPU requests and site paths as provenance. They are not universal launchers. Adapt site-specific paths/resources in a separate copy before submission. This repository performs no scheduler submission automatically. Re-running a production configuration is a new simulation, not a bitwise reproduction or continuation of the deposited trajectory; restart checkpoints are not included.

## File integrity

`provenance/source_manifest.json` records source and published SHA-256 hashes and transformations. All copied final figures and numerical tables remain byte-identical to their source files. Compressed CMS/PSE contents recover their source hashes. The original source directories remain unchanged. A publication-time validation report is recorded in `provenance/validation_report.json`.
