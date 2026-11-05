# TGase stability and substrate positioning in aqueous ethanol

Molecular dynamics of microbial transglutaminase (1IU4 numbering) with **Cbz-Gln-Gly (ZQG)** in water, nominal **20% v/v ethanol**, and nominal **50% v/v ethanol**. This repository organizes the simulation inputs, analysis code, numerical results and final figures from the Han Liu collaboration.

The study contains **one 200 ns production trajectory per condition**, or 600 ns in total. It examines protein fluctuations, substrate positioning, active-site geometry, local solvent composition and residue–substrate interaction energies. It does **not** establish enzyme activity, catalytic rates, binding free energies or replicated statistical differences between solvents.

![Representative enzyme–substrate conformations](figures/Figure_6_Representative_conformations.svg)

## Study design

| Condition | Water molecules | Ethanol molecules | Ethanol solvent number fraction | Total atoms | Production |
|---|---:|---:|---:|---:|---:|
| Water | 14,493 | 0 | 0 | 48,688 | 200 ns |
| Nominal 20% v/v ethanol | 27,885 | 2,112 | 0.070407 | 107,872 | 200 ns |
| Nominal 50% v/v ethanol | 19,415 | 5,584 | 0.223369 | 113,710 | 200 ns |

Each system contains one enzyme, one ZQG molecule and one chloride ion. The recorded production protocol uses Desmond/Schrödinger 2023-4, OPLS_2005, SPC water, NPT at 298 K and 1.01325 bar, a 9 Å cutoff and 2/2/6 fs multiple time steps. Production configurations record seed 2007; different solvent conditions are not independent replicate simulations of the same condition.

## What the results show

The following values describe the **50–200 ns portion of these particular trajectories**. “Near” means Cys64 SG–ZQG C28 ≤5 Å, a proximity criterion rather than proof of a catalytic pose.

| Condition | Near frames, dense sampling | Near frames, energy sampling | Cys64 NH–O22 proxy among sampled near frames |
|---|---:|---:|---:|
| Water | 368/751 (49.0%) | 73/150 (48.7%) | 65/73 (89.0%) |
| 20% ethanol | 22/751 (2.9%) | 4/150 (2.7%) | 0/4; very small subset |
| 50% ethanol | 695/751 (92.5%) | 136/150 (90.7%) | 0/136 |

The 20% trajectory samples substantial substrate displacement. The 50% trajectory retains ZQG near Cys64 much more often, but its near frames lack the monitored Cys64 backbone NH–O22 hydrogen-bond geometry. **Residence near the catalytic cysteine and productive substrate orientation are different observations.** Protein RMSD/RMSF provide a separate description of conformational stability and flexibility.

Dense and energy-sampled percentages differ because they use different frame grids. Energy analysis uses every fifth saved frame, approximately 1 ns apart. The dense late window contains 751 frames, while the sampled late window contains 150 frames. Neither count represents independent experimental or simulation replicates.

## Repository guide

| Location | Contents |
|---|---|
| [Methods](docs/METHODS.md) | Preparation, production settings, definitions, sampling and limitations |
| [Reproduction guide](docs/REPRODUCIBILITY.md) | CSV-only checks/plots, trajectory analysis, topology restoration and PyMOL |
| [Data dictionary](docs/DATA_DICTIONARY.md) | What the tables and principal columns mean |
| [Figure legends](docs/FIGURE_LEGENDS.md) | Supplied legends for the six final figures |
| [figures/](figures/) | Six final SVG figures and supplementary contact diagrams |
| [data/tables/](data/tables/) | Frame-level observables, pair energies and summary tables |
| [analysis/](analysis/) | Portable analysis scripts; unchanged source versions in `original/` |
| [simulation/](simulation/) | Preparation dependencies, historical Slurm scripts, production configurations and compressed input/output CMS files |
| [visualization/](visualization/) | Aligned PDBs, PyMOL script, rendered panels and compressed editable sessions |
| [data/representatives/](data/representatives/) | State-specific MAEGZ snapshots and their selection metadata |
| [provenance/](provenance/) | Source hashes, portability patches, completion evidence and external trajectory inventory |

## Quick reproduction from the published tables

No GPU or Schrödinger license is needed for these checks and plots:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-analysis.txt
python scripts/validate_publication.py
python scripts/reproduce_from_tables.py --output reproduced
```

The scripts verify data consistency and recompute the state/solvent and residue-pair summary tables before generating plots. They preserve the published data and final artwork. See the reproduction guide for the distinction between regenerated analysis figures and the supplied final layouts.

## Data availability and scope

Compressed production input/output CMS files, configurations, representative structures, analysis CSVs and final figures are included. **Full DTR trajectories and MD checkpoints remain on HPC and are not included in this Git repository.** Their locations and file sizes are recorded in [the external trajectory inventory](provenance/external_trajectory_inventory.json). Full trajectory-level reproduction requires those files and licensed Schrödinger software; the deposited tables support the CSV-only reproduction above.

Original source files were preserved. Portability edits change path handling only and are recorded as patches and SHA-256 hashes. The existing repository MIT license is retained; it does not license Schrödinger software or other third-party tools. See [LICENSE](LICENSE).
