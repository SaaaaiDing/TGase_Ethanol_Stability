#!/usr/bin/env python3
"""Re-evaluate a small set of explicit-solvent Desmond energy groups."""

import argparse
import csv
import math
import os
import tempfile
from pathlib import Path

from schrodinger.application.desmond.packages import energygroup, topo, traj
from schrodinger.utils import fileutils
from schrodinger.utils import sea


parser = argparse.ArgumentParser()
parser.add_argument("--system", choices=("ethanol20", "ethanol50"), required=True)
args = parser.parse_args()

analysis_dir = Path(os.environ["ANALYSIS_DIR"]).expanduser().resolve()
root = analysis_dir.parent
systems = {
    "ethanol20": (
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/1IU4_EtOH20_md_55774175-out.cms",
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/production_trj",
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/1IU4_EtOH20_md_55774175-out.cfg",
        2112,
        27885,
    ),
    "ethanol50": (
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/1IU4_EtOH50_md_55831542-out.cms",
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/production_trj",
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/1IU4_EtOH50_md_55831542-out.cfg",
        5584,
        19415,
    ),
}
cms_path, trj_path, cfg_path, ethanol_count, water_count = systems[args.system]
for path in (cms_path, trj_path, cfg_path):
    if not path.exists():
        raise FileNotFoundError(path)
cfg = sea.Map(cfg_path.read_text())
if "ORIG_CFG" not in cfg:
    raise RuntimeError(f"ORIG_CFG is missing from {cfg_path}")

_, cms = topo.read_cms(str(cms_path))
ct = cms.fsys_ct
protein_aids = set(cms.select_atom("protein"))
groups = {}
for name, number in (("CYS64", 64), ("VAL65", 65)):
    groups[name] = {
        aid for aid in protein_aids
        if ct.atom[aid].chain.strip() == "A"
        and ct.atom[aid].resnum == number
    }
groups["ZQG"] = set(cms.select_atom("res.ptype LIG1"))
groups["ethanol"] = set(cms.select_atom("res.ptype UNK"))
groups["water"] = set(cms.select_atom("res.ptype SPC"))
if len(groups["ZQG"]) != 42:
    raise RuntimeError("ZQG selection is not 42 atoms")
if len(groups["ethanol"]) != ethanol_count * 9:
    raise RuntimeError("Ethanol atom count is inconsistent")
if len(groups["water"]) != water_count * 3:
    raise RuntimeError("Water atom count is inconsistent")
if not all(groups.values()):
    raise RuntimeError("One or more energy groups are empty")
names = list(groups)
if sum(len(group) for group in groups.values()) != len(set.union(*groups.values())):
    raise RuntimeError("Energy groups overlap")
print(f"System: {args.system}", flush=True)
print(f"Groups: {[(name, len(groups[name])) for name in names]}", flush=True)
frames = traj.read_traj(str(trj_path))
selected = [frame for frame in frames if 50.0 <= frame.time / 1000.0 < 55.0]
selected = selected[::5]
if not 4 <= len(selected) <= 7:
    raise RuntimeError(f"Unexpected number of selected frames: {len(selected)}")
print(
    f"Selected frame times (ns): {[round(frame.time / 1000, 3) for frame in selected]}",
    flush=True,
)


def report_backend_files(previous):
    temp_root = Path(fileutils.get_directory_path(fileutils.TEMP))
    new_dirs = [
        path for path in temp_root.iterdir()
        if path not in previous and path.is_dir() and (path / "energy.cfg").exists()
    ]
    for folder in sorted(new_dirs):
        print(f"Retained backend directory: {folder}", flush=True)
        print(f"Backend files: {sorted(path.name for path in folder.iterdir())}", flush=True)
        for name in ("energy.cfg", "output.engrp"):
            path = folder / name
            if not path.is_file():
                print(f"Missing backend file: {path}", flush=True)
                continue
            print(f"{name}: {path.stat().st_size} bytes", flush=True)
            if name == "output.engrp":
                lines = path.read_text(errors="replace").splitlines()
                print(f"First backend energy lines: {lines[:12]}", flush=True)
                print(f"Last backend energy lines: {lines[-8:]}", flush=True)
            else:
                lines = path.read_text(errors="replace").splitlines()
                relevant = [
                    line for line in lines if "vrun.frameset" in line
                    or "energy_groups.first" in line
                    or "energy_groups.interval" in line
                ]
                print(f"Energy frame settings: {relevant[:12]}", flush=True)
        for path in sorted(folder.glob("*.log")):
            tail = path.read_text(errors="replace").splitlines()[-20:]
            print(f"Last lines of {path.name}: {tail}", flush=True)


options = set(energygroup.InteractionEnergy.DEFAULT_OPTIONS)
os.environ["SCHRODINGER_JOB_DEBUG"] = "1"
with tempfile.TemporaryDirectory(prefix="energy_pilot_frames_") as temp_dir:
    subset_path = Path(temp_dir) / "selected_trj"
    traj.write_traj(selected, str(subset_path))
    temp_root = Path(fileutils.get_directory_path(fileutils.TEMP))
    previous = set(temp_root.iterdir())
    try:
        result = energygroup.get_energies(
            str(cfg_path),
            str(subset_path),
            cms,
            energygroup.SliceParams(0, 0, 0),
            [groups[name] for name in names],
            options=options,
        )
    except Exception:
        report_backend_files(previous)
        raise
    if len(result) != 2 or not result[0]:
        report_backend_files(previous)
        raise RuntimeError(
            "Desmond returned no parsed energy frames; inspect the retained backend "
            "directory shown above before another run"
        )
    times, components = result
    if len(times) != len(components):
        raise RuntimeError("Energy time and component counts disagree")
print(f"Computed {len(times)} frames", flush=True)

output = analysis_dir / f"{args.system}_energy_pilot_50_55ns.csv"
rows = []
for time_ps, frame in zip(times, components):
    for residue in ("CYS64", "VAL65"):
        for partner in ("ZQG", "ethanol", "water"):
            pair = frozenset((names.index(residue), names.index(partner)))
            if pair not in frame:
                raise RuntimeError(f"Missing energy pair: {residue}, {partner}")
            energy = frame[pair]
            electrostatic = float(energy.elec)
            vdw = float(energy.vdw)
            if not (math.isfinite(electrostatic) and math.isfinite(vdw)):
                raise RuntimeError(f"Nonfinite energy: {residue}, {partner}")
            rows.append({
                "system": args.system,
                "time_ns": float(time_ps) / 1000.0,
                "residue": residue,
                "partner": partner,
                "electrostatic_kcal_per_mol": electrostatic,
                "vdw_kcal_per_mol": vdw,
                "nonbonded_sum_kcal_per_mol": electrostatic + vdw,
            })

with output.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(f"Saved: {output}", flush=True)
print("Pair energies are force-field interaction terms, not binding free energies.")
print("Ethanol and water totals include all solvent molecules in their groups.")
