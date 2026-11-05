#!/usr/bin/env python3
"""Sample explicit-solvent Desmond pair energies by ligand proximity state.

Run with Schrödinger 2023.4: $SCHRODINGER/run python3 step12_desmond_pair_energies.py --system ethanol50
These are nonbonded force-field interaction terms, not MM/GBSA or binding free energies.
"""

import argparse
import csv
import math
import os
import statistics
import tempfile
from pathlib import Path

from schrodinger.application.desmond.packages import energygroup, topo, traj
from schrodinger.utils import sea


SYSTEMS = {
    "water": (
        "1IU4_0/desmond_md_job_4/1IU4_md_batch-out.cms",
        "1IU4_0/desmond_md_job_4/1IU4_md_batch_trj",
        "1IU4_0/desmond_md_job_4/1IU4_md_batch-out.cfg",
        0, 14493,
    ),
    "ethanol20": (
        "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/1IU4_EtOH20_md_55774175-out.cms",
        "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/production_trj",
        "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/1IU4_EtOH20_md_55774175-out.cfg",
        2112, 27885,
    ),
    "ethanol50": (
        "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/1IU4_EtOH50_md_55831542-out.cms",
        "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/production_trj",
        "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/1IU4_EtOH50_md_55831542-out.cfg",
        5584, 19415,
    ),
}

RESIDUES = {
    62: "TYR", 64: "CYS", 65: "VAL", 75: "TYR", 208: "ARG",
    238: "ARG", 252: "VAL", 253: "ASN", 254: "PHE", 255: "ASP",
    274: "HIS", 276: "ASN", 277: "HIS", 278: "TYR",
}

FIELDS = (
    "system", "frame_index", "time_ns", "period_ns", "SG_C28_state",
    "SG_C28_A", "Cys64_NH_O22_proxy", "Val65_NH_O22_proxy",
    "residue", "partner", "electrostatic_kcal_per_mol",
    "vdw_kcal_per_mol", "nonbonded_sum_kcal_per_mol",
)


def read_reference(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1002:
        raise RuntimeError(f"Expected 1002 reference frames in {path}; found {len(rows)}")
    return rows


def build_groups(cms, ethanol_count, water_count):
    ct = cms.fsys_ct
    protein = set(cms.select_atom("protein"))
    groups = {}
    for number, expected in RESIDUES.items():
        aids = {
            aid for aid in protein
            if ct.atom[aid].chain.strip() == "A" and ct.atom[aid].resnum == number
        }
        observed = {ct.atom[aid].pdbres.strip() for aid in aids}
        if not aids or observed != {expected}:
            raise RuntimeError(f"Residue A:{number}: expected {expected}, got {observed}")
        groups[f"{expected}{number}"] = aids
    groups["ZQG"] = set(cms.select_atom("res.ptype LIG1"))
    if len(groups["ZQG"]) != 42:
        raise RuntimeError("Expected 42 ZQG atoms")
    if ethanol_count:
        groups["ethanol"] = set(cms.select_atom("res.ptype UNK"))
    groups["water"] = set(cms.select_atom("res.ptype SPC"))
    if len(groups.get("ethanol", ())) != 9 * ethanol_count:
        raise RuntimeError("Unexpected ethanol atom count")
    if len(groups["water"]) != 3 * water_count:
        raise RuntimeError("Unexpected water atom count")
    if sum(map(len, groups.values())) != len(set.union(*groups.values())):
        raise RuntimeError("Energy groups overlap")
    return groups


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", required=True, choices=SYSTEMS)
    parser.add_argument("--stride", type=int, default=5,
                        help="Trajectory frame stride; 5 is about 1 ns")
    parser.add_argument("--batch-size", type=int, default=25)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.stride < 1 or args.batch_size < 1:
        parser.error("Stride and batch size must be positive")

    analysis_dir = Path(os.environ["ANALYSIS_DIR"]).expanduser().resolve()
    paths = SYSTEMS[args.system]
    cms_path, trj_path, cfg_path = (analysis_dir.parent / x for x in paths[:3])
    ethanol_count, water_count = paths[3:]
    ref_path = analysis_dir / f"{args.system}_carbonyl_backbone_hbonds.csv"
    for path in (cms_path, trj_path, cfg_path, ref_path):
        if not path.exists():
            raise FileNotFoundError(path)
    if "ORIG_CFG" not in sea.Map(cfg_path.read_text()):
        raise RuntimeError(f"ORIG_CFG missing from {cfg_path}")

    _, cms = topo.read_cms(str(cms_path))
    groups = build_groups(cms, ethanol_count, water_count)
    names = list(groups)
    residue_names = names[:len(RESIDUES)]
    partners = [name for name in ("ZQG", "ethanol", "water") if name in groups]
    frames = traj.read_traj(str(trj_path))
    reference = read_reference(ref_path)
    if len(frames) != len(reference):
        raise RuntimeError("Trajectory and reference frame counts differ")
    selected_indices = list(range(0, len(frames), args.stride))
    for i in selected_indices:
        if abs(frames[i].time / 1000 - float(reference[i]["time_ns"])) > 1e-5:
            raise RuntimeError(f"Frame {i} time disagrees with reference")

    print(f"System: {args.system}; frames: {len(selected_indices)}; "
          f"groups: {[(k, len(v)) for k, v in groups.items()]}", flush=True)
    print(f"Pairs per frame: {len(residue_names) * len(partners)}", flush=True)
    if args.dry_run:
        return

    output = analysis_dir / f"{args.system}_pair_energies_stride{args.stride}.csv"
    summary_path = analysis_dir / f"{args.system}_ligand_pair_summary_stride{args.stride}.csv"
    if output.exists() or summary_path.exists():
        raise FileExistsError(f"Output exists: {output} or {summary_path}")

    name_to_index = {name: i for i, name in enumerate(names)}
    options = set(energygroup.InteractionEnergy.DEFAULT_OPTIONS)
    all_rows = []
    for start in range(0, len(selected_indices), args.batch_size):
        indices = selected_indices[start:start + args.batch_size]
        with tempfile.TemporaryDirectory(prefix=f"{args.system}_energy_") as temp_dir:
            subset = str(Path(temp_dir) / "selected_trj")
            traj.write_traj([frames[i] for i in indices], subset)
            result = energygroup.get_energies(
                str(cfg_path), subset, cms,
                energygroup.SliceParams(0, 0, 0),
                [groups[name] for name in names], options=options,
            )
            if len(result) != 2:
                raise RuntimeError(f"Backend returned no pair energies for frames {indices}")
            times, components = result
            if len(times) != len(indices) or len(components) != len(indices):
                raise RuntimeError(f"Backend energy frame count differs for {indices}")
            for i, time_ps, energy_frame in zip(indices, times, components):
                expected_ns = frames[i].time / 1000.0
                if abs(time_ps / 1000.0 - expected_ns) > 0.01:
                    raise RuntimeError(f"Backend time mismatch at frame {i}")
                ref = reference[i]
                distance = float(ref["SG_C28_A"])
                common = {
                    "system": args.system,
                    "frame_index": i,
                    "time_ns": expected_ns,
                    "period_ns": "0-50" if expected_ns < 50 else "50-200",
                    "SG_C28_state": "near" if distance <= 5 else "far",
                    "SG_C28_A": distance,
                    "Cys64_NH_O22_proxy": int(ref["Cys64_NH_O22_proxy"]),
                    "Val65_NH_O22_proxy": int(ref["Val65_NH_O22_proxy"]),
                }
                for residue in residue_names:
                    for partner in partners:
                        pair = frozenset((name_to_index[residue], name_to_index[partner]))
                        if pair not in energy_frame:
                            raise RuntimeError(f"Backend missing {residue}-{partner} at {i}")
                        energy = energy_frame[pair]
                        elec, vdw = float(energy.elec), float(energy.vdw)
                        if not all(map(math.isfinite, (elec, vdw))):
                            raise RuntimeError(f"Nonfinite energy for {residue}-{partner}")
                        all_rows.append({
                            **common, "residue": residue, "partner": partner,
                            "electrostatic_kcal_per_mol": elec,
                            "vdw_kcal_per_mol": vdw,
                            "nonbonded_sum_kcal_per_mol": elec + vdw,
                        })
        print(f"Computed frames {indices[0]}..{indices[-1]}; "
              f"{len(all_rows)} energy rows total", flush=True)

    if len(all_rows) != len(selected_indices) * len(residue_names) * len(partners):
        raise RuntimeError("Wrong number of pair-energy rows")
    write_csv(output, all_rows, FIELDS)

    summary = []
    for residue in residue_names:
        for period in ("0-50", "50-200"):
            for state in ("near", "far"):
                values = [
                    r for r in all_rows if r["residue"] == residue
                    and r["partner"] == "ZQG" and r["period_ns"] == period
                    and r["SG_C28_state"] == state
                ]
                if not values:
                    continue
                energies = [r["nonbonded_sum_kcal_per_mol"] for r in values]
                summary.append({
                    "system": args.system, "residue": residue,
                    "period_ns": period, "SG_C28_state": state,
                    "sampled_frames": len(values),
                    "interpret_with_caution": int(len(values) < 10),
                    "mean_pair_kcal_per_mol": statistics.mean(energies),
                    "median_pair_kcal_per_mol": statistics.median(energies),
                    "sd_pair_kcal_per_mol": statistics.stdev(energies)
                    if len(values) > 1 else "",
                })
    write_csv(summary_path, summary, list(summary[0]))
    print(f"Saved: {output}", flush=True)
    print(f"Saved: {summary_path}", flush=True)
    print("Residue-ZQG terms are pair interactions, not binding free energies or MM/GBSA decomposition.")
    print("Water and ethanol totals contain every molecule in each group; do not compare raw totals across systems.")


if __name__ == "__main__":
    main()
