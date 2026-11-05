#!/usr/bin/env python3
"""Count water and ethanol hydroxyl oxygens around Cys64 SG in MD."""

import csv
import os
from pathlib import Path

import numpy as np
from schrodinger.application.desmond.packages import topo, traj


analysis_dir = Path(os.environ["ANALYSIS_DIR"]).expanduser().resolve()
root = Path(os.environ.get("TGASE_DATA_ROOT", str(analysis_dir.parent))).expanduser().resolve()
radii = (4, 6, 8)
systems = {
    "water": (
        root / "1IU4_0/desmond_md_job_4/1IU4_md_batch-out.cms",
        root / "1IU4_0/desmond_md_job_4/1IU4_md_batch_trj",
        14493, 0,
    ),
    "ethanol20": (
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/1IU4_EtOH20_md_55774175-out.cms",
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/production_trj",
        27885, 2112,
    ),
    "ethanol50": (
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/1IU4_EtOH50_md_55831542-out.cms",
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/production_trj",
        19415, 5584,
    ),
}


def solvent_oxygen_aids(ct):
    water_aids, ethanol_aids = [], []
    for molecule in ct.molecule:
        size = len(molecule.atom)
        if size not in (3, 9):
            continue
        atoms = list(molecule.atom)
        elements = sorted(atom.element for atom in atoms)
        oxygens = [atom.index for atom in atoms if atom.element == "O"]
        if elements == ["H", "H", "O"]:
            water_aids.extend(oxygens)
        elif elements == ["C", "C", "H", "H", "H", "H", "H", "H", "O"]:
            ethanol_aids.extend(oxygens)
    return water_aids, ethanol_aids


summary = []
for label, (cms_path, trj_path, expected_water, expected_ethanol) in systems.items():
    if not cms_path.is_file() or not trj_path.exists():
        raise FileNotFoundError(f"Missing input: {cms_path} or {trj_path}")
    _, cms = topo.read_cms(str(cms_path))
    ct = cms.fsys_ct
    sg = ct.atom[541]
    if (sg.pdbres.strip(), sg.resnum, sg.pdbname.strip()) != ("CYS", 64, "SG"):
        raise RuntimeError(f"{label}: AID 541 is not Cys64 SG")

    water_aids, ethanol_aids = solvent_oxygen_aids(ct)
    if (len(water_aids), len(ethanol_aids)) != (expected_water, expected_ethanol):
        raise RuntimeError(
            f"{label}: solvent count mismatch: "
            f"water={len(water_aids)}, ethanol={len(ethanol_aids)}"
        )
    water_gids = topo.aids2gids(cms, water_aids, include_pseudoatoms=False)
    ethanol_gids = topo.aids2gids(cms, ethanol_aids, include_pseudoatoms=False)
    sg_gid = topo.aids2gids(cms, [541], include_pseudoatoms=False)[0]
    bulk_ethanol_fraction = expected_ethanol / (expected_water + expected_ethanol)

    rows = []
    for frame in traj.read_traj(str(trj_path)):
        box = np.asarray(frame.box, dtype=float)
        if box.shape != (3, 3) or not np.allclose(box, np.diag(np.diag(box)), atol=1e-3):
            raise RuntimeError(f"{label}: an orthorhombic box is required")
        lengths = np.diag(box)
        if np.any(lengths <= 0):
            raise RuntimeError(f"{label}: invalid box lengths")
        center = np.asarray(frame.pos(sg_gid), dtype=float).reshape(3)
        row = {"time_ns": frame.time / 1000.0}

        for kind, gids in (("water", water_gids), ("ethanol", ethanol_gids)):
            if gids:
                xyz = np.asarray(frame.pos(gids), dtype=float)
                delta = xyz - center
                delta -= np.rint(delta / lengths) * lengths
                distance_squared = np.einsum("ij,ij->i", delta, delta)
            else:
                distance_squared = np.empty(0, dtype=float)
            for radius in radii:
                row[f"{kind}_O_{radius}A"] = int(
                    np.count_nonzero(distance_squared <= radius * radius)
                )
        rows.append(row)

    if len(rows) != 1002:
        raise RuntimeError(f"{label}: expected 1002 frames, found {len(rows)}")
    times = np.asarray([row["time_ns"] for row in rows])
    if not np.all(np.diff(times) > 0) or not 199.9 < times[-1] < 200.1:
        raise RuntimeError(f"{label}: invalid trajectory times")

    frame_path = analysis_dir / f"{label}_Cys64_solvent_shell.csv"
    with frame_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"\n{label}: {len(rows)} frames; {expected_water} water, "
        f"{expected_ethanol} ethanol; bulk ethanol mole fraction "
        f"{bulk_ethanol_fraction:.4f}"
    )
    for start, stop in ((0, 50), (50, 200.01)):
        selected = [row for row in rows if start <= row["time_ns"] < stop]
        for radius in radii:
            w = np.asarray([row[f"water_O_{radius}A"] for row in selected])
            e = np.asarray([row[f"ethanol_O_{radius}A"] for row in selected])
            total_oxygen = int(w.sum() + e.sum())
            local_fraction = float(e.sum() / total_oxygen) if total_oxygen else float("nan")
            enrichment = (
                local_fraction / bulk_ethanol_fraction
                if bulk_ethanol_fraction > 0 and total_oxygen else float("nan")
            )
            entry = {
                "system": label,
                "start_ns": start,
                "end_ns": min(stop, 200),
                "radius_A": radius,
                "frames": len(selected),
                "mean_water_O_per_frame": float(w.mean()),
                "mean_ethanol_O_per_frame": float(e.mean()),
                "fraction_frames_with_ethanol_O": float(np.mean(e > 0)),
                "bulk_ethanol_mole_fraction": bulk_ethanol_fraction,
                "local_ethanol_O_fraction": local_fraction,
                "local_to_bulk_ethanol_fraction_ratio": enrichment,
            }
            summary.append(entry)
            if radius == 6:
                print(
                    f"  {start:3d}-{min(stop, 200):3.0f} ns, 6 A: "
                    f"mean water O={w.mean():.2f}, ethanol O={e.mean():.2f}; "
                    f"local ethanol fraction={local_fraction:.3f}; "
                    f"local/bulk ratio={enrichment:.2f}"
                )
    print(f"Saved: {frame_path}")

summary_path = analysis_dir / "Cys64_solvent_shell_summary.csv"
with summary_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
    writer.writeheader()
    writer.writerows(summary)
print(f"Saved: {summary_path}")
