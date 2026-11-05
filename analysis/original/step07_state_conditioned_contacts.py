#!/usr/bin/env python3
"""Analyze ZQG-protein contacts and solvent around Cys64 by reactive proximity."""

import csv
import os
from pathlib import Path

import numpy as np
from schrodinger.application.desmond.packages import topo, traj


analysis_dir = Path(os.environ["ANALYSIS_DIR"]).expanduser().resolve()
root = analysis_dir.parent
contact_cutoff_A = 4.5
near_cutoff_A = 5.0

systems = {
    "water": (
        root / "1IU4_0/desmond_md_job_4/1IU4_md_batch-out.cms",
        root / "1IU4_0/desmond_md_job_4/1IU4_md_batch_trj",
    ),
    "ethanol20": (
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/1IU4_EtOH20_md_55774175-out.cms",
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/production_trj",
    ),
    "ethanol50": (
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/1IU4_EtOH50_md_55831542-out.cms",
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/production_trj",
    ),
}


def load_csv(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1002:
        raise RuntimeError(f"Expected 1002 frames in {path}, found {len(rows)}")
    return rows


for label, (cms_path, trj_path) in systems.items():
    if not cms_path.is_file() or not trj_path.exists():
        raise FileNotFoundError(f"Missing input: {cms_path} or {trj_path}")
    for suffix in ("reactive_geometry", "Cys64_solvent_shell"):
        input_csv = analysis_dir / f"{label}_{suffix}.csv"
        if not input_csv.is_file():
            raise FileNotFoundError(f"Missing input: {input_csv}")

contact_summary = []
solvent_summary = []

for label, (cms_path, trj_path) in systems.items():

    geometry = load_csv(analysis_dir / f"{label}_reactive_geometry.csv")
    solvent = load_csv(analysis_dir / f"{label}_Cys64_solvent_shell.csv")
    geometry_times = np.array([float(row["time_ns"]) for row in geometry])
    solvent_times = np.array([float(row["time_ns"]) for row in solvent])
    if not np.allclose(geometry_times, solvent_times, atol=1e-5, rtol=0):
        raise RuntimeError(f"{label}: geometry and solvent frame times differ")
    near = np.array([
        float(row["Cys64_SG_to_ZQG_C28_A"]) <= near_cutoff_A
        for row in geometry
    ])

    _, cms = topo.read_cms(str(cms_path))
    ct = cms.fsys_ct
    protein_aids = [
        aid for aid in cms.select_atom("protein")
        if ct.atom[aid].element != "H"
    ]
    ligand_aids = [
        atom.index for atom in ct.atom
        if atom.pdbres.strip() == "LIG1" and atom.element != "H"
    ]
    if len(ligand_aids) != 24 or not 2000 <= len(protein_aids) <= 4000:
        raise RuntimeError(
            f"{label}: unexpected heavy atom counts: "
            f"protein={len(protein_aids)}, ligand={len(ligand_aids)}"
        )

    residue_keys = [
        (ct.atom[aid].chain.strip(), ct.atom[aid].resnum,
         ct.atom[aid].inscode.strip(), ct.atom[aid].pdbres.strip())
        for aid in protein_aids
    ]
    unique_residues = sorted(set(residue_keys), key=lambda x: (x[0], x[1], x[2]))
    residue_lookup = {key: index for index, key in enumerate(unique_residues)}
    residue_index = np.array([residue_lookup[key] for key in residue_keys])
    if len(unique_residues) != 331:
        raise RuntimeError(f"{label}: expected 331 protein residues, found {len(unique_residues)}")

    protein_gids = topo.aids2gids(cms, protein_aids, include_pseudoatoms=False)
    ligand_gids = topo.aids2gids(cms, ligand_aids, include_pseudoatoms=False)
    tr = traj.read_traj(str(trj_path))
    if len(tr) != 1002:
        raise RuntimeError(f"{label}: expected 1002 trajectory frames")

    contacts = np.zeros((len(tr), len(unique_residues)), dtype=bool)
    for frame_number, frame in enumerate(tr):
        if abs(frame.time / 1000.0 - geometry_times[frame_number]) > 1e-5:
            raise RuntimeError(f"{label}: trajectory and geometry frame times differ")
        box = np.asarray(frame.box, dtype=float)
        if box.shape != (3, 3) or not np.allclose(box, np.diag(np.diag(box)), atol=1e-3):
            raise RuntimeError(f"{label}: an orthorhombic box is required")
        lengths = np.diag(box)
        if np.any(lengths <= 0):
            raise RuntimeError(f"{label}: invalid box lengths")

        protein_xyz = np.asarray(frame.pos(protein_gids), dtype=float)
        ligand_xyz = np.asarray(frame.pos(ligand_gids), dtype=float)
        delta = protein_xyz[:, None, :] - ligand_xyz[None, :, :]
        delta -= np.rint(delta / lengths) * lengths
        minimum_squared = np.min(np.einsum("ijk,ijk->ij", delta, delta), axis=1)
        contacted_residues = residue_index[minimum_squared <= contact_cutoff_A ** 2]
        contacts[frame_number, contacted_residues] = True

        if frame_number and frame_number % 250 == 0:
            print(f"{label}: processed {frame_number}/{len(tr)} frames", flush=True)

    periods = (("0-50", 0.0, 50.0), ("50-200", 50.0, 200.01))
    for period_name, begin, end in periods:
        in_period = (geometry_times >= begin) & (geometry_times < end)
        for state_name, state_mask in (("near", near), ("far", ~near)):
            selection = in_period & state_mask
            nframes = int(selection.sum())
            if not nframes:
                continue

            occupancy = np.mean(contacts[selection], axis=0)
            for residue, frequency in zip(unique_residues, occupancy):
                if frequency == 0:
                    continue
                chain, number, inscode, name = residue
                contact_summary.append({
                    "system": label,
                    "period_ns": period_name,
                    "SG_C28_state": state_name,
                    "frames_in_state": nframes,
                    "chain": chain,
                    "residue_number": number,
                    "insertion_code": inscode,
                    "residue_name": name,
                    "contact_fraction_4p5A": float(frequency),
                })

            w = np.array([
                int(solvent[i]["water_O_6A"]) for i in np.flatnonzero(selection)
            ])
            e = np.array([
                int(solvent[i]["ethanol_O_6A"]) for i in np.flatnonzero(selection)
            ])
            oxygen_total = int(w.sum() + e.sum())
            local_fraction = float(e.sum() / oxygen_total) if oxygen_total else float("nan")
            solvent_summary.append({
                "system": label,
                "period_ns": period_name,
                "SG_C28_state": state_name,
                "frames_in_state": nframes,
                "fraction_of_period": nframes / int(in_period.sum()),
                "mean_water_O_within_6A": float(w.mean()),
                "mean_ethanol_O_within_6A": float(e.mean()),
                "fraction_frames_any_ethanol_O_within_6A": float(np.mean(e > 0)),
                "local_ethanol_O_fraction_within_6A": local_fraction,
            })

            top = sorted(
                zip(unique_residues, occupancy),
                key=lambda item: item[1], reverse=True,
            )[:10]
            print(
                f"{label} {period_name} ns {state_name}: "
                f"{nframes} frames; mean water O={w.mean():.2f}; "
                f"ethanol O={e.mean():.2f}; local ethanol fraction="
                f"{local_fraction:.3f}", flush=True,
            )
            print(
                "  Top contacts: " + ", ".join(
                    f"{key[3]}{key[1]} {frequency:.1%}"
                    for key, frequency in top if frequency > 0
                ), flush=True,
            )

    del tr, contacts, cms, ct

contact_path = analysis_dir / "ZQG_residue_contacts_by_state.csv"
with contact_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(contact_summary[0]))
    writer.writeheader()
    writer.writerows(contact_summary)

solvent_path = analysis_dir / "Cys64_solvent_by_ZQG_state.csv"
with solvent_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(solvent_summary[0]))
    writer.writeheader()
    writer.writerows(solvent_summary)

print(f"Saved: {contact_path}")
print(f"Saved: {solvent_path}")
