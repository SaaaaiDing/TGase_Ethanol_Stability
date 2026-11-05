#!/usr/bin/env python3
"""Compare aligned protein C-alpha RMSF and export state-specific MD snapshots."""

import csv
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from schrodinger.application.desmond.packages import topo, traj


analysis_dir = Path(os.environ["ANALYSIS_DIR"]).expanduser().resolve()
root = analysis_dir.parent
systems = {
    "water": (
        root / "1IU4_0/desmond_md_job_4/1IU4_md_batch-out.cms",
        root / "1IU4_0/desmond_md_job_4/1IU4_md_batch_trj",
        2679,
    ),
    "ethanol20": (
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/1IU4_EtOH20_md_55774175-out.cms",
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/production_trj",
        5168,
    ),
    "ethanol50": (
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/1IU4_EtOH50_md_55831542-out.cms",
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/production_trj",
        5168,
    ),
}
colors = {"water": "#333333", "ethanol20": "#dc7624", "ethanol50": "#2678a4"}
names = {"water": "Water", "ethanol20": "20% ethanol", "ethanol50": "50% ethanol"}


def residue_key(atom):
    return atom.chain.strip(), atom.resnum, atom.inscode.strip()


def minimum_image(delta, lengths):
    return delta - np.rint(delta / lengths) * lengths


def unwrap_ca(xyz, lengths):
    steps = minimum_image(np.diff(xyz, axis=0), lengths)
    if np.any(np.linalg.norm(steps, axis=1) > 8.0):
        raise RuntimeError("Noncontiguous C-alpha trace; inspect chain or residue mapping")
    return np.vstack((xyz[0], xyz[0] + np.cumsum(steps, axis=0)))


def fit(mobile, reference):
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    left, _, right = np.linalg.svd(
        (mobile - mobile_center).T @ (reference - reference_center)
    )
    sign = np.linalg.det(left @ right)
    rotation = left @ np.diag([1.0, 1.0, np.sign(sign)]) @ right
    translation = reference_center - mobile_center @ rotation
    return rotation, translation


def write_representative(cms, frame, ca_aids, protein_aids, ligand_aids,
                         carbon_aid, output):
    selected_aids = sorted(set(protein_aids) | set(ligand_aids))
    all_gids = topo.aids2gids(cms, selected_aids, include_pseudoatoms=False)
    if len(all_gids) != len(selected_aids):
        raise RuntimeError("Representative AID/GID mapping failed")
    raw = np.asarray(frame.pos(all_gids), dtype=float)
    indices = {aid: i for i, aid in enumerate(selected_aids)}
    box = np.asarray(frame.box, dtype=float)
    lengths = np.diag(box)
    if box.shape != (3, 3) or not np.allclose(box, np.diag(lengths), atol=1e-3):
        raise RuntimeError("Orthorhombic trajectory required")
    ca_raw = raw[[indices[aid] for aid in ca_aids]]
    ca_unwrapped = unwrap_ca(ca_raw, lengths)
    ca_by_residue = {
        residue_key(cms.fsys_ct.atom[aid]): i for i, aid in enumerate(ca_aids)
    }
    positions = np.empty_like(raw)
    for aid in protein_aids:
        i = indices[aid]
        ca_i = ca_by_residue[residue_key(cms.fsys_ct.atom[aid])]
        positions[i] = ca_unwrapped[ca_i] + minimum_image(
            raw[i] - ca_raw[ca_i], lengths
        )
    sg_i = indices[541]
    carbon_i = indices[carbon_aid]
    carbon_pos = positions[sg_i] + minimum_image(
        raw[carbon_i] - raw[sg_i], lengths
    )
    for aid in ligand_aids:
        i = indices[aid]
        positions[i] = carbon_pos + minimum_image(
            raw[i] - raw[carbon_i], lengths
        )
    st = cms.extract(selected_aids)
    if st.atom_total != len(selected_aids):
        raise RuntimeError("Extracted snapshot has an unexpected atom count")
    st.setXYZ(positions)
    st.write(str(output))


fig, ax = plt.subplots(figsize=(11, 4.6), constrained_layout=True)
all_rmsf = {}
reference_keys = None
representative_records = []
representative_dir = analysis_dir / "representatives"
representative_dir.mkdir(exist_ok=True)

for label, (cms_path, trj_path, carbon_aid) in systems.items():
    state_path = analysis_dir / f"{label}_carbonyl_backbone_hbonds.csv"
    for path in (cms_path, trj_path, state_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing input: {path}")
    with state_path.open(newline="") as handle:
        state_rows = list(csv.DictReader(handle))
    if len(state_rows) != 1002:
        raise RuntimeError(f"{label}: expected 1002 state records")

    _, cms = topo.read_cms(str(cms_path))
    ct = cms.fsys_ct
    protein_aids = sorted(cms.select_atom("protein"))
    protein_aid_set = set(protein_aids)
    ca_aids = sorted(
        (aid for aid in protein_aids if ct.atom[aid].pdbname.strip() == "CA"),
        key=lambda aid: residue_key(ct.atom[aid]),
    )
    keys = [residue_key(ct.atom[aid]) for aid in ca_aids]
    if len(keys) != 331 or len(set(keys)) != len(keys):
        raise RuntimeError(f"{label}: invalid protein C-alpha mapping")
    if reference_keys is None:
        reference_keys = keys
    elif keys != reference_keys:
        raise RuntimeError(f"{label}: protein residues differ from the water system")
    ligand_aids = [
        atom.index for atom in ct.atom if atom.pdbres.strip() == "LIG1"
    ]
    ligand_heavy = [aid for aid in ligand_aids if ct.atom[aid].element != "H"]
    if len(protein_aids) != 5166 or len(ligand_heavy) != 24:
        raise RuntimeError(f"{label}: protein or ZQG atom selection failed")
    if carbon_aid not in ligand_heavy or 541 not in protein_aid_set:
        raise RuntimeError(f"{label}: reactive atom selection failed")
    ca_gids = topo.aids2gids(cms, ca_aids, include_pseudoatoms=False)
    heavy_gids = topo.aids2gids(cms, ligand_heavy, include_pseudoatoms=False)
    sg_gid, carbon_gid = topo.aids2gids(
        cms, [541, carbon_aid], include_pseudoatoms=False
    )
    frames = traj.read_traj(str(trj_path))
    if len(frames) != len(state_rows):
        raise RuntimeError(f"{label}: trajectory/state frame count mismatch")

    times = np.empty(len(frames))
    protein_positions = np.empty((len(frames), len(ca_aids), 3))
    ligand_positions = np.empty((len(frames), len(ligand_heavy), 3))
    sg_ca_i = keys.index(("A", 64, ""))
    for i, frame in enumerate(frames):
        times[i] = frame.time / 1000.0
        if abs(times[i] - float(state_rows[i]["time_ns"])) > 1e-5:
            raise RuntimeError(f"{label}: time mismatch at frame {i}")
        box = np.asarray(frame.box, dtype=float)
        if box.shape != (3, 3) or not np.allclose(box, np.diag(np.diag(box)), atol=1e-3):
            raise RuntimeError(f"{label}: expected orthorhombic box")
        lengths = np.diag(box)
        if np.any(lengths <= 0):
            raise RuntimeError(f"{label}: invalid box dimensions")
        ca_raw = np.asarray(frame.pos(ca_gids), dtype=float)
        ca_unwrapped = unwrap_ca(ca_raw, lengths)
        sg_raw = np.asarray(frame.pos(sg_gid), dtype=float).reshape(3)
        carbon_raw = np.asarray(frame.pos(carbon_gid), dtype=float).reshape(3)
        sg_unwrapped = ca_unwrapped[sg_ca_i] + minimum_image(
            sg_raw - ca_raw[sg_ca_i], lengths
        )
        carbon_unwrapped = sg_unwrapped + minimum_image(
            carbon_raw - sg_raw, lengths
        )
        ligand_raw = np.asarray(frame.pos(heavy_gids), dtype=float)
        protein_positions[i] = ca_unwrapped
        ligand_positions[i] = carbon_unwrapped + minimum_image(
            ligand_raw - carbon_raw, lengths
        )
        distance = np.linalg.norm(carbon_unwrapped - sg_unwrapped)
        if abs(distance - float(state_rows[i]["SG_C28_A"])) > 0.02:
            raise RuntimeError(f"{label}: frame {i} disagrees with Step 08")

    late = (times >= 50.0) & (times <= 200.01)
    if late.sum() != 751:
        raise RuntimeError(f"{label}: expected 751 frames from 50 to 200 ns")
    ref = protein_positions[np.flatnonzero(late)[0]].copy()
    for _ in range(3):
        aligned = []
        for p in protein_positions[late]:
            rotation, translation = fit(p, ref)
            aligned.append(p @ rotation + translation)
        aligned = np.asarray(aligned)
        ref = aligned.mean(axis=0)
    aligned_ca = np.empty_like(protein_positions)
    aligned_ligand = np.empty_like(ligand_positions)
    for i, p in enumerate(protein_positions):
        rotation, translation = fit(p, ref)
        aligned_ca[i] = p @ rotation + translation
        aligned_ligand[i] = ligand_positions[i] @ rotation + translation
    mean_ca = aligned_ca[late].mean(axis=0)
    rmsf = np.sqrt(np.mean(np.sum((aligned_ca[late] - mean_ca) ** 2, axis=2), axis=0))
    all_rmsf[label] = rmsf
    ax.plot([key[1] for key in keys], rmsf, color=colors[label], label=names[label], lw=1.1)
    print(f"{label}: 50-200 ns C-alpha median RMSF = {np.median(rmsf):.2f} A", flush=True)

    near = np.array([float(row["SG_C28_A"]) <= 5.0 for row in state_rows])
    hb = np.array([int(row["either_backbone_NH_O22_proxy"]) == 1 for row in state_rows])
    groups = {
        "early_near": (times < 50.0) & near,
        "late_near": late & near,
        "late_far": late & ~near,
        "late_near_hbond": late & near & hb,
    }
    for group_name, mask in groups.items():
        candidates = np.flatnonzero(mask)
        if len(candidates) < 20:
            continue
        group_ca = aligned_ca[candidates]
        group_ligand = aligned_ligand[candidates]
        ca_center = group_ca.mean(axis=0)
        ligand_center = group_ligand.mean(axis=0)
        ca_distance2 = np.mean(np.sum((group_ca - ca_center) ** 2, axis=2), axis=1)
        ligand_distance2 = np.mean(
            np.sum((group_ligand - ligand_center) ** 2, axis=2), axis=1
        )
        scores = ligand_distance2 + 0.25 * ca_distance2
        frame_i = int(candidates[np.argmin(scores)])
        output = representative_dir / f"{label}_{group_name}_{frame_i:04d}.maegz"
        write_representative(
            cms, frames[frame_i], ca_aids, protein_aids, ligand_aids,
            carbon_aid, output,
        )
        representative_records.append({
            "system": label,
            "state": group_name,
            "frames_in_state": len(candidates),
            "frame_index_zero_based": frame_i,
            "time_ns": times[frame_i],
            "SG_C28_A": state_rows[frame_i]["SG_C28_A"],
            "carbonyl_backbone_hbond_proxy": int(hb[frame_i]),
            "snapshot_maegz": str(output),
        })
        print(
            f"{label} {group_name}: n={len(candidates)}, "
            f"representative frame={frame_i}, t={times[frame_i]:.2f} ns; {output}",
            flush=True,
        )

output_csv = analysis_dir / "protein_CA_RMSF_50_200ns.csv"
with output_csv.open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["chain", "resnum", "insertion_code", "water_A", "ethanol20_A", "ethanol50_A"])
    for i, (chain, number, insertion) in enumerate(reference_keys):
        writer.writerow([chain, number, insertion] + [
            float(all_rmsf[label][i]) for label in systems
        ])

metadata_csv = representative_dir / "representative_frames.csv"
with metadata_csv.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(representative_records[0]))
    writer.writeheader()
    writer.writerows(representative_records)

ax.set_xlabel("Protein residue number (chain A)")
ax.set_ylabel("C-alpha RMSF (A), 50-200 ns")
ax.set_xlim(min(key[1] for key in reference_keys), max(key[1] for key in reference_keys))
ax.legend(frameon=False, ncol=3)
ax.spines[["top", "right"]].set_visible(False)
fig.savefig(analysis_dir / "protein_CA_RMSF_50_200ns.svg")
fig.savefig(analysis_dir / "protein_CA_RMSF_50_200ns.png", dpi=300)
print(f"Saved: {output_csv}")
print(f"Saved: {metadata_csv}")
