#!/usr/bin/env python3
"""Measure substrate carbonyl O22 contacts with Cys64 and Val65 backbone NH."""

import csv
import os
from pathlib import Path

import numpy as np
from schrodinger.application.desmond.packages import topo, traj


analysis_dir = Path(os.environ["ANALYSIS_DIR"]).expanduser().resolve()
root = analysis_dir.parent
systems = {
    "water": (
        root / "1IU4_0/desmond_md_job_4/1IU4_md_batch-out.cms",
        root / "1IU4_0/desmond_md_job_4/1IU4_md_batch_trj",
        2679, 2680,
    ),
    "ethanol20": (
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/1IU4_EtOH20_md_55774175-out.cms",
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/production_trj",
        5168, 5169,
    ),
    "ethanol50": (
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/1IU4_EtOH50_md_55831542-out.cms",
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/production_trj",
        5168, 5169,
    ),
}


def displacement(first, second, lengths):
    delta = first - second
    return delta - np.rint(delta / lengths) * lengths


def donor_geometry(n_xyz, h_xyz, o_xyz, lengths):
    n_to_o = displacement(o_xyz, n_xyz, lengths)
    h_to_n = displacement(n_xyz, h_xyz, lengths)
    h_to_o = displacement(o_xyz, h_xyz, lengths)
    distance = float(np.linalg.norm(n_to_o))
    cosine = np.dot(h_to_n, h_to_o) / (
        np.linalg.norm(h_to_n) * np.linalg.norm(h_to_o)
    )
    angle = float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
    proxy = distance <= 3.5 and angle >= 120.0
    return distance, angle, int(proxy)


for label, (cms_path, trj_path, carbon_aid, oxygen_aid) in systems.items():
    geometry_path = analysis_dir / f"{label}_reactive_geometry.csv"
    if not cms_path.is_file() or not trj_path.exists() or not geometry_path.is_file():
        raise FileNotFoundError(
            f"Missing input for {label}: {cms_path}, {trj_path}, or {geometry_path}"
        )

    with geometry_path.open(newline="") as handle:
        reference = list(csv.DictReader(handle))
    if len(reference) != 1002:
        raise RuntimeError(f"{label}: expected 1002 reference geometry frames")

    _, cms = topo.read_cms(str(cms_path))
    ct = cms.fsys_ct
    sg, carbon, oxygen = (ct.atom[aid] for aid in (541, carbon_aid, oxygen_aid))
    if (sg.pdbres.strip(), sg.resnum, sg.pdbname.strip()) != ("CYS", 64, "SG"):
        raise RuntimeError(f"{label}: Cys64 SG selection failed")
    if (carbon.pdbres.strip(), carbon.pdbname.strip()) != ("LIG1", "C28"):
        raise RuntimeError(f"{label}: reactive carbon selection failed")
    if (oxygen.pdbres.strip(), oxygen.pdbname.strip()) != ("LIG1", "O22"):
        raise RuntimeError(f"{label}: carbonyl oxygen selection failed")

    donor_aids = []
    for residue_name, number in (("CYS", 64), ("VAL", 65)):
        nitrogens = [
            atom for atom in ct.atom
            if atom.chain.strip() == "A"
            and atom.pdbres.strip() == residue_name
            and atom.resnum == number
            and atom.pdbname.strip() == "N"
        ]
        if len(nitrogens) != 1:
            raise RuntimeError(f"{label}: expected one {residue_name}{number} backbone N")
        nitrogen = nitrogens[0]
        hydrogens = [atom for atom in nitrogen.bonded_atoms if atom.element == "H"]
        if len(hydrogens) != 1:
            raise RuntimeError(
                f"{label}: expected one hydrogen bonded to "
                f"{residue_name}{number} backbone N, found {len(hydrogens)}"
            )
        donor_aids.extend((nitrogen.index, hydrogens[0].index))
        print(
            f"{label}: {residue_name}{number} backbone N/H "
            f"AIDs={nitrogen.index}/{hydrogens[0].index}", flush=True,
        )

    aids = [sg.index, carbon.index, oxygen.index] + donor_aids
    gids = topo.aids2gids(cms, aids, include_pseudoatoms=False)
    if len(gids) != 7:
        raise RuntimeError(f"{label}: AID to GID mapping failed")

    frames = traj.read_traj(str(trj_path))
    if len(frames) != len(reference):
        raise RuntimeError(f"{label}: trajectory and reference frame counts differ")
    rows = []
    for index, frame in enumerate(frames):
        time_ns = frame.time / 1000.0
        if abs(time_ns - float(reference[index]["time_ns"])) > 1e-5:
            raise RuntimeError(f"{label}: frame {index} time does not match reference")
        box = np.asarray(frame.box, dtype=float)
        if box.shape != (3, 3) or not np.allclose(box, np.diag(np.diag(box)), atol=1e-3):
            raise RuntimeError(f"{label}: orthorhombic box required")
        lengths = np.diag(box)
        if np.any(lengths <= 0):
            raise RuntimeError(f"{label}: invalid box lengths")

        sg_xyz, c_xyz, o_xyz, cys_n, cys_h, val_n, val_h = np.asarray(
            frame.pos(gids), dtype=float
        )
        sg_c28 = float(np.linalg.norm(displacement(sg_xyz, c_xyz, lengths)))
        ref_distance = float(reference[index]["Cys64_SG_to_ZQG_C28_A"])
        if abs(sg_c28 - ref_distance) > 0.02:
            raise RuntimeError(f"{label}: frame {index} SG-C28 distance disagrees with step04")

        cys_distance, cys_angle, cys_proxy = donor_geometry(
            cys_n, cys_h, o_xyz, lengths
        )
        val_distance, val_angle, val_proxy = donor_geometry(
            val_n, val_h, o_xyz, lengths
        )
        rows.append({
            "time_ns": time_ns,
            "SG_C28_A": sg_c28,
            "Cys64_N_O22_A": cys_distance,
            "Cys64_N_H_O22_deg": cys_angle,
            "Cys64_NH_O22_proxy": cys_proxy,
            "Val65_N_O22_A": val_distance,
            "Val65_N_H_O22_deg": val_angle,
            "Val65_NH_O22_proxy": val_proxy,
            "either_backbone_NH_O22_proxy": int(bool(cys_proxy or val_proxy)),
        })

    output = analysis_dir / f"{label}_carbonyl_backbone_hbonds.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {output}", flush=True)

    for start, stop in ((0.0, 50.0), (50.0, 200.01)):
        for state_name, predicate in (
            ("near", lambda d: d <= 5.0),
            ("far", lambda d: d > 5.0),
        ):
            group = [
                row for row in rows
                if start <= row["time_ns"] < stop
                and predicate(row["SG_C28_A"])
            ]
            if not group:
                continue
            fractions = [
                np.mean([row[key] for row in group])
                for key in (
                    "Cys64_NH_O22_proxy",
                    "Val65_NH_O22_proxy",
                    "either_backbone_NH_O22_proxy",
                )
            ]
            print(
                f"{label} {start:.0f}-{min(stop, 200):.0f} ns "
                f"{state_name}: n={len(group)}; "
                f"Cys64 NH-O22={fractions[0]:.1%}; "
                f"Val65 NH-O22={fractions[1]:.1%}; "
                f"either={fractions[2]:.1%}", flush=True,
            )
