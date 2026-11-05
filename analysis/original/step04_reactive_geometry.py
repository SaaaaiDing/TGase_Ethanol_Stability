#!/usr/bin/env python3
"""Measure Cys64 SG to ZQG Gln side-chain carbonyl geometry in Desmond MD."""

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
        (541, 2679, 2680),
    ),
    "ethanol20": (
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/1IU4_EtOH20_md_55774175-out.cms",
        root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/production_trj",
        (541, 5168, 5169),
    ),
    "ethanol50": (
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/1IU4_EtOH50_md_55831542-out.cms",
        root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/production_trj",
        (541, 5168, 5169),
    ),
}


def minimum_image(displacement, lengths):
    return displacement - np.rint(displacement / lengths) * lengths


for label, (cms_path, trj_path, aids) in systems.items():
    if not cms_path.is_file() or not trj_path.exists():
        raise FileNotFoundError(f"Missing input for {label}: {cms_path} or {trj_path}")

    _, cms = topo.read_cms(str(cms_path))
    sg, carbon, oxygen = (cms.fsys_ct.atom[aid] for aid in aids)
    assert (sg.pdbres.strip(), sg.resnum, sg.pdbname.strip()) == ("CYS", 64, "SG")
    assert (carbon.pdbres.strip(), carbon.pdbname.strip(), carbon.element) == ("LIG1", "C28", "C")
    assert (oxygen.pdbres.strip(), oxygen.pdbname.strip(), oxygen.element) == ("LIG1", "O22", "O")
    assert oxygen.index in {neighbor.index for neighbor in carbon.bonded_atoms}
    amide_nitrogens = [
        neighbor for neighbor in carbon.bonded_atoms
        if neighbor.element == "N"
        and sum(n.element == "H" for n in neighbor.bonded_atoms) == 2
    ]
    assert len(amide_nitrogens) == 1, "C28 is not the expected primary amide carbonyl"

    gids = topo.aids2gids(cms, list(aids), include_pseudoatoms=False)
    assert len(gids) == 3
    output = analysis_dir / f"{label}_reactive_geometry.csv"
    rows = []
    frames = traj.read_traj(str(trj_path))
    for frame in frames:
        box = np.asarray(frame.box, dtype=float)
        if box.shape != (3, 3) or not np.allclose(box, np.diag(np.diag(box)), atol=1e-3):
            raise RuntimeError(f"{label}: expected an orthorhombic frame box")
        lengths = np.diag(box)
        if np.any(lengths <= 0):
            raise RuntimeError(f"{label}: invalid box lengths")

        xyz_sg, xyz_c, xyz_o = np.asarray(frame.pos(gids), dtype=float)
        c_to_sg = minimum_image(xyz_sg - xyz_c, lengths)
        c_to_o = minimum_image(xyz_o - xyz_c, lengths)
        distance = float(np.linalg.norm(c_to_sg))
        oxygen_length = float(np.linalg.norm(c_to_o))
        if not 0.9 < oxygen_length < 1.8 or distance == 0:
            raise RuntimeError(f"{label}: unexpected carbonyl bond or SG position")
        cosine = np.dot(c_to_o, c_to_sg) / (oxygen_length * distance)
        angle = float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
        rows.append((frame.time / 1000.0, distance, angle))

    if len(rows) != 1002:
        raise RuntimeError(f"{label}: expected 1002 frames, found {len(rows)}")
    times = np.array([row[0] for row in rows])
    if not np.all(np.diff(times) > 0) or not 199.9 < times[-1] < 200.1:
        raise RuntimeError(f"{label}: unexpected trajectory times")

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_ns", "Cys64_SG_to_ZQG_C28_A", "O22_C28_SG_angle_deg"])
        writer.writerows(rows)

    print(f"\n{label}: {len(rows)} frames; SG/C28/O22 AIDs={aids}")
    for left, right in ((0, 50), (50, 100), (100, 150), (150, 200.01)):
        selected = [row for row in rows if left <= row[0] < right]
        distances = np.array([row[1] for row in selected])
        angles = np.array([row[2] for row in selected])
        close = distances <= 5.0
        angle_text = f"{np.median(angles[close]):.1f}" if np.any(close) else "NA"
        print(
            f"  {left:3d}-{min(right, 200):3.0f} ns: "
            f"median SG-C28={np.median(distances):5.2f} A; "
            f"<=4 A={np.mean(distances <= 4.0):5.1%}; "
            f"<=5 A={np.mean(close):5.1%}; "
            f"median O22-C28-SG angle when <=5 A={angle_text} deg"
        )
    print(f"Saved: {output}")
