import csv
from pathlib import Path

import numpy as np
from schrodinger.application.desmond.packages import topo, traj

ROOT = Path(
    "/storage/home/sqd5856/Desktop/default/HanLiu/prep_files/yaml/1IU4/"
    "boltz_results_1IU4/predictions/1IU4"
)
OUTPUT = ROOT / "md_analysis_0_20_50"

SYSTEMS = {
    "water": (
        ROOT / "1IU4_0/desmond_md_job_4",
        "1IU4_md_batch-out.cms",
        "1IU4_md_batch_trj",
    ),
    "ethanol20": (
        ROOT / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175",
        "1IU4_EtOH20_md_55774175-out.cms",
        "production_trj",
    ),
    "ethanol50": (
        ROOT / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/"
        "npt_500ps/production_55831542",
        "1IU4_EtOH50_md_55831542-out.cms",
        "production_trj",
    ),
}

for label, (directory, cms_name, trj_name) in SYSTEMS.items():
    print(f"\nAnalyzing {label}...", flush=True)

    _, cms_model = topo.read_cms(str(directory / cms_name))
    system = cms_model.fsys_ct
    protein_aids = set(cms_model.select_atom("protein"))

    cys_sg_aids = [
        atom.index for atom in system.atom
        if atom.index in protein_aids
        and atom.resnum == 64
        and atom.pdbres.strip() == "CYS"
        and atom.pdbname.strip() == "SG"
    ]
    ligand_aids = [
        atom.index for atom in system.atom
        if atom.pdbres.strip() == "LIG1" and atom.element != "H"
    ]

    if len(cys_sg_aids) != 1 or len(ligand_aids) != 24:
        raise RuntimeError(
            f"{label}: expected one Cys64 SG and 24 ZQG heavy atoms"
        )

    sg_gid = topo.aids2gids(
        cms_model, cys_sg_aids, include_pseudoatoms=False
    )[0]
    ligand_gids = topo.aids2gids(
        cms_model, ligand_aids, include_pseudoatoms=False
    )

    frames = traj.read_traj(str(directory / trj_name))
    rmsd_csv = OUTPUT / f"{label}_rmsd.csv"
    rmsd_data = np.loadtxt(rmsd_csv, delimiter=",", skiprows=1)

    if rmsd_data.shape != (len(frames), 3):
        raise RuntimeError(f"{label}: RMSD CSV and trajectory do not match")

    times_ns = []
    distances = []

    for frame in frames:
        box = np.asarray(frame.box, dtype=float)
        inverse_box = np.linalg.inv(box)

        sg_xyz = np.asarray(frame.pos(sg_gid), dtype=float).reshape(3)
        ligand_xyz = np.asarray(frame.pos(ligand_gids), dtype=float)

        displacement = ligand_xyz - sg_xyz
        fractional = displacement @ inverse_box
        minimum_image = (
            displacement - np.rint(fractional) @ box
        )
        distance = np.linalg.norm(minimum_image, axis=1).min()

        times_ns.append(frame.time / 1000.0)
        distances.append(float(distance))

    times_ns = np.asarray(times_ns)
    distances = np.asarray(distances)

    if not np.allclose(times_ns, rmsd_data[:, 0], atol=0.001):
        raise RuntimeError(f"{label}: time axes do not match")
    if not np.isfinite(distances).all():
        raise RuntimeError(f"{label}: invalid distances")

    output_csv = OUTPUT / f"{label}_active_site_proximity.csv"
    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "time_ns",
            "Cys64_SG_to_ZQG_nearest_heavy_A",
            "protein_CA_RMSD_A",
            "ZQG_heavy_RMSD_protein_aligned_A",
        ])
        for time, distance, row in zip(times_ns, distances, rmsd_data):
            writer.writerow([
                f"{time:.6f}",
                f"{distance:.6f}",
                f"{row[1]:.6f}",
                f"{row[2]:.6f}",
            ])

    for start, stop in ((0.0, 50.0), (50.0, 200.001)):
        use = (times_ns >= start) & (times_ns < stop)
        selected = distances[use]
        print(
            f"  {start:g}-{stop:g} ns: "
            f"median distance={np.median(selected):.2f} A; "
            f"within 5 A={np.mean(selected <= 5.0):.1%}; "
            f"within 8 A={np.mean(selected <= 8.0):.1%}"
        )

    print(f"  Saved: {output_csv}", flush=True)

print("\nACTIVE-SITE PROXIMITY ANALYSIS COMPLETED")
