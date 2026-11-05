import csv
from pathlib import Path

import numpy as np
from schrodinger.application.desmond.packages import analysis, topo, traj

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

OUTPUT.mkdir(parents=True, exist_ok=True)

for label, (directory, cms_name, trj_name) in SYSTEMS.items():
    print(f"Analyzing {label}...", flush=True)

    msys_model, cms_model = topo.read_cms(str(directory / cms_name))
    system = cms_model.fsys_ct

    protein_aids = cms_model.select_atom("protein")
    ca_aids = sorted(
        aid for aid in protein_aids
        if system.atom[aid].pdbname.strip() == "CA"
    )
    ligand_aids = [
        atom.index for atom in system.atom
        if atom.pdbres.strip() == "LIG1" and atom.element != "H"
    ]

    if len(ca_aids) != 331:
        raise RuntimeError(f"{label}: expected 331 protein CA atoms")
    if len(ligand_aids) != 24:
        raise RuntimeError(
            f"{label}: expected 24 ZQG heavy atoms, found {len(ligand_aids)}"
        )

    topo.make_glued_topology(msys_model, cms_model)

    frames = traj.read_traj(str(directory / trj_name))
    first = frames[0]

    ca_gids = topo.aids2gids(
        cms_model, ca_aids, include_pseudoatoms=False
    )
    ligand_gids = topo.aids2gids(
        cms_model, ligand_aids, include_pseudoatoms=False
    )

    ca_reference = first.pos(ca_gids).copy()
    ligand_reference = first.pos(ligand_gids).copy()

    protein_rmsd = analysis.RMSD(
        msys_model,
        cms_model,
        ca_aids,
        ca_reference,
        fit_aids=ca_aids,
        fit_ref_pos=ca_reference,
    )
    ligand_rmsd = analysis.LigandRMSD(
        msys_model,
        cms_model,
        ligand_aids,
        ligand_reference,
        fit_aids=ca_aids,
        fit_ref_pos=ca_reference,
    )

    protein_values, ligand_values = analysis.analyze(
        frames, protein_rmsd, ligand_rmsd
    )

    protein_values = np.asarray(protein_values, dtype=float).reshape(-1)
    ligand_values = np.asarray(ligand_values, dtype=float).reshape(-1)
    times_ns = np.asarray(
        [frame.time / 1000.0 for frame in frames], dtype=float
    )

    if not (
        len(times_ns) == len(protein_values) == len(ligand_values)
    ):
        raise RuntimeError(f"{label}: RMSD output length mismatch")

    if not (
        np.isfinite(protein_values).all()
        and np.isfinite(ligand_values).all()
    ):
        raise RuntimeError(f"{label}: non-finite RMSD values")

    output_csv = OUTPUT / f"{label}_rmsd.csv"
    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "time_ns",
            "protein_CA_RMSD_A",
            "ZQG_heavy_RMSD_protein_aligned_A",
        ])
        for time, protein, ligand in zip(
            times_ns, protein_values, ligand_values
        ):
            writer.writerow([
                f"{time:.6f}",
                f"{protein:.6f}",
                f"{ligand:.6f}",
            ])

    late = times_ns >= 50.0
    print(f"  Frames: {len(frames)}")
    print(f"  Protein CA atoms: {len(ca_aids)}")
    print(f"  ZQG heavy atoms: {len(ligand_aids)}")
    print(
        f"  Mean protein CA RMSD, 50-200 ns: "
        f"{protein_values[late].mean():.3f} A"
    )
    print(
        f"  Mean protein-aligned ZQG RMSD, 50-200 ns: "
        f"{ligand_values[late].mean():.3f} A"
    )
    print(f"  Saved: {output_csv}", flush=True)

print("RMSD ANALYSIS COMPLETED")
