#!/usr/bin/env python3
"""Analyze existing MD medoids and draw compact atom-level contact schematics.

Run under Schrödinger 2023.4 from the MD analysis directory:
    $SCHRODINGER/run python3 step13_representative_atom_contacts.py

The five snapshots must have been exported by step09_rmsf_representatives.py.
"""

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from schrodinger import structure


CASES = (
    ("water_late_near_hbond_0580", "Water", "Late near + carbonyl H-bond",
     116.0, [64, 65, 75, 62, 254, 277, 278], "#0072B2"),
    ("ethanol20_late_far_0627", "20% ethanol", "Late displaced",
     125.4, [208, 238, 237, 239, 75, 253], "#D55E00"),
    ("ethanol50_late_near_0709", "50% ethanol", "Late near, alternative pose",
     141.8, [64, 65, 252, 253, 254, 276, 277, 278], "#009E73"),
    ("water_late_far_0746", "Water", "Late displaced",
     149.2, [75, 69, 62, 74, 77, 64], "#0072B2"),
    ("ethanol20_early_near_0129", "20% ethanol", "Early near",
     25.8, [64, 65, 252, 253, 254, 277, 278], "#D55E00"),
)
MAIN = (0, 1, 2)
SUPPLEMENT = (3, 4)
CONTACT_CUTOFF = 4.5
HBOND_DA_CUTOFF = 3.5
HBOND_ANGLE_CUTOFF = 120.0


def atom_id(atom):
    return f"{atom.pdbres.strip()}{atom.resnum} {atom.pdbname.strip()}"


def norm(a, b):
    return float(np.linalg.norm(np.asarray(a.xyz, dtype=float) -
                                np.asarray(b.xyz, dtype=float)))


def dha_angle(donor, hydrogen, acceptor):
    d = np.asarray(donor.xyz, dtype=float)
    h = np.asarray(hydrogen.xyz, dtype=float)
    a = np.asarray(acceptor.xyz, dtype=float)
    u, v = d - h, a - h
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v) /
        (np.linalg.norm(u) * np.linalg.norm(v)), -1, 1))))


def donor_hydrogens(atom):
    if atom.element not in ("N", "O"):
        return []
    return [neighbor for neighbor in atom.bonded_atoms if neighbor.element == "H"]


def geometric_hbond(protein_atom, ligand_atom):
    """Return best D-H...O angle; nitrogen acceptor assignments are excluded."""
    trials = []
    if ligand_atom.element == "O":
        trials.extend((protein_atom, h, ligand_atom, "protein donor")
                      for h in donor_hydrogens(protein_atom))
    if protein_atom.element == "O":
        trials.extend((ligand_atom, h, protein_atom, "ligand donor")
                      for h in donor_hydrogens(ligand_atom))
    valid = []
    for donor, hydrogen, acceptor, direction in trials:
        distance = norm(donor, acceptor)
        angle = dha_angle(donor, hydrogen, acceptor)
        if distance <= HBOND_DA_CUTOFF and angle >= HBOND_ANGLE_CUTOFF:
            valid.append((angle, direction))
    return max(valid, default=None)


def unique_atom(atoms, residue, name):
    matches = [a for a in atoms if a.resnum == residue and
               a.pdbname.strip() == name]
    if len(matches) != 1:
        raise ValueError(f"Expected one {residue} {name}, found {len(matches)}")
    return matches[0]


def analyze_case(directory, item):
    tag, system, state, time, focus, color = item
    path = directory / f"{tag}.maegz"
    if not path.is_file():
        raise FileNotFoundError(path)
    st = structure.StructureReader.read(str(path))
    ligand = [atom for atom in st.atom if atom.pdbres.strip() == "LIG1"]
    protein = [atom for atom in st.atom if atom.pdbres.strip() != "LIG1"
               and atom.chain.strip() == "A"]
    if len(protein) != 5166 or len(ligand) != 42 or st.atom_total != 5208:
        raise ValueError(f"{tag}: expected 5166 protein and 42 ligand atoms")
    protein_heavy = [a for a in protein if a.element != "H"]
    ligand_heavy = [a for a in ligand if a.element != "H"]
    if len(ligand_heavy) != 24:
        raise ValueError(f"{tag}: expected 24 ligand heavy atoms")
    cys_sg = unique_atom(protein, 64, "SG")
    cys_n = unique_atom(protein, 64, "N")
    zqg_c28 = unique_atom(ligand, ligand[0].resnum, "C28")
    zqg_o22 = unique_atom(ligand, ligand[0].resnum, "O22")
    sg_c28 = norm(cys_sg, zqg_c28)
    cys_h = donor_hydrogens(cys_n)
    if len(cys_h) != 1:
        raise ValueError(f"{tag}: Cys64 backbone N needs one explicit H")
    nh_o22_A = norm(cys_n, zqg_o22)
    nh_o22_deg = dha_angle(cys_n, cys_h[0], zqg_o22)
    oxyanion_proxy = nh_o22_A <= HBOND_DA_CUTOFF and nh_o22_deg >= HBOND_ANGLE_CUTOFF
    expected_near = "late_near" in tag or "early_near" in tag
    if expected_near != (sg_c28 <= 5.0):
        raise ValueError(f"{tag}: reactive state disagrees with snapshot geometry")
    if "late_near_hbond" in tag and not oxyanion_proxy:
        raise ValueError(f"{tag}: expected Cys64 NH--O22 geometry")

    heavy_xyz = np.asarray([a.xyz for a in protein_heavy], dtype=float)
    ligand_xyz = np.asarray([a.xyz for a in ligand_heavy], dtype=float)
    distances = np.linalg.norm(heavy_xyz[:, None, :] - ligand_xyz[None, :, :], axis=2)
    all_contacts = []
    for i, j in np.argwhere(distances <= CONTACT_CUTOFF):
        p, l = protein_heavy[i], ligand_heavy[j]
        hb = geometric_hbond(p, l)
        all_contacts.append({
            "case": tag, "system": system, "state": state,
            "time_ns": time, "protein_residue": f"{p.pdbres.strip()}{p.resnum}",
            "protein_atom": p.pdbname.strip(), "ligand_atom": l.pdbname.strip(),
            "distance_A": float(distances[i, j]),
            "geometric_DH_O": int(hb is not None),
            "D_H_O_angle_deg": hb[0] if hb is not None else "",
            "donor_side": hb[1] if hb is not None else "",
        })

    chosen = []
    for residue_number in focus:
        residue_contacts = [r for r in all_contacts if
                            int(r["protein_residue"][3:]) == residue_number]
        if not residue_contacts:
            continue
        hbond_contacts = [r for r in residue_contacts if r["geometric_DH_O"]]
        polar_contacts = [r for r in residue_contacts if
                          r["protein_atom"].startswith(("N", "O")) and
                          r["ligand_atom"].startswith(("N", "O")) and
                          r["distance_A"] <= 3.5]
        candidates = hbond_contacts or polar_contacts or residue_contacts
        chosen.append(min(candidates, key=lambda r: r["distance_A"]))
    if not chosen:
        raise ValueError(f"{tag}: no focused contacts within {CONTACT_CUTOFF} A")
    return {
        "tag": tag, "system": system, "state": state, "time_ns": time,
        "color": color, "SG_C28_A": sg_c28,
        "Cys64_N_O22_A": nh_o22_A,
        "Cys64_N_H_O22_angle_deg": nh_o22_deg,
        "Cys64_NH_O22_proxy": oxyanion_proxy,
        "contacts": chosen,
        "all_contacts": all_contacts,
    }


def draw_panel(ax, result, label):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    color = result["color"]
    ax.add_patch(plt.Rectangle((0.05, 9.3), 9.9, .58,
                               facecolor=color, edgecolor="none"))
    ax.text(.26, 9.58, f"{label}  {result['system']} | {result['state']}",
            color="white", weight="bold", size=11, va="center")
    ax.text(.27, 8.92, f"Representative frame: {result['time_ns']:.1f} ns",
            size=9, color="#333333")
    ax.text(.27, 8.42,
            f"Cys64 SG–C28: {result['SG_C28_A']:.1f} Å   |   "
            f"Cys64 NH–O22 geometry: {'yes' if result['Cys64_NH_O22_proxy'] else 'no'}",
            size=8.6, color="#333333")
    ax.text(.27, 7.86, "PROTEIN ATOM", size=7.5, color="#666666", weight="bold")
    ax.text(7.25, 7.86, "ZQG ATOM", size=7.5, color="#666666", weight="bold")
    contacts = result["contacts"][:8]
    for idx, row in enumerate(contacts):
        y = 7.25 - idx * .78
        left = f"{row['protein_residue']} {row['protein_atom']}"
        right = row["ligand_atom"]
        hb = bool(row["geometric_DH_O"])
        line_color = "#0072B2" if hb else "#7D8991"
        ax.text(.27, y, left, ha="left", va="center", size=9.1)
        ax.plot([4.15, 7.1], [y, y], ls="--", lw=1.4 if hb else 1.0,
                color=line_color)
        ax.text(5.62, y + .17, f"{row['distance_A']:.1f} Å",
                ha="center", va="bottom", size=8, color="#333333",
                bbox=dict(facecolor="white", edgecolor="none", pad=.2))
        ax.text(7.25, y, right, ha="left", va="center", size=9.1,
                weight="bold" if hb else "normal")
        if hb:
            ax.text(9.65, y, "H", ha="right", va="center", color=line_color,
                    weight="bold", size=8)
    ax.text(.27, .22, "H: geometric donor–H···O contact", size=7.5,
            color="#555555")


def write_figure(results, output, selected, basename):
    n = len(selected)
    fig, axes = plt.subplots(1, n, figsize=(n * 5.05, 6.1),
                             layout="constrained")
    if n == 1:
        axes = [axes]
    for label, index, ax in zip("ABCDE", selected, axes):
        draw_panel(ax, results[index], label)
    fig.suptitle("1IU4–ZQG: representative atom contacts",
                 fontsize=15, weight="bold")
    for ext in ("svg", "png"):
        fig.savefig(output / f"{basename}.{ext}", dpi=400,
                    bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", type=Path,
                        default=Path.cwd())
    args = parser.parse_args()
    directory = args.analysis_dir.expanduser().resolve()
    source = directory / "representatives"
    output = directory / "representative_contacts"
    output.mkdir(exist_ok=True)
    results = [analyze_case(source, item) for item in CASES]
    rows = [row for result in results for row in result["all_contacts"]]
    with (output / "all_atom_contacts_4p5A.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (output / "representative_contact_summary.json").open("w") as handle:
        json.dump([{k: v for k, v in result.items() if k != "all_contacts"}
                   for result in results], handle, indent=2)
    write_figure(results, output, MAIN, "three_condition_atom_contacts")
    write_figure(results, output, SUPPLEMENT, "state_transition_atom_contacts")
    for result in results:
        print(f"{result['tag']}: SG-C28 {result['SG_C28_A']:.2f} A; "
              f"NH-O22 {result['Cys64_N_O22_A']:.2f} A; "
              f"NH-O22 proxy {result['Cys64_NH_O22_proxy']}; "
              f"focus contacts {len(result['contacts'])}", flush=True)
    print(f"Saved contact figures and tables: {output}", flush=True)


if __name__ == "__main__":
    main()
