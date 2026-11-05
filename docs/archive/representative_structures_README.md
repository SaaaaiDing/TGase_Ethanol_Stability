# 1IU4–ZQG: representative MD conformations

## Contents

- `1IU4_three_states_3D.png` and `.pdf`: actual aligned coordinates of three representative frames; local coordinate renderer shows the protein Cα trace, ligand and highlighted residues as atom-resolved sticks, plus measured interaction lines.
- `1IU4_three_states_2D_contacts.svg`, `.pdf`, `.png`: ZQG chemical connectivity and atom names taken from each Maestro structure, with specific protein–ligand hydrogen bonds. SVG/PDF are vector artwork.
- `water_116ns.pdb`, `ethanol20_125p4ns.pdb`, `ethanol50_141p8ns.pdb`: protein and ZQG coordinates, aligned on 331 protein Cα atoms. Solvent is omitted only from these figure structures.
- `*_pymol.pml`: PyMOL scripts for a full cartoon-and-sticks ray-traced image at 2400 × 1800 pixels. Run from this folder with `pymol -cq water_116ns_pymol.pml` and likewise for the other conditions. These scripts were prepared but could not be executed in the figure-creation environment, which lacks PyMOL.
- `render_1IU4_pymol.py`: preferred unified PyMOL script. Run `pymol -cq render_1IU4_pymol.py` from this folder. It reuses one overview camera and one close-up camera for all three prealigned states, validates the six displayed distances, writes three full-enzyme cartoons and three contact images, saves editable PyMOL sessions, and (if Pillow is available) assembles `pymol_figures/three_state_contact_panel.png`. This script was syntax checked and its input coordinates were independently verified, but PyMOL is absent here, so the rendering itself remains to be checked on the machine where it is run.
- `render_representatives.py`: reproducible parsing, measurements, alignment, and artwork generation; rerun beside an `upload/` folder holding the three input MAEGZ snapshots.

## Frame-level evidence

| System | Frame (ns) | Cys64 SG–ZQG C28 (Å) | Geometry-verified hydrogen bonds (heavy-atom distance, Å; D–H···A angle) |
|---|---:|---:|---|
| Water | 116.0 | 3.47 | Cys64 backbone N–H···O22 2.73, 167.7°; Tyr75 OH–H···O24 2.74, 158.9° |
| 20% ethanol | 125.4 | 12.70 | Arg208 NH2–H···O24 2.77, 168.2°; Arg238 NH2–H···O21⁻ 2.58, 167.8° |
| 50% ethanol | 141.8 | 4.45 | Asn253 ND2–H···O24 2.68, 156.7°; ZQG N36–H···Asn276 O 3.06, 162.1° |

The 2D ligand drawing is a chemically connected layout, not a projection of the 3D pose. Protein label positions and dashed line lengths in the 2D diagram are schematic; all annotated distances are calculated from the corresponding 3D snapshots. The 3D panels use a common viewing rotation after aligning the protein backbones. The 20% frame shows a displaced ligand; the SG–C28 dashed line represents a nonbonded separation, not an interaction or hydrogen bond.

These selected conformations illustrate mechanisms consistent with the time-series analyses. A single frame does not establish hydrogen-bond occupancy, catalytic rate, or binding free energy. Do not present the Cα trace artwork as an actual PyMOL cartoon rendering; for that visual style, use the included PyMOL scene scripts and inspect the generated PNGs.
