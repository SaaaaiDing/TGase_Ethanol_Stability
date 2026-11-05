#!/usr/bin/env python3
"""Neutralize packed 1IU4 in ethanol/water and prepare 160 ps GUI relaxation.

Run with $SCHRODINGER/run python3 from a new working directory.
"""

import argparse
import itertools
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from schrodinger.application.desmond import cms, constants
from schrodinger.application.matsci import desmondutils
from schrodinger.structutils import analyze


for candidate in (Path(__file__).resolve().parent,
                  Path(__file__).resolve().parent.parent):
    if (candidate / "expand_1IU4_EtOH20.py").is_file():
        sys.path.insert(0, str(candidate))
        break
else:
    raise RuntimeError("Cannot find expand_1IU4_EtOH20.py")

import expand_1IU4_EtOH20 as helper
import workflow_1IU4_EtOH20 as reference


PACKED = dict(water=19416, ethanol=5584, chloride=0, solute=2)
NEUTRAL = dict(water=19415, ethanol=5584, chloride=1, solute=2)


def check(st, expected, charge, atoms):
    groups = helper.groups(st)
    counts = {key: len(value) for key, value in groups.items()}
    helper.require(counts == expected, f"Unexpected molecules: {counts}")
    helper.require(st.formal_charge == charge,
                   f"Unexpected formal charge: {st.formal_charge}")
    helper.require(st.atom_total == atoms,
                   f"Unexpected atom count: {st.atom_total}")
    solutes = sorted((len(ids), sum(st.atom[i].formal_charge for i in ids))
                     for ids in groups["solute"])
    helper.require(solutes == [(42, -1), (5166, 2)],
                   f"Unexpected protein or ligand: {solutes}")
    print("Composition:", counts, flush=True)
    print("Atoms:", st.atom_total, "Formal charge:", st.formal_charge,
          flush=True)
    return groups


def find_ion_site(st, groups):
    lengths = helper.box(st)
    waters = groups["water"]
    solute_ids = [i for mol in groups["solute"] for i in mol
                  if st.atom[i].element != "H"]
    solute_xyz = np.array([[st.atom[i].x, st.atom[i].y, st.atom[i].z]
                           for i in solute_ids])
    oxygens = [next(i for i in mol if st.atom[i].element == "O")
               for mol in waters]
    oxygen_xyz = np.array([[st.atom[i].x, st.atom[i].y, st.atom[i].z]
                           for i in oxygens])
    solute_tree = cKDTree(np.mod(solute_xyz, lengths), boxsize=lengths)
    solute_distance = solute_tree.query(np.mod(oxygen_xyz, lengths))[0]
    all_atoms = list(st.atom)
    all_xyz = np.array([[a.x, a.y, a.z] for a in all_atoms])
    all_tree = cKDTree(np.mod(all_xyz, lengths), boxsize=lengths)
    near_d, near_i = all_tree.query(np.mod(oxygen_xyz, lengths), k=8)
    choices = []
    for n, water in enumerate(waters):
        if solute_distance[n] < 8.0:
            continue
        own = set(water)
        clearance = min((float(d) for d, j in zip(near_d[n], near_i[n])
                         if all_atoms[int(j)].index not in own), default=0.0)
        if clearance >= 2.5:
            choices.append((clearance, n))
    helper.require(bool(choices), "No suitable bulk-water site for chloride.")
    clearance, selected = max(choices)
    print("Selected water:", waters[selected], flush=True)
    print("Ion distance to solute (A):", float(solute_distance[selected]),
          flush=True)
    print("Nearest retained atom (A):", clearance, flush=True)
    return set(waters[selected]), oxygen_xyz[selected]


def make_clean_structure(st, removed, ion_xyz):
    xyz = np.array([[a.x, a.y, a.z] for a in st.atom])
    clean = helper.clean_copy(st,
                              [a.index for a in st.atom if a.index not in removed],
                              xyz)
    ion = clean.addAtom("Cl", *ion_xyz)
    ion.formal_charge = -1
    ion.pdbres = "CL"
    ion.pdbname = "CL"
    ion.resnum = 1
    ion.chain = "I"
    clean.title = "1IU4 in approximately 50 percent v/v ethanol, neutralized"
    for axis in "abc":
        for coord in "xyz":
            field = f"r_chorus_box_{axis}{coord}"
            clean.property[field] = st.property[field]
    clean.retype()
    clean.write("neutralized_clean.mae")
    check(clean, NEUTRAL, 0, 113710)
    return clean


def assign_forcefield(clean, original_box):
    built = desmondutils.run_system_builder(
        clean, "neutralized_50", forcefield="OPLS_2005",
        water_fftype="SPC", split_components=True, rezero_system=False,
        check_infinite=False, copy_output=False, set_incorporation=False)
    helper.require(bool(built) and Path(built).is_file(),
                   "Force-field assignment failed.")
    shutil.copy2(built, "neutralized.cms")
    model = cms.Cms(file="neutralized.cms")
    roles = dict(water=constants.CT_TYPE.VAL.SOLVENT,
                 ethanol=constants.CT_TYPE.VAL.COSOLVENT,
                 chloride=constants.CT_TYPE.VAL.ION,
                 solute=constants.CT_TYPE.VAL.SOLUTE)
    for component in model.comp_ct:
        present = [key for key, molecules in helper.groups(component).items()
                   if molecules]
        helper.require(len(present) == 1,
                       "A force-field component contains different species.")
        component.property[constants.CT_TYPE] = roles[present[0]]
    model.write("neutralized.cms")
    final = helper.read_full("neutralized.cms")
    groups = check(final, NEUTRAL, 0, 113710)
    helper.require(np.allclose(helper.box(final), original_box, atol=1e-4),
                   "Force-field assignment changed the box.")
    helper.forcefield_audit("neutralized.cms")
    return final, groups


def make_relaxation(gui_msj, st, groups):
    source = Path(gui_msj).read_text()
    positions = list(re.finditer(r"(?m)^\s*simulate\s*\{", source))
    helper.require(len(positions) == 6, "Expected five GUI relaxation stages and production.")
    durations = []
    for i in range(5):
        part = source[positions[i].start():positions[i + 1].start()]
        match = re.search(r"(?m)^\s*time\s*=\s*([0-9.]+)\b", part)
        helper.require(match is not None, f"Missing time in GUI stage {i + 1}.")
        durations.append(float(match.group(1)))
    helper.require(durations == [100., 12., 12., 12., 24.],
                   f"Unexpected GUI relaxation durations: {durations}")
    text = source[:positions[5].start()]
    ids = sorted(i for mol in groups["solute"] for i in mol)
    ranges = []
    for _, segment in itertools.groupby(enumerate(ids),
                                        lambda pair: pair[1] - pair[0]):
        indices = [pair[1] for pair in segment]
        ranges.append(f"atom.num {indices[0]}-{indices[-1]}")
    asl = "(" + " or ".join(ranges) + ") and not atom.ele H"
    heavy = {i for i in ids if st.atom[i].element != "H"}
    helper.require(set(analyze.evaluate_asl(st, asl)) == heavy,
                   "Solute restraint selection mismatch.")
    text, count = re.subn(r"\batom\s*=\s*(?:\"[^\"]*\"|'[^']*'|[^\s}]+)",
                          lambda _: f'atom = "{asl}"', text)
    helper.require(count == 4,
                   f"Expected four GUI restraint selections; found {count}.")
    updated_positions = list(re.finditer(r"(?m)^\s*simulate\s*\{", text))
    helper.require(len(updated_positions) == 5, "Unexpected GUI stages after editing.")
    final_stage = text[updated_positions[4].start():]
    for name in ("jobname", "dir", "trajectory", "maeff_output.trjdir"):
        helper.require(not re.search(r"(?m)^\s*" + re.escape(name)
                                     + r"\s*=", final_stage),
                       f"GUI stage already sets {name}.")
    closing = text.rfind("}")
    helper.require(closing >= updated_positions[4].start(), "Invalid GUI MSJ.")
    text = (text[:closing]
            + '\n jobname = "$MAINJOBNAME"\n dir = "."\n'
            + ' trajectory = { name = "equil_trj" interval = 1.0 }\n'
            + ' maeff_output.trjdir = "equil_trj"\n'
            + text[closing:])
    Path("equilibration_160ps.msj").write_text(text)
    print("Relaxation durations (ps):", durations, "Total:", sum(durations),
          flush=True)


def prepare(args):
    st = helper.read_full(args.input_cms)
    groups = check(st, PACKED, 1, 113712)
    original_box = helper.box(st)
    if Path("neutralized.cms").exists():
        final = helper.read_full("neutralized.cms")
        final_groups = check(final, NEUTRAL, 0, 113710)
        helper.require(np.allclose(helper.box(final), original_box, atol=1e-4),
                       "Existing neutralized CMS has a different box.")
        helper.forcefield_audit("neutralized.cms")
        print("Reusing validated neutralized.cms", flush=True)
    else:
        if Path("neutralized_clean.mae").exists():
            clean = helper.read_full("neutralized_clean.mae")
            check(clean, NEUTRAL, 0, 113710)
            helper.require(np.allclose(helper.box(clean), original_box, atol=1e-4),
                           "Existing clean MAE has a different box.")
            print("Reusing validated neutralized_clean.mae", flush=True)
        else:
            removed, ion_xyz = find_ion_site(st, groups)
            clean = make_clean_structure(st, removed, ion_xyz)
        final, final_groups = assign_forcefield(clean, original_box)
    make_relaxation(args.gui_msj, final, final_groups)
    os.environ["REFERENCE_CFG"] = str(Path(args.reference_cfg).resolve())
    os.environ["PRODUCTION_NS"] = "200"
    settings = reference.make_production_config()
    helper.require(settings["production_ns"] == 200.0,
                   "Expected 200 ns production.")
    Path("production.msj").write_text('''task { task = "desmond:auto" }
simulate {
 title = "Production MD"
 cfg_file = "production.cfg"
 jobname = "$MAINJOBNAME"
 dir = "."
 compress = ""
}
''')
    print("PREPARATION PASSED. Start the 160 ps GUI relaxation.", flush=True)


def check_relaxation(_args):
    final = helper.read_full("equilibrated.cms")
    groups = check(final, NEUTRAL, 0, 113710)
    helper.forcefield_audit("equilibrated.cms")
    image_gap, diameter = helper.solute_geometry(final, groups)
    lengths = helper.box(final)
    density = float(final.total_weight * 1.66053906660 / np.prod(lengths))
    print("Final box (A):", lengths.tolist(), flush=True)
    print("Final density (g/cm3):", density, flush=True)
    print("Solute periodic-image gap (A):", image_gap, flush=True)
    print("Solute diameter (A):", diameter, flush=True)
    helper.require(image_gap > 15.0, "The periodic-image gap is too small.")
    print("GEOMETRY PASSED. Review the final density and volume trend before MD.",
          flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)
    p = subcommands.add_parser("prepare")
    p.add_argument("--input-cms", required=True)
    p.add_argument("--gui-msj", required=True)
    p.add_argument("--reference-cfg", required=True)
    subcommands.add_parser("check")
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args)
    else:
        check_relaxation(args)
