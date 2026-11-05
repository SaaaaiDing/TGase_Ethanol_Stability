#!/usr/bin/env python3
"""Expand the existing neutral 1IU4 ethanol-water system with Schrodinger 2023-4.

Run with $SCHRODINGER/run python3. The build command never launches MD.
"""
import argparse
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist
from schrodinger import structure
from schrodinger.application.desmond import cms, constants
from schrodinger.application.matsci import desmondutils
from schrodinger.structutils import analyze

OLD_COUNTS = dict(water=9295, ethanol=704, chloride=1, solute=2)
NEW_COUNTS = dict(water=27885, ethanol=2112, chloride=1, solute=2)


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def read_full(path):
    with structure.StructureReader(str(path)) as reader:
        return max(reader, key=lambda st: st.atom_total)


def groups(st):
    result = {key: [] for key in OLD_COUNTS}
    for mol in st.molecule:
        atoms = list(mol.atom)
        formula = Counter(a.element for a in atoms)
        if formula == {'H': 2, 'O': 1}:
            kind = 'water'
        elif formula == {'H': 6, 'C': 2, 'O': 1}:
            kind = 'ethanol'
        elif formula == {'Cl': 1} and atoms[0].formal_charge == -1:
            kind = 'chloride'
        else:
            kind = 'solute'
        result[kind].append([a.index for a in atoms])
    return result


def box(st):
    matrix = np.array([[st.property[f'r_chorus_box_{a}{c}'] for c in 'xyz']
                       for a in 'abc'], dtype=float)
    require(np.allclose(matrix, np.diag(np.diag(matrix)), atol=1e-5),
            'Expected an orthorhombic box.')
    lengths = np.diag(matrix).copy()
    require(np.all(lengths > 0), 'Invalid box dimensions.')
    return lengths


def audit(st, expected):
    g = groups(st)
    counts = {key: len(value) for key, value in g.items()}
    require(counts == expected, f'Unexpected composition: {counts}')
    solutes = sorted((len(ids), sum(st.atom[i].formal_charge for i in ids))
                     for ids in g['solute'])
    require(solutes == [(42, -1), (5166, 2)], f'Unexpected solutes: {solutes}')
    require(st.formal_charge == 0, 'The system is not neutral.')
    print('Composition:', counts, flush=True)
    print('Atoms:', st.atom_total, 'Formal charge:', st.formal_charge, flush=True)
    return g


def whole_centered(st, g):
    lengths = box(st)
    raw = np.array([[a.x, a.y, a.z] for a in st.atom])
    xyz = np.full_like(raw, np.nan)
    adj = [[] for _ in st.atom]
    for b in st.bond:
        i, j = b.atom1.index - 1, b.atom2.index - 1
        adj[i].append(j)
        adj[j].append(i)
    all_molecules = [ids for value in g.values() for ids in value]
    for ids in all_molecules:
        root = ids[0] - 1
        xyz[root] = raw[root]
        stack = [root]
        while stack:
            i = stack.pop()
            for j in adj[i]:
                delta = raw[j] - raw[i]
                delta -= lengths * np.rint(delta / lengths)
                candidate = xyz[i] + delta
                if np.isnan(xyz[j, 0]):
                    xyz[j] = candidate
                    stack.append(j)
                else:
                    require(np.linalg.norm(xyz[j] - candidate) < 1e-3,
                            'Inconsistent bond unwrapping.')
    protein = np.array(max(g['solute'], key=len)) - 1
    center = 0.5 * (xyz[protein].min(0) + xyz[protein].max(0))
    for ids in all_molecules:
        idx = np.array(ids) - 1
        shift = lengths * np.rint((xyz[idx].mean(0) - center) / lengths)
        xyz[idx] -= shift + center
    require(np.isfinite(xyz).all(), 'Invalid coordinates after unwrapping.')
    return xyz


def clean_copy(st, ids, xyz):
    out = structure.create_new_structure()
    mapping = {}
    for i in sorted(ids):
        old = st.atom[i]
        new = out.addAtom(old.element, *xyz[i - 1])
        for field in ('formal_charge', 'pdbres', 'pdbname', 'resnum', 'chain', 'inscode'):
            setattr(new, field, getattr(old, field))
        mapping[i] = new.index
    for b in st.bond:
        i, j = b.atom1.index, b.atom2.index
        if i in mapping and j in mapping:
            out.addBond(mapping[i], mapping[j], b.order)
    out.retype()
    return out


def solute_geometry(st, g):
    xyz = whole_centered(st, g)
    ids = sorted(i for mol in g['solute'] for i in mol if st.atom[i].element != 'H')
    points = xyz[np.array(ids) - 1]
    lengths = box(st)
    tree = cKDTree(points)
    gap = min(float(tree.query(points + np.array(s) * lengths)[0].min())
              for s in itertools.product((-1, 0, 1), repeat=3) if s != (0, 0, 0))
    return gap, float(pdist(points).max())


def forcefield_audit(path):
    model = cms.Cms(file=str(path))
    require(desmondutils.is_opls2005_model(model), 'Expected OPLS_2005.')
    charge = 0.0
    for ct in model.comp_ct:
        sites = list(ct.ffio.site)
        size = sum(s.type == 'atom' for s in sites)
        require(size > 0 and ct.atom_total % size == 0, 'Invalid component template.')
        charge += (ct.atom_total // size) * sum(s.charge for s in sites)
    require(abs(charge) < 1e-3, f'Nonzero force-field charge: {charge}')
    print('Force-field net charge:', charge, flush=True)


def write_md_files(source_dir, st, g):
    text = (source_dir / 'equilibration_2023_4.msj').read_text()
    require(not re.search(r'(?m)^\s*minimize\s*\{', text), 'Legacy minimize stage found.')
    require(len(re.findall(r'(?m)^simulate\s*\{', text)) == 6 and 'Brownie' in text,
            'Expected the previously successful six-stage relaxation protocol.')
    ids = sorted(i for mol in g['solute'] for i in mol)
    ranges = []
    for _, block in itertools.groupby(enumerate(ids), lambda pair: pair[1] - pair[0]):
        values = [pair[1] for pair in block]
        ranges.append(f'atom.num {values[0]}-{values[-1]}')
    asl = '(' + ' or '.join(ranges) + ') and not atom.ele H'
    expected = {i for i in ids if st.atom[i].element != 'H'}
    require(set(analyze.evaluate_asl(st, asl)) == expected, 'Restraint selection mismatch.')
    text, n = re.subn(r'\batom\s*=\s*"[^"]*"', lambda m: f'atom = "{asl}"', text)
    require(n == 5, f'Expected five restraint selections; found {n}.')
    Path('equilibration_large.msj').write_text(text)
    for name in ('production.msj', 'production.cfg', 'settings.json'):
        shutil.copy2(source_dir / name, name)


def build(args):
    source = Path(args.input_cms).resolve()
    source_dir = source.parent
    require(source.is_file(), f'Missing input: {source}')
    for name in ('equilibration_2023_4.msj', 'production.msj', 'production.cfg', 'settings.json'):
        require((source_dir / name).is_file(), f'Missing reference file: {source_dir / name}')
    require(not any(Path(name).exists() for name in
                    ('fixed_old_system.maegz', 'expanded.cms', 'expanded_mix_system-out.cms')),
            'Build outputs already exist. Use a new empty working directory.')
    st = read_full(source)
    g = audit(st, OLD_COUNTS)
    forcefield_audit(source)
    xyz = whole_centered(st, g)
    retained = clean_copy(st, range(1, st.atom_total + 1), xyz)
    retained.title = 'Existing neutral 1IU4 ethanol-water system'
    retained.write('fixed_old_system.maegz')
    with structure.StructureWriter('solvent_templates.maegz') as writer:
        for kind in ('ethanol', 'water'):
            ids = g[kind][0]
            local = xyz.copy()
            local[np.array(ids) - 1] -= local[np.array(ids) - 1].mean(0)
            template = clean_copy(st, ids, local)
            template.title = kind
            writer.append(template)
    old_gap, diameter = solute_geometry(st, g)
    solvent_mass = sum(st.atom[i].atomic_weight for kind in ('water', 'ethanol')
                       for mol in g[kind] for i in mol)
    target_mass = st.total_weight + 2 * solvent_mass
    estimated_edge = (target_mass * 1.66053906660 / 1.004055683309238) ** (1 / 3)
    require(estimated_edge - diameter > 20.0,
            'The proposed solvent count does not provide the planned geometric margin.')
    plan = dict(input_cms=str(source), old_counts=OLD_COUNTS, target_counts=NEW_COUNTS,
                add_ethanol=1408, add_water=18590, solute_diameter_A=diameter,
                estimated_equilibrated_cubic_edge_A=estimated_edge,
                initial_density_setting_g_cm3=0.5,
                note='The box estimate is not a substitute for post-equilibration checks.')
    Path('expansion_plan.json').write_text(json.dumps(plan, indent=2))
    print('Expansion plan:', plan, flush=True)
    runner = Path(os.environ['SCHRODINGER']) / 'run'
    command = [str(runner), 'disordered_system_builder_gui_dir/disordered_system_builder_driver.py',
               '-molecules', '19998', '-composition', '1408:18590',
               '-substrate', 'fixed_old_system.maegz', '-substrate_mode', 'immersed',
               '-pbc', 'new', '-pack', '-density', '0.500', '-scale_vdw', '0.50',
               '-obey', 'density', '-tries_per_mol', '50', '-tries_per_dc', '20',
               '-seed', '2345', '-preserve_res_info', '-water_fftype', 'SPC',
               '-split_components', '-forcefield', 'OPLS_2005',
               'solvent_templates.maegz', 'expanded_mix']
    print('Starting Disordered System Builder; detailed output: packing.log', flush=True)
    with open('packing.log', 'w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    result = Path('expanded_mix_system-out.cms')
    require(result.is_file(), 'No CMS produced. Inspect packing.log before continuing.')
    model = cms.Cms(file=str(result))
    roles = dict(water=constants.CT_TYPE.VAL.SOLVENT,
                 ethanol=constants.CT_TYPE.VAL.COSOLVENT,
                 chloride=constants.CT_TYPE.VAL.ION,
                 solute=constants.CT_TYPE.VAL.SOLUTE)
    for ct in model.comp_ct:
        present = [kind for kind, mols in groups(ct).items() if mols]
        require(len(present) == 1, 'Mixed molecule types in a force-field component.')
        ct.property[constants.CT_TYPE] = roles[present[0]]
    model.write('expanded.cms')
    final = read_full('expanded.cms')
    fg = audit(final, NEW_COUNTS)
    require(final.atom_total == 107872, 'Unexpected atom count.')
    forcefield_audit('expanded.cms')
    initial_gap, new_diameter = solute_geometry(final, fg)
    old_ids = sorted(i for mol in g['solute'] for i in mol)
    new_ids = sorted(i for mol in fg['solute'] for i in mol)
    p = xyz[np.array(old_ids) - 1]
    q = whole_centered(final, fg)[np.array(new_ids) - 1]
    p = p - p.mean(0)
    q = q - q.mean(0)
    u, _, vt = np.linalg.svd(p.T @ q)
    correction = np.diag([1.0, 1.0, np.linalg.det(u @ vt)])
    rmsd = float(np.sqrt(np.mean(np.sum((p @ (u @ correction @ vt) - q) ** 2, axis=1))))
    require(rmsd < 0.01, f'Solute geometry/order changed during packing: RMSD={rmsd}')
    require(initial_gap > 20.0, f'Insufficient initial image gap: {initial_gap}')
    print('Initial box (angstrom):', box(final).tolist(), flush=True)
    print('Initial solute image distance (angstrom):', initial_gap, flush=True)
    print('Solute RMSD after rigid alignment (angstrom):', rmsd, flush=True)
    write_md_files(source_dir, final, fg)
    print('BUILD PASSED. expanded.cms and equilibration_large.msj are ready.', flush=True)
    print('Re-equilibrate the expanded system before production MD.', flush=True)


def check(args):
    from schrodinger.application.desmond.packages import traj
    st = read_full('equilibrated.cms')
    g = audit(st, NEW_COUNTS)
    forcefield_audit('equilibrated.cms')
    frames = traj.read_traj('equil_trj')
    require(len(frames) >= 101, 'Equilibration trajectory is too short.')
    selected = [f for f in frames if f.time >= frames[-1].time - 2000.0]
    t = np.array([f.time for f in selected])
    v = np.array([abs(np.linalg.det(f.box)) for f in selected])
    require(len(v) >= 100 and np.isfinite(v).all() and np.all(v > 0), 'Invalid volume data.')
    rho = st.total_weight * 1.66053906660 / v
    mid = len(v) // 2
    block = abs(v[:mid].mean() - v[mid:].mean()) / v.mean()
    drift = abs(np.polyfit(t-t[0], v, 1)[0]) * (t[-1]-t[0]) / v.mean()
    gap, diameter = solute_geometry(st, g)
    settings = json.loads(Path('settings.json').read_text())
    report = dict(mean_density_g_cm3=float(rho.mean()), volume_block_change_fraction=float(block),
                  volume_drift_fraction=float(drift), final_solute_image_distance_A=gap,
                  solute_diameter_A=diameter, box_A=box(st).tolist(),
                  note='Volume screen plus a final-frame geometry check; not proof of full equilibration.')
    Path('expanded_equilibration_check.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)
    require(0.85 < rho.mean() < 1.30 and block < 0.02 and drift < 0.02,
            'Density/volume requires further equilibration.')
    require(gap > max(float(settings['cutoff_A']), 15.0),
            'Final image separation is below the planned minimum margin.')
    print('BASIC EQUILIBRATION CHECKS PASSED.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('build')
    p.add_argument('--input-cms', required=True)
    commands.add_parser('check')
    args = parser.parse_args()
    if args.command == 'build':
        build(args)
    else:
        check(args)
