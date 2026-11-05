from schrodinger.application.desmond.config import DEFAULT_SETTING
#!/usr/bin/env python3
"""Run with the Schrodinger Python interpreter on the allocated compute node."""
import itertools
import json
import math
import os
import shutil
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from schrodinger import structure
from schrodinger.application.desmond import cms, constants
from schrodinger.application.matsci import desmondutils
from schrodinger.structutils import analyze
from schrodinger.utils import sea


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def read_full(path):
    with structure.StructureReader(str(path)) as reader:
        return max(reader, key=lambda st: st.atom_total)


def classify(st):
    groups = {name: [] for name in ('water', 'ethanol', 'chloride', 'solute')}
    for mol in st.molecule:
        atoms = list(mol.atom)
        elements = Counter(atom.element for atom in atoms)
        if elements == {'H': 2, 'O': 1}:
            name = 'water'
        elif elements == {'C': 2, 'H': 6, 'O': 1}:
            name = 'ethanol'
        elif elements == {'Cl': 1} and atoms[0].formal_charge == -1:
            name = 'chloride'
        else:
            name = 'solute'
        groups[name].append(atoms)
    return groups


def box_matrix(st):
    return np.array([[st.property[f'r_chorus_box_{axis}{coord}']
                      for coord in 'xyz'] for axis in 'abc'], dtype=float)


def audit(path, neutral):
    st = read_full(path)
    groups = classify(st)
    expected = {'water': 9295 if neutral else 9296,
                'ethanol': 704, 'chloride': 1 if neutral else 0, 'solute': 2}
    actual = {key: len(value) for key, value in groups.items()}
    print('Composition:', actual, flush=True)
    require(actual == expected, f'Unexpected composition: {actual}')
    solutes = sorted((len(atoms), sum(a.formal_charge for a in atoms))
                     for atoms in groups['solute'])
    require(solutes == [(42, -1), (5166, 2)],
            f'Unexpected protein/ligand atom counts or charges: {solutes}')
    require(st.formal_charge == (0 if neutral else 1),
            f'Unexpected formal charge: {st.formal_charge}')
    require(st.atom_total == (39430 if neutral else 39432),
            f'Unexpected total atom count: {st.atom_total}')
    box = box_matrix(st)
    require(np.allclose(box, np.diag(np.diag(box)), atol=1e-6),
            'This workflow expects an orthorhombic box.')
    require(np.all(np.diag(box) > 0), 'Invalid box dimensions.')
    print('Atoms:', st.atom_total, 'Formal charge:', st.formal_charge)
    print('Box (angstrom):', np.diag(box).tolist(), flush=True)
    return st, groups, box


def audit_forcefield(path):
    model = cms.Cms(file=str(path))
    total_charge = 0.0
    require(desmondutils.is_opls2005_model(model), 'Expected OPLS_2005.')
    for ct in model.comp_ct:
        sites = list(ct.ffio.site)
        n_atoms = sum(site.type == 'atom' for site in sites)
        require(n_atoms > 0 and ct.atom_total % n_atoms == 0,
                'Inconsistent component/template atom counts in rebuilt CMS.')
        total_charge += (ct.atom_total // n_atoms) * sum(s.charge for s in sites)
    print('Force-field net charge:', total_charge, flush=True)
    require(abs(total_charge) < 1e-3, 'Force-field net charge is not zero.')


def make_production_config():
    text = Path(os.environ['REFERENCE_CFG']).read_text()
    cfg = sea.Map(text)
    temperature = cfg['temperature'].val
    while isinstance(temperature, (list, tuple)):
        temperature = temperature[0]
    temperature = float(temperature)
    pressure = cfg['pressure'].val
    while isinstance(pressure, (list, tuple)):
        pressure = pressure[0]
    pressure = float(pressure)
    require(cfg['ensemble']['class'].val == 'NPT',
            'Expected an NPT pure-water reference configuration.')
    require(not cfg['annealing'].val, 'Annealing reference CFG is not supported.')
    if 'restrain' in cfg:
        require(cfg['restrain'].val == 'none',
                'Reference production CFG contains restraints; review it first.')
    if 'restraints' in cfg:
        require(not cfg['restraints']['new'].val,
                'Reference CFG contains atom-specific production restraints.')
    require(200.0 < temperature < 400.0 and 0.0 < pressure < 10.0,
            'Reference temperature/pressure requires manual review.')
    if os.environ['PRODUCTION_NS']:
        duration_ns = float(os.environ['PRODUCTION_NS'])
        require(math.isfinite(duration_ns) and duration_ns > 0,
                'Production duration must be a positive finite number.')
        cfg['time'].val = duration_ns * 1000.0
    else:
        duration_ns = float(cfg['time'].val) / 1000.0
    # Bind the model and output paths to this run, not the old GUI job.
    for key in ('model_file', 'cfg_file'):
        if key in cfg:
            del cfg[key]
    cfg['trajectory']['name'].val = 'production_trj'
    cfg['maeff_output']['name'].val = '$JOBNAME-out.cms'
    cfg['maeff_output']['trjdir'] = 'production_trj'
    cfg['eneseq']['name'].val = '$JOBNAME.ene'
    Path('production.cfg').write_text(str(cfg))
    settings = {'temperature_K': temperature, 'pressure_bar': pressure,
                'production_ns': duration_ns, 'cutoff_A': float(cfg['cutoff_radius'].val if 'cutoff_radius' in cfg else DEFAULT_SETTING.MD['cutoff_radius'].val),
                'reference_cfg': os.environ['REFERENCE_CFG']}
    Path('settings.json').write_text(json.dumps(settings, indent=2))
    print('Production settings:', settings, flush=True)
    return settings


def prepare():
    settings = make_production_config()
    original, groups, box = audit(os.environ['INPUT_CMS'], neutral=False)
    lengths = np.diag(box)
    solute_heavy = np.array([[a.x, a.y, a.z] for mol in groups['solute']
                            for a in mol if a.element != 'H'])
    oxygen_atoms = [next(a for a in mol if a.element == 'O')
                    for mol in groups['water']]
    oxygen_xyz = np.array([[a.x, a.y, a.z] for a in oxygen_atoms])
    solute_tree = cKDTree(np.mod(solute_heavy, lengths), boxsize=lengths)
    solute_dist = solute_tree.query(np.mod(oxygen_xyz, lengths))[0]
    all_atoms = list(original.atom)
    all_xyz = np.array([[a.x, a.y, a.z] for a in all_atoms])
    all_tree = cKDTree(np.mod(all_xyz, lengths), boxsize=lengths)
    near_d, near_i = all_tree.query(np.mod(oxygen_xyz, lengths), k=8)
    candidates = []
    for i, mol in enumerate(groups['water']):
        if solute_dist[i] < 8.0:
            continue
        own_atoms = {a.index for a in mol}
        clearance = min(d for d, j in zip(near_d[i], near_i[i])
                        if all_atoms[int(j)].index not in own_atoms)
        if clearance >= 2.5:
            candidates.append((clearance, i))
    require(candidates, 'No suitable bulk-water site for chloride was found.')
    clearance, selected = max(candidates)
    removed = {a.index for a in groups['water'][selected]}
    ion_xyz = oxygen_xyz[selected]
    print('Replacing water with Cl-:', sorted(removed))
    print('Ion-to-solute distance (angstrom):', float(solute_dist[selected]))
    print('Nearest retained atom (angstrom):', float(clearance), flush=True)

    # Reconstruct atoms and bonds in a fresh CT. Do not copy FFIO/component metadata.
    clean = structure.create_new_structure()
    index_map = {}
    for old in original.atom:
        if old.index in removed:
            continue
        new = clean.addAtom(old.element, old.x, old.y, old.z)
        for field in ('formal_charge', 'pdbres', 'pdbname', 'resnum', 'chain', 'inscode'):
            setattr(new, field, getattr(old, field))
        index_map[old.index] = new.index
    for bond in original.bond:
        a, b = bond.atom1.index, bond.atom2.index
        if a in index_map and b in index_map:
            clean.addBond(index_map[a], index_map[b], bond.order)
    ion = clean.addAtom('Cl', *ion_xyz)
    ion.formal_charge = -1
    ion.pdbres = 'CL'
    ion.pdbname = 'CL'
    ion.resnum = 1
    ion.chain = 'I'
    clean.title = '1IU4 in approximately 20 percent v/v ethanol, neutralized'
    for axis in 'abc':
        for coord in 'xyz':
            key = f'r_chorus_box_{axis}{coord}'
            clean.property[key] = original.property[key]
    clean.retype()
    clean.write('neutralized_clean.mae')
    audit('neutralized_clean.mae', neutral=True)

    # Use the same Materials Science system-building API as the disorder workflow.
    built = desmondutils.run_system_builder(
        clean, 'neutralized', forcefield='OPLS_2005', water_fftype='SPC',
        split_components=True, rezero_system=False, check_infinite=False,
        copy_output=False, set_incorporation=False)
    require(built and Path(built).is_file(), 'Force-field assignment failed.')
    shutil.copy2(built, 'neutralized.cms')
    # The Materials Science builder labels components as solvent by default.
    # Restore biological component roles without changing their force fields.
    model = cms.Cms(file='neutralized.cms')
    for component in model.comp_ct:
        component_groups = classify(component)
        present = [key for key, value in component_groups.items() if value]
        require(len(present) == 1, 'A rebuilt component contains mixed molecule types.')
        roles = {'water': constants.CT_TYPE.VAL.SOLVENT,
                 'ethanol': constants.CT_TYPE.VAL.COSOLVENT,
                 'chloride': constants.CT_TYPE.VAL.ION,
                 'solute': constants.CT_TYPE.VAL.SOLUTE}
        component.property[constants.CT_TYPE] = roles[present[0]]
    model.write('neutralized.cms')
    st, final_groups, final_box = audit('neutralized.cms', neutral=True)
    require(np.allclose(final_box, box, atol=1e-4), 'The initial box changed unexpectedly.')
    audit_forcefield('neutralized.cms')

    # Restrain only the protein and bound ligand, never the ethanol molecules.
    indices = sorted(a.index for mol in final_groups['solute'] for a in mol)
    ranges = []
    for _, group in itertools.groupby(enumerate(indices), lambda pair: pair[1] - pair[0]):
        values = [pair[1] for pair in group]
        ranges.append(f'atom.num {values[0]}-{values[-1]}')
    asl = '(' + ' or '.join(ranges) + ') and not atom.ele H'
    expected_heavy = sum(a.element != 'H' for mol in final_groups['solute'] for a in mol)
    require(len(analyze.evaluate_asl(st, asl)) == expected_heavy, 'Restraint selection mismatch.')

    def restraints(fc):
        return 'none' if fc == 0 else f'{{ atom = "{asl}" force_constant = {fc} }}'

    stages = ['task { task = "desmond:auto" }',
              f'minimize {{ max_steps = 5000 steepest_descent_steps = 1000 '
              f'restrain = {restraints(50)} }}']
    # Duration is in ps; force constants are in kcal/mol/angstrom^2.
    schedule = [
        (100, 10.0, 'NVT', 'Brownie', 50),
        (100, 10.0, 'NVT', 'Langevin', 50),
        (100, settings['temperature_K'], 'NVT', 'Langevin', 50),
        (1000, settings['temperature_K'], 'NPT', 'Langevin', 5),
        (1000, settings['temperature_K'], 'NPT', 'Langevin', 1),
        (5000, settings['temperature_K'], 'NPT', 'MTK', 0),
    ]
    for i, (duration, temp, ensemble, method, fc) in enumerate(schedule):
        final = i == len(schedule) - 1
        step = '[0.002 0.002 0.006]' if final else '[0.001 0.001 0.003]'
        thermostat_tau = 1.0 if final else 0.1
        extra = ''
        if method == 'Brownie':
            extra += 'ensemble.brownie.delta_max = 0.1\n'
        if final:
            extra += 'jobname = "$MAINJOBNAME"\ndir = "."\n'
            extra += 'trajectory = { name = "equil_trj" interval = 10.0 }\n'
            extra += 'maeff_output.trjdir = "equil_trj"\n'
        stages.append(f'''simulate {{
 title = "Equilibration stage {i + 1}"
 time = {duration}
 timestep = {step}
 temperature = {temp}
 pressure = [{settings['pressure_bar']} isotropic]
 ensemble = {{ class = {ensemble} method = {method}
               thermostat.tau = {thermostat_tau} barostat.tau = 2.0 }}
 restrain = {restraints(fc)}
 eneseq.interval = 1.0
 checkpt.write_last_step = true
 compress = ""
 {extra}
}}''')
    Path('equilibration.msj').write_text('\n\n'.join(stages) + '\n')
    Path('production.msj').write_text('''task { task = "desmond:auto" }
simulate {
 title = "Production MD"
 cfg_file = "production.cfg"
 jobname = "$MAINJOBNAME"
 dir = "."
 compress = ""
}
''')


def check_equilibration():
    from schrodinger.application.desmond.packages import traj
    st, groups, box = audit('equilibrated.cms', neutral=True)
    audit_forcefield('equilibrated.cms')
    frames = traj.read_traj('equil_trj')
    require(len(frames) >= 101, 'Equilibration trajectory is missing or too short.')
    end = frames[-1].time
    selected = [frame for frame in frames if frame.time >= end - 2000.0]
    times = np.array([frame.time for frame in selected])
    volumes = np.array([abs(np.linalg.det(frame.box)) for frame in selected])
    require(len(selected) >= 100 and np.all(np.isfinite(volumes)) and np.all(volumes > 0),
            'Insufficient or invalid volume data.')
    densities = st.total_weight * 1.66053906660 / volumes
    midpoint = len(volumes) // 2
    block_change = abs(volumes[:midpoint].mean() - volumes[midpoint:].mean()) / volumes.mean()
    drift = abs(np.polyfit(times - times[0], volumes, 1)[0]) * (times[-1] - times[0]) / volumes.mean()
    settings = json.loads(Path('settings.json').read_text())
    solute_xyz = np.array([[a.x, a.y, a.z] for mol in groups['solute']
                          for a in mol if a.element != 'H'])
    tree = cKDTree(solute_xyz)
    image_distance = min(float(tree.query(solute_xyz + np.array(shift) @ box)[0].min())
                         for shift in itertools.product((-1, 0, 1), repeat=3)
                         if shift != (0, 0, 0))
    report = {'mean_density_g_cm3': float(densities.mean()),
              'volume_block_change_fraction': float(block_change),
              'volume_drift_fraction': float(drift),
              'solute_image_distance_A': image_distance,
              'note': 'A basic density/volume screen, not proof of full equilibration.'}
    Path('equilibration_check.json').write_text(json.dumps(report, indent=2))
    print('Equilibration check:', report, flush=True)
    require(0.85 < densities.mean() < 1.30 and block_change < 0.02 and drift < 0.02,
            'Density/volume has not stabilized. Equilibrated CMS is retained; extend NPT before production.')
    require(image_distance > settings['cutoff_A'],
            'Solute is within the short-range cutoff of its periodic image; enlarge the solvent box.')
    print('Basic checks passed. Proceeding to production MD.', flush=True)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Prepare and check the existing 1IU4 ethanol-water system.')
    commands = parser.add_subparsers(dest='command', required=True)
    prep_parser = commands.add_parser('prepare')
    prep_parser.add_argument('--input-cms', required=True)
    prep_parser.add_argument('--reference-cfg', required=True)
    prep_parser.add_argument('--production-ns', type=float, default=None)
    commands.add_parser('check')
    args = parser.parse_args()

    if args.command == 'prepare':
        for filename in (args.input_cms, args.reference_cfg):
            require(Path(filename).is_file(), f'Missing input: {filename}')
        generated = ('production.cfg', 'settings.json', 'neutralized_clean.mae',
                     'neutralized.cms', 'neutralized_system-out.cms',
                     'equilibration.msj', 'production.msj')
        require(not any(Path(name).exists() for name in generated),
                'Preparation outputs already exist. Use a new working directory.')
        os.environ['INPUT_CMS'] = str(Path(args.input_cms).resolve())
        os.environ['REFERENCE_CFG'] = str(Path(args.reference_cfg).resolve())
        os.environ['PRODUCTION_NS'] = (
            '' if args.production_ns is None else str(args.production_ns))
        prepare()
    else:
        check_equilibration()
