"""Validate publication integrity and reproduce key summaries without running MD."""
import ast
import gzip
import hashlib
import importlib.metadata
import importlib.util
import json
import platform
import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / 'data/tables'
def sha(data): return hashlib.sha256(data).hexdigest()
manifest = json.loads((ROOT/'provenance/source_manifest.json').read_text())
for item in manifest:
    path = ROOT / item['published_path']; data = path.read_bytes()
    assert sha(data) == item['published_sha256'], f'Published hash changed: {path}'
    if path.name.endswith('.gz'):
        assert sha(gzip.decompress(data)) == item['source_sha256'], f'Compressed content changed: {path}'
    if item['transformation'] == 'none': assert item['source_sha256'] == item['published_sha256']
for path in ROOT.rglob('*.py'): ast.parse(path.read_text(), filename=str(path))

spec = importlib.util.spec_from_file_location('comparison', ROOT/'analysis/make_1IU4_md_comparison.py')
comparison = importlib.util.module_from_spec(spec); spec.loader.exec_module(comparison)
energies, frames, solvent = comparison.load_inputs(TABLES)
state, pairs = comparison.build_tables(energies, frames, solvent)
for actual, name in [(state, '1IU4_state_and_solvent_comparison.csv'), (pairs, '1IU4_residue_pair_energy_comparison.csv')]:
    expected = pd.read_csv(TABLES/name)
    pd.testing.assert_frame_equal(actual.reset_index(drop=True), expected, check_dtype=False, rtol=1e-9, atol=1e-10)

summary = pd.read_csv(TABLES/'reactive_geometry_summary_50_200ns.csv').set_index('system')
records = []
for system, expected_near in [('water',368), ('ethanol20',22), ('ethanol50',695)]:
    geometry = pd.read_csv(TABLES/f'{system}_reactive_geometry.csv')
    rmsd = pd.read_csv(TABLES/f'{system}_rmsd.csv')
    assert len(geometry) == len(rmsd) == 1002
    assert np.all(np.diff(geometry.time_ns) > 0) and abs(geometry.time_ns.iloc[-1]-200) < .01
    assert np.allclose(geometry.time_ns, rmsd.time_ns, atol=1e-6)
    late = geometry.loc[geometry.time_ns >= 50]
    near = late.Cys64_SG_to_ZQG_C28_A.le(5)
    assert len(late) == 751 and int(near.sum()) == expected_near
    assert np.isclose(near.mean(), summary.loc[system,'fraction_SG_C28_le_5A_50_to_200_ns'])
    sampled = state.loc[(state.system==system) & (state.period_ns=='50-200') & (state.SG_C28_state=='near')].iloc[0]
    assert int(sampled.n_period_frames) == 150
    records.append({'system':system, 'frames':len(geometry), 'terminal_time_ns':float(geometry.time_ns.iloc[-1]),
                    'dense_late_frames':len(late), 'dense_near_frames':int(near.sum()),
                    'sampled_late_near_frames':int(sampled.n_sampled_frames), 'sampled_near_hbond_count':int(sampled.Cys64_NH_O22_count)})

for document in [ROOT/'README.md', ROOT/'docs/METHODS.md', ROOT/'docs/REPRODUCIBILITY.md']:
    prose = re.sub(r'```.*?```|`[^`]*`', '', document.read_text(), flags=re.S)
    for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', prose):
        if '://' in target or target.startswith('#'): continue
        assert (document.parent/target.split('#')[0]).exists(), f'Broken link in {document}: {target}'

secrets = re.compile(r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----|\bgh[pousr]_[A-Za-z0-9]{30,}|\bgithub_pat_[A-Za-z0-9_]{30,}|\bsk-[A-Za-z0-9]{35,}')
for path in ROOT.rglob('*'):
    if path.is_file() and path.suffix in ['.py','.md','.json','.slurm','.cfg','.msj','.txt']:
        assert not secrets.search(path.read_text(errors='replace')), f'Credential-like content: {path}'
report = {'status':'passed', 'imported_files_verified':len(manifest), 'key_summary_tables_reproduced':2,
          'trajectory_observables':records, 'python':platform.python_version(),
          'packages':{name:importlib.metadata.version(name) for name in ['numpy','pandas','matplotlib','Pillow']},
          'scope':'File integrity, code syntax and deposited-table consistency. No new MD, trajectory recalculation or PyMOL rerendering.'}
(ROOT/'provenance/validation_report.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
