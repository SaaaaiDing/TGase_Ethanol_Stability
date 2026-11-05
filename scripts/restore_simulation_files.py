"""Restore deposited CMS/configuration files into a fresh raw-data layout."""
import argparse
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--output-root', required=True, type=Path)
a = p.parse_args()
out = a.output_root.expanduser().resolve()
if out.exists() and any(out.iterdir()):
    raise SystemExit('Choose a new or empty output directory; existing data will not be overwritten.')
source = ROOT / 'simulation/raw_layout'
for path in sorted(source.rglob('*')):
    if not path.is_file(): continue
    relative = path.relative_to(source)
    data = path.read_bytes()
    if path.name.endswith('.cms.gz'):
        relative = relative.with_suffix('')
        data = gzip.decompress(data)
    target = out / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
print(f'Restored simulation files to {out}. Full trajectory directories must be supplied separately.')
