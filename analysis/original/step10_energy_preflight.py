#!/usr/bin/env python3
"""Inspect the installed Desmond energy API and the three simulation systems."""

import inspect
import os
from collections import Counter
from pathlib import Path

from schrodinger.application.desmond.packages import topo


analysis_dir = Path(os.environ["ANALYSIS_DIR"]).expanduser().resolve()
root = analysis_dir.parent
systems = {
    "water": root / "1IU4_0/desmond_md_job_4/1IU4_md_batch-out.cms",
    "ethanol20": root / "1IU4_0_20Eoth/expanded_box_iI2MQA/production_55774175/1IU4_EtOH20_md_55774175-out.cms",
    "ethanol50": root / "1IU4_0_20Eoth/1IU4_EtOH50/neutralization_50/npt_500ps/production_55831542/1IU4_EtOH50_md_55831542-out.cms",
}

try:
    from schrodinger.application.desmond.packages import energygroup
except ImportError as exc:
    print(f"Desmond energygroup import failed: {exc}")
else:
    print(f"Desmond energygroup module: {energygroup.__file__}")
    for name in ("EnergyFacade", "InteractionEnergy", "SliceParams", "get_energies"):
        obj = getattr(energygroup, name, None)
        print(f"{name}: {inspect.signature(obj) if obj is not None else 'UNAVAILABLE'}")

for label, path in systems.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    _, cms = topo.read_cms(str(path))
    ct = cms.fsys_ct
    residue_atoms = Counter(atom.pdbres.strip() for atom in ct.atom)
    print(f"\n{label}: {path}")
    print(f"Total atoms: {ct.atom_total}")
    print(f"Residue names by atom count: {residue_atoms.most_common(12)}")
    print(f"Protein ASL atoms: {len(cms.select_atom('protein'))}")
    for name in ("LIG1", "UNK", "HOH", "WAT", "SPC"):
        try:
            found = cms.select_atom(f"res.ptype {name}")
        except Exception as exc:
            print(f"ASL res.ptype {name}: ERROR {type(exc).__name__}: {exc}")
        else:
            print(f"ASL res.ptype {name}: {len(found)} atoms")
    cfgs = sorted(path.parent.glob("*.cfg"))
    print("Configuration files:")
    for cfg in cfgs[:15]:
        print(f"  {cfg.name}: {cfg.stat().st_size} bytes")
    if not cfgs:
        print("  NONE")

print("\nRead-only preflight complete; no energies were calculated.")
