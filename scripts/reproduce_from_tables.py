"""Rebuild analysis plots and summary tables in a fresh directory."""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, default=Path('reproduced')); a=p.parse_args()
out=a.output.expanduser().resolve()
if out.exists() and any(out.iterdir()): raise SystemExit('Choose a new or empty output directory.')
out.mkdir(parents=True,exist_ok=True)
for path in (ROOT/'data/tables').glob('*.csv'): shutil.copy2(path,out/path.name)
env=os.environ.copy(); env['ANALYSIS_DIR']=str(out); env['MPLBACKEND']='Agg'
for name in ['step03_plot_timeseries.py','step05_reactive_geometry_figure.py']:
    subprocess.run([sys.executable,str(ROOT/'analysis'/name)],env=env,check=True)
subprocess.run([sys.executable,str(ROOT/'analysis/make_1IU4_md_comparison.py'),
                '--input-dir',str(ROOT/'data/tables'),'--output-dir',str(out)],env=env,check=True)
for name in ['reactive_geometry_summary_50_200ns.csv','1IU4_state_and_solvent_comparison.csv','1IU4_residue_pair_energy_comparison.csv']:
    pd.testing.assert_frame_equal(pd.read_csv(out/name),pd.read_csv(ROOT/'data/tables'/name),
                                  check_dtype=False,rtol=1e-9,atol=1e-10)
print(f'Reproduction passed. Three summary tables match; analysis-stage plots are in {out}.')
