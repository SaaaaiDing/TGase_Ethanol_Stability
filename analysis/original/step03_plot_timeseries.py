from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(
    "/storage/home/sqd5856/Desktop/default/HanLiu/prep_files/yaml/1IU4/"
    "boltz_results_1IU4/predictions/1IU4/md_analysis_0_20_50"
)

CONDITIONS = {
    "water": ("Water", "#333333"),
    "ethanol20": ("20% ethanol", "#D55E00"),
    "ethanol50": ("50% ethanol", "#0072B2"),
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "svg.fonttype": "none",
    "axes.linewidth": 0.8,
})

fig, axes = plt.subplots(
    3, 1, figsize=(8.0, 7.2), sharex=True,
    constrained_layout=True
)

for key, (label, color) in CONDITIONS.items():
    data = np.loadtxt(
        ROOT / f"{key}_active_site_proximity.csv",
        delimiter=",",
        skiprows=1,
    )

    time = data[:, 0]
    cys_distance = data[:, 1]
    protein_rmsd = data[:, 2]
    ligand_rmsd = data[:, 3]

    axes[0].plot(time, protein_rmsd, color=color, lw=0.8, label=label)
    axes[1].plot(time, ligand_rmsd, color=color, lw=0.8)
    axes[2].plot(time, cys_distance, color=color, lw=0.8)

    print(f"\n{label}")
    for start in (0, 50, 100, 150):
        selected = (time >= start) & (time < start + 50)
        print(
            f"  {start:3d}-{start + 50:3d} ns: "
            f"median ligand RMSD={np.median(ligand_rmsd[selected]):6.2f} A; "
            f"median Cys64 distance={np.median(cys_distance[selected]):6.2f} A; "
            f"distance <=8 A={np.mean(cys_distance[selected] <= 8):6.1%}"
        )

axes[0].set_ylabel("Protein CA RMSD (A)")
axes[1].set_ylabel("ZQG heavy-atom RMSD (A)")
axes[2].set_ylabel("Cys64 SG to ZQG\nnearest heavy atom (A)")
axes[2].set_xlabel("Production time (ns)")

for ax in axes:
    ax.set_xlim(0, 200)
    ax.set_ylim(bottom=0)
    ax.axvline(50, color="#888888", ls="--", lw=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes[0].legend(frameon=False, ncol=3, loc="upper left")
axes[2].set_xticks(np.arange(0, 201, 25))

svg_path = ROOT / "protein_ZQG_active_site_timeseries.svg"
png_path = ROOT / "protein_ZQG_active_site_timeseries.png"
fig.savefig(svg_path)
fig.savefig(png_path, dpi=600)
plt.close(fig)

print(f"\nSaved: {svg_path}")
print(f"Saved: {png_path}")
