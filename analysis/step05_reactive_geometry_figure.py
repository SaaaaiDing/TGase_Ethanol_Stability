#!/usr/bin/env python3
"""Summarize and plot Cys64 SG to ZQG C28 reactive-site geometry."""

import csv
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


analysis_dir = Path(os.environ["ANALYSIS_DIR"]).expanduser().resolve()
labels = ["water", "ethanol20", "ethanol50"]
titles = ["Water", "20% ethanol", "50% ethanol"]
colors = ["#333333", "#df6c14", "#1f78b4"]


def read_rows(path):
    with path.open(newline="") as handle:
        data = list(csv.DictReader(handle))
    if len(data) != 1002:
        raise RuntimeError(f"Expected 1002 frames in {path}, found {len(data)}")
    times = np.array([float(row["time_ns"]) for row in data])
    distances = np.array([float(row["Cys64_SG_to_ZQG_C28_A"]) for row in data])
    angles = np.array([float(row["O22_C28_SG_angle_deg"]) for row in data])
    if not np.all(np.isfinite(distances)) or not np.all(np.isfinite(angles)):
        raise RuntimeError(f"Nonfinite distances or angles in {path}")
    if not np.all(np.diff(times) > 0) or not 199.9 < times[-1] < 200.1:
        raise RuntimeError(f"Invalid frame times in {path}")
    return times, distances, angles


fig, axes = plt.subplots(
    3, 2, figsize=(11.4, 8.8),
    gridspec_kw={"width_ratios": [3.4, 1.2]},
    sharex="col", constrained_layout=True,
)
fig.suptitle("Cys64 approach to the ZQG glutamine carbonyl", fontsize=13)
summary = []

for row_num, (label, title, color) in enumerate(zip(labels, titles, colors)):
    times, distances, angles = read_rows(
        analysis_dir / f"{label}_reactive_geometry.csv"
    )
    after_50 = (times >= 50.0) & (times <= 200.01)
    close = distances <= 5.0
    close_after_50 = after_50 & close
    selected_angles = angles[close_after_50]
    if len(selected_angles):
        angle_median = float(np.median(selected_angles))
        angle_q25, angle_q75 = np.percentile(selected_angles, [25, 75])
    else:
        angle_median = angle_q25 = angle_q75 = float("nan")

    record = {
        "system": label,
        "frames_50_to_200_ns": int(np.sum(after_50)),
        "median_SG_C28_A_50_to_200_ns": float(np.median(distances[after_50])),
        "fraction_SG_C28_le_4A_50_to_200_ns": float(
            np.mean(distances[after_50] <= 4.0)
        ),
        "fraction_SG_C28_le_5A_50_to_200_ns": float(
            np.mean(distances[after_50] <= 5.0)
        ),
        "median_angle_deg_given_distance_le_5A_50_to_200_ns": angle_median,
        "angle_q25_deg_given_distance_le_5A_50_to_200_ns": float(angle_q25),
        "angle_q75_deg_given_distance_le_5A_50_to_200_ns": float(angle_q75),
    }
    summary.append(record)

    ax = axes[row_num, 0]
    ax.plot(times, distances, color=color, linewidth=0.85)
    ax.axhline(5.0, color="#777777", linestyle="--", linewidth=0.9)
    ax.axvline(50.0, color="#aaaaaa", linestyle=":", linewidth=0.9)
    ax.set_xlim(0, 200)
    ax.set_ylim(0, max(10.0, np.ceil(np.max(distances) / 5.0) * 5.0))
    ax.set_ylabel(f"{title}\nSG–C28 (Å)")
    ax.text(
        0.98, 0.94,
        f"50–200 ns: {100 * record['fraction_SG_C28_le_5A_50_to_200_ns']:.1f}% ≤5 Å",
        transform=ax.transAxes, ha="right", va="top", fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82},
    )

    angle_ax = axes[row_num, 1]
    if len(selected_angles):
        angle_ax.hist(
            selected_angles, bins=np.arange(0, 185, 5),
            color=color, alpha=0.85, density=True,
        )
        angle_ax.axvline(angle_median, color="black", linewidth=0.9)
    angle_ax.set_xlim(0, 180)
    angle_ax.set_ylabel("Density")
    angle_ax.set_title(title, fontsize=10)

axes[2, 0].set_xlabel("Production time (ns)")
axes[2, 1].set_xlabel("O22–C28–SG angle (degrees)\n50–200 ns, SG–C28 ≤5 Å")

for extension in ("svg", "png"):
    output = analysis_dir / f"Cys64_ZQG_C28_geometry.{extension}"
    fig.savefig(output, dpi=300)
    print(f"Saved: {output}")
plt.close(fig)

summary_path = analysis_dir / "reactive_geometry_summary_50_200ns.csv"
with summary_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
    writer.writeheader()
    writer.writerows(summary)

for row in summary:
    print(
        f"{row['system']}: median SG-C28 "
        f"{row['median_SG_C28_A_50_to_200_ns']:.2f} A; "
        f"<=4 A {100 * row['fraction_SG_C28_le_4A_50_to_200_ns']:.1f}%; "
        f"<=5 A {100 * row['fraction_SG_C28_le_5A_50_to_200_ns']:.1f}%; "
        f"angle median [Q1,Q3] conditional on <=5 A: "
        f"{row['median_angle_deg_given_distance_le_5A_50_to_200_ns']:.1f} "
        f"[{row['angle_q25_deg_given_distance_le_5A_50_to_200_ns']:.1f}, "
        f"{row['angle_q75_deg_given_distance_le_5A_50_to_200_ns']:.1f}] degrees"
    )
print(f"Saved: {summary_path}")
