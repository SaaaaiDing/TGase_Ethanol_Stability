#!/usr/bin/env python3
"""Make descriptive figures and tables from the three Desmond pair-energy CSVs.

Usage: python make_1IU4_md_comparison.py --input-dir PATH --output-dir PATH
This script uses sampled frames only; it does not estimate independent replicates.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


SYSTEMS = ("water", "ethanol20", "ethanol50")
LABELS = {"water": "Water", "ethanol20": "20% ethanol", "ethanol50": "50% ethanol"}
COLORS = {"water": "#0072B2", "ethanol20": "#D55E00", "ethanol50": "#009E73"}
ETHANOL_COUNT = {"water": 0, "ethanol20": 2112, "ethanol50": 5584}
WATER_COUNT = {"water": 14493, "ethanol20": 27885, "ethanol50": 19415}
RESIDUES = (
    "CYS64", "VAL65", "HIS274", "ASP255", "TYR62", "TYR75",
    "VAL252", "ASN253", "PHE254", "ASN276", "HIS277", "TYR278",
    "ARG208", "ARG238",
)
PERIODS = ("0-50", "50-200")
STATES = ("near", "far")
RAW_COLS = (
    "frame_index", "time_ns", "period_ns", "SG_C28_state", "SG_C28_A",
    "Cys64_NH_O22_proxy", "Val65_NH_O22_proxy", "residue", "partner",
    "electrostatic_kcal_per_mol", "vdw_kcal_per_mol",
    "nonbonded_sum_kcal_per_mol",
)


def load_inputs(directory):
    energies, frames, solvent = {}, {}, {}
    for system in SYSTEMS:
        raw = pd.read_csv(directory / f"{system}_pair_energies_stride5.csv")
        missing = set(RAW_COLS) - set(raw.columns)
        if missing or len(raw) != 201 * 14 * (2 if system == "water" else 3):
            raise ValueError(f"{system}: missing columns {missing} or wrong row count")
        if not np.allclose(
            raw.electrostatic_kcal_per_mol + raw.vdw_kcal_per_mol,
            raw.nonbonded_sum_kcal_per_mol, atol=1e-6,
        ):
            raise ValueError(f"{system}: elec + vdw does not match total")
        q = raw.loc[(raw.residue == "CYS64") & (raw.partner == "ZQG")].copy()
        q = q.sort_values("frame_index").reset_index(drop=True)
        if len(q) != 201 or not np.array_equal(q.frame_index, np.arange(0, 1001, 5)):
            raise ValueError(f"{system}: unexpected frame indices")
        if not np.array_equal(
            raw.groupby("frame_index").size().to_numpy(),
            np.full(201, 28 if system == "water" else 42),
        ):
            raise ValueError(f"{system}: duplicated or missing pairs")
        sh = pd.read_csv(directory / f"{system}_Cys64_solvent_shell.csv")
        if len(sh) != 1002 or not np.allclose(
            sh.iloc[q.frame_index].time_ns.to_numpy(), q.time_ns.to_numpy(), atol=1e-5
        ):
            raise ValueError(f"{system}: solvent shell frame mapping failed")
        energies[system] = raw
        frames[system] = q
        solvent[system] = sh.iloc[q.frame_index].reset_index(drop=True)
    return energies, frames, solvent


def build_tables(energies, frames, solvent):
    state_rows, pair_rows = [], []
    for system in SYSTEMS:
        frame = frames[system]
        shell = solvent[system]
        bulk = ETHANOL_COUNT[system] / (ETHANOL_COUNT[system] + WATER_COUNT[system])
        for period in PERIODS:
            mask_period = frame.period_ns.eq(period)
            n_period = int(mask_period.sum())
            for state in STATES:
                mask = mask_period & frame.SG_C28_state.eq(state)
                n = int(mask.sum())
                if n == 0:
                    continue
                water6 = shell.loc[mask, "water_O_6A"]
                ethanol6 = shell.loc[mask, "ethanol_O_6A"]
                oxygen_sum = int(water6.sum() + ethanol6.sum())
                hbond_count = int(frame.loc[mask, "Cys64_NH_O22_proxy"].sum())
                state_rows.append({
                    "system": system, "period_ns": period, "SG_C28_state": state,
                    "n_sampled_frames": n, "n_period_frames": n_period,
                    "state_fraction": n / n_period,
                    "median_SG_C28_A": frame.loc[mask, "SG_C28_A"].median(),
                    "Cys64_NH_O22_count": hbond_count,
                    "Cys64_NH_O22_fraction": hbond_count / n,
                    "mean_water_O_6A": water6.mean(),
                    "mean_ethanol_O_6A": ethanol6.mean(),
                    "local_ethanol_O_fraction_6A":
                        ethanol6.sum() / oxygen_sum if oxygen_sum else np.nan,
                    "bulk_ethanol_O_number_fraction": bulk,
                    "small_sample_n_lt_10": n < 10,
                })
                selected = energies[system].loc[
                    energies[system].period_ns.eq(period)
                    & energies[system].SG_C28_state.eq(state)
                    & energies[system].partner.eq("ZQG")
                ]
                for residue in RESIDUES:
                    rows = selected.loc[selected.residue.eq(residue)]
                    if len(rows) != n:
                        raise ValueError(f"{system}: missing {residue} in {period} {state}")
                    pair_rows.append({
                        "system": system, "period_ns": period,
                        "SG_C28_state": state, "residue": residue,
                        "n_sampled_frames": n,
                        "mean_electrostatic_kcal_per_mol":
                            rows.electrostatic_kcal_per_mol.mean(),
                        "mean_vdw_kcal_per_mol": rows.vdw_kcal_per_mol.mean(),
                        "mean_pair_kcal_per_mol": rows.nonbonded_sum_kcal_per_mol.mean(),
                        "median_pair_kcal_per_mol": rows.nonbonded_sum_kcal_per_mol.median(),
                        "frame_sd_pair_kcal_per_mol":
                            rows.nonbonded_sum_kcal_per_mol.std(ddof=1),
                        "small_sample_n_lt_10": n < 10,
                    })
    return pd.DataFrame(state_rows), pd.DataFrame(pair_rows)


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelsize": 9, "axes.titlesize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "figure.dpi": 120, "savefig.dpi": 400,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
    })


def label_panel(ax, letter):
    ax.text(-0.12, 1.05, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", ha="left", va="bottom")


def make_state_figure(frames, state_table, output):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.2), layout="constrained")
    ax = axes[0, 0]
    for system in SYSTEMS:
        f = frames[system]
        color = COLORS[system]
        ax.scatter(f.time_ns, f.SG_C28_A, s=8, alpha=.32, color=color,
                   linewidths=0, rasterized=True)
        smooth = f.SG_C28_A.rolling(9, center=True, min_periods=1).median()
        ax.plot(f.time_ns, smooth, color=color, lw=1.65, label=LABELS[system])
    ax.axhline(5, color="#555555", ls="--", lw=1)
    ax.axvline(50, color="#B6B6B6", ls=":", lw=1)
    ax.set(xlim=(0, 200), ylim=(0, 45), xlabel="Trajectory time (ns)",
           ylabel="Cys64 SG–ZQG C28 (Å)",
           title="Reactive-carbon distance (9-frame median)")
    ax.legend(frameon=False, loc="upper left", ncol=1, fontsize=8)
    label_panel(ax, "A")

    ax = axes[0, 1]
    periods = np.arange(2)
    offsets = (-.26, 0, .26)
    for system, offset in zip(SYSTEMS, offsets):
        rows = state_table.loc[
            state_table.system.eq(system) & state_table.SG_C28_state.eq("near")
        ].set_index("period_ns")
        y = np.array([rows.loc[p, "state_fraction"] * 100 for p in PERIODS])
        ax.bar(periods + offset, y, width=.23, color=COLORS[system],
               label=LABELS[system])
        for x, val, p in zip(periods + offset, y, PERIODS):
            n = int(rows.loc[p, "n_sampled_frames"])
            ax.text(x, val + 1.6, f"{n}", ha="center", va="bottom", fontsize=7)
    ax.set(xticks=periods, xticklabels=("0–50", "50–200"), ylim=(0, 109),
           xlabel="Trajectory window (ns)", ylabel="Reactive-near frames (%)",
           title="SG–C28 ≤5 Å (counts above bars)")
    label_panel(ax, "B")

    ax = axes[1, 0]
    for system, offset in zip(SYSTEMS, offsets):
        rows = state_table.loc[
            state_table.system.eq(system) & state_table.SG_C28_state.eq("near")
        ].set_index("period_ns")
        for j, p in enumerate(PERIODS):
            row = rows.loc[p]
            x, n = j + offset, int(row.n_sampled_frames)
            if n < 10:
                ax.plot(x, 0, marker="x", ms=8, color=COLORS[system])
                ax.text(x, 5, f"n={n}", ha="center", fontsize=7,
                        color=COLORS[system])
            else:
                y = row.Cys64_NH_O22_fraction * 100
                ax.bar(x, y, width=.23, color=COLORS[system])
                ax.text(x, y + 1.6, f"{y:.0f}%", ha="center", fontsize=7)
    ax.set(xticks=periods, xticklabels=("0–50", "50–200"), ylim=(0, 107),
           xlabel="Trajectory window (ns)",
           ylabel="Cys64 NH–O22 geometry (%)",
           title="Conditional on reactive-near frames")
    label_panel(ax, "C")

    ax = axes[1, 1]
    entries = (
        ("ethanol20", "0-50", "near", "20%\nearly near"),
        ("ethanol20", "50-200", "far", "20%\nlate far"),
        ("ethanol50", "0-50", "near", "50%\nearly near"),
        ("ethanol50", "50-200", "near", "50%\nlate near"),
    )
    for i, (system, period, state, label) in enumerate(entries):
        row = state_table.loc[
            state_table.system.eq(system) & state_table.period_ns.eq(period)
            & state_table.SG_C28_state.eq(state)
        ].iloc[0]
        local = row.local_ethanol_O_fraction_6A * 100
        bulk = row.bulk_ethanol_O_number_fraction * 100
        ax.plot([i, i], [local, bulk], color=COLORS[system], lw=1.4, alpha=.7)
        ax.scatter(i, local, color=COLORS[system], s=50, marker="o", zorder=3)
        ax.scatter(i, bulk, color=COLORS[system], s=55, marker="D",
                   facecolors="none", lw=1.5, zorder=3)
        ax.text(i, local + 1.1, f"{local:.1f}%", ha="center", fontsize=7)
    ax.plot([], [], "ko", label="Cys64 6 Å shell")
    ax.plot([], [], marker="D", ls="", markerfacecolor="none", color="k",
            label="Whole-system O number fraction")
    ax.set(xticks=np.arange(4), xticklabels=[x[3] for x in entries],
           ylim=(0, 29), ylabel="Ethanol O / (ethanol O + water O) (%)",
           title="Local solvent composition at Cys64")
    ax.legend(frameon=False, loc="upper left", fontsize=7)
    label_panel(ax, "D")

    fig.suptitle("1IU4–ZQG: reactive geometry and Cys64 solvation",
                 fontsize=13, fontweight="bold")
    for suffix in ("svg", "png"):
        fig.savefig(output / f"1IU4_reactive_state_and_solvent.{suffix}",
                    bbox_inches="tight")
    plt.close(fig)


def make_energy_figure(pair_table, output):
    columns = (
        ("water", "0-50", "Water\n0–50 ns", 38),
        ("ethanol20", "0-50", "20% EtOH\n0–50 ns", 45),
        ("ethanol50", "0-50", "50% EtOH\n0–50 ns", 49),
        ("water", "50-200", "Water\n50–200 ns", 73),
        ("ethanol50", "50-200", "50% EtOH\n50–200 ns", 136),
    )
    data = np.empty((len(RESIDUES), len(columns)), dtype=float)
    for i, residue in enumerate(RESIDUES):
        for j, (system, period, _, expected_n) in enumerate(columns):
            r = pair_table.loc[
                pair_table.system.eq(system) & pair_table.period_ns.eq(period)
                & pair_table.SG_C28_state.eq("near")
                & pair_table.residue.eq(residue)
            ].iloc[0]
            if int(r.n_sampled_frames) != expected_n:
                raise ValueError("Unexpected figure sample count")
            data[i, j] = r.mean_pair_kcal_per_mol
    fig, ax = plt.subplots(figsize=(8.5, 7.7), layout="constrained")
    image = ax.imshow(data, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-20, vcenter=0, vmax=20),
                      interpolation="none", aspect="auto")
    ax.set_xticks(np.arange(len(columns)),
                  [f"{c[2]}\n(n={c[3]})" for c in columns])
    ax.set_yticks(np.arange(len(RESIDUES)), RESIDUES)
    ax.tick_params(length=0, axis="both", pad=7)
    ax.axvline(2.5, color="white", linewidth=7)
    ax.set_title("Residue–ZQG nonbonded pair energy | SG–C28 ≤5 Å",
                 fontsize=12, fontweight="bold", pad=16)
    for i in range(len(RESIDUES)):
        for j in range(len(columns)):
            color = "white" if abs(data[i, j]) > 10 else "#202020"
            ax.text(j, i, f"{data[i, j]:.1f}", ha="center", va="center",
                    color=color, fontsize=8)
    cb = fig.colorbar(image, ax=ax, fraction=.042, pad=.035)
    cb.set_label("Mean electrostatic + vdW energy (kcal/mol)", fontsize=9)
    ax.text(0, -0.10, "20% ethanol late-near: n=4; excluded from this comparison.",
            transform=ax.transAxes, fontsize=8, color="#444444", va="top")
    ax.text(0, -0.15, "Descriptive interaction terms; neither binding free energies nor residue ΔG.",
            transform=ax.transAxes, fontsize=8, color="#444444", va="top")
    for suffix in ("svg", "png"):
        fig.savefig(output / f"1IU4_residue_ZQG_pair_energy.{suffix}",
                    bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    energies, frames, solvent = load_inputs(args.input_dir)
    state_table, pair_table = build_tables(energies, frames, solvent)
    state_table.to_csv(args.output_dir / "1IU4_state_and_solvent_comparison.csv", index=False)
    pair_table.to_csv(args.output_dir / "1IU4_residue_pair_energy_comparison.csv", index=False)
    style()
    make_state_figure(frames, state_table, args.output_dir)
    make_energy_figure(pair_table, args.output_dir)
    print(state_table.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"Saved figures and tables to {args.output_dir}")


if __name__ == "__main__":
    main()
