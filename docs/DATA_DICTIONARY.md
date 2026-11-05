# Data dictionary

Distances are in Å, times in ns, angles in degrees and energies in kcal/mol unless the column name specifies otherwise. State occupancies are fractions from 0 to 1, not percentages. Window/state conditioning must be matched before comparing rows.

| File pattern | Content |
|---|---|
| `*_rmsd.csv` | Time, enzyme Cα RMSD, protein-aligned ZQG heavy-atom RMSD |
| `*_active_site_proximity.csv` | Minimum Cys64 SG–ligand-heavy-atom separation plus the RMSD time series |
| `*_reactive_geometry.csv` | SG–C28 distance and O22–C28–SG angle for all 1,002 saved frames |
| `reactive_geometry_summary_50_200ns.csv` | Dense late-window proximity fractions and conditional angle quantiles |
| `protein_CA_RMSF_50_200ns.csv` | Per-residue Cα fluctuations for the three conditions |
| `*_carbonyl_backbone_hbonds.csv` | Per-frame Cys64/Val65 backbone NH–O22 geometry and binary proxies |
| `*_Cys64_solvent_shell.csv` | Per-frame local water/ethanol oxygen counts at several shell radii |
| `Cys64_solvent_shell_summary.csv`, `Cys64_solvent_by_ZQG_state.csv` | Window-level and state-conditioned solvent descriptors |
| `ZQG_residue_contacts_by_state.csv` | Residue contact occupancies by time window and proximity state |
| `*_pair_energies_stride5.csv` | Explicit-solvent pair terms for 201 sampled frames; multiple residue/partner rows per frame |
| `*_ligand_pair_summary_stride5.csv` | Per-residue ZQG pair-energy summaries by window and state |
| `1IU4_state_and_solvent_comparison.csv` | Principal sampled-frame state, hydrogen-bond and local-solvent summary |
| `1IU4_residue_pair_energy_comparison.csv` | Principal state-conditioned residue–ZQG energy summary |
| `representative_frames_portable.csv` | Selected state representatives, exact time/index and portable snapshot path |
| `ethanol50_energy_pilot_50_55ns.csv` | Historical pilot; not the final full-trajectory comparison |

## Principal comparison columns

| Column | Definition |
|---|---|
| `system` | `water`, `ethanol20` or `ethanol50` |
| `period_ns` | Timestamp-based early/late window, labeled `0-50` or `50-200` |
| `SG_C28_state` | `near` for SG–C28 ≤5 Å; otherwise `far` |
| `n_sampled_frames` | Number of sampled frames in this window/state |
| `n_period_frames` | All sampled frames in the window, before state conditioning |
| `state_fraction` | `n_sampled_frames / n_period_frames` |
| `Cys64_NH_O22_count` | Frames meeting the backbone hydrogen-bond proxy within the selected state |
| `Cys64_NH_O22_fraction` | Proxy count divided by `n_sampled_frames` |
| `mean_water_O_6A`, `mean_ethanol_O_6A` | Mean oxygen counts within 6 Å of Cys64 SG |
| `local_ethanol_O_fraction_6A` | Sum ethanol O counts / sum total solvent O counts; not the mean of per-frame ratios |
| `bulk_ethanol_O_number_fraction` | Ethanol molecules / (ethanol + water molecules), not volume fraction |
| `small_sample_n_lt_10` | Flag for a subset containing fewer than ten sampled frames |
| `mean_electrostatic_kcal_per_mol` | Mean residue–ZQG electrostatic pair term |
| `mean_vdw_kcal_per_mol` | Mean residue–ZQG van der Waals pair term |
| `mean_pair_kcal_per_mol` | Mean sum of the two force-field pair terms |
| `median_pair_kcal_per_mol` | Median of the summed pair terms |
| `frame_sd_pair_kcal_per_mol` | SD across correlated sampled frames; not a replicate uncertainty |

Raw pair-energy files also contain residue–water and residue–ethanol group interactions. Those groups contain every corresponding solvent molecule, and the system sizes differ. Their raw totals should not be treated as directly comparable binding-affinity measures.
