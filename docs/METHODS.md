# Methods and interpretation

This description consolidates the supplied methods document, analysis scripts and production configurations. The original [Word document](TG_ZQG_Methods_Han.docx), its [text extraction](archive/methods_document_text.md), the [earlier analysis notes](archive/analysis_methods_and_legends.md) and the [final figure legends](FIGURE_LEGENDS.md) are preserved. The Word document contains the study's reference list. No new MD was run to create this repository.

## Model and preparation

The supplied methods identify the enzyme as microbial transglutaminase from *Streptomyces mobaraensis*, using PDB 1IU4 sequence/residue numbering, and the complex prediction method as Boltz-2. ZQG is N-benzyloxycarbonyl-L-glutaminylglycine (Cbz-Gln-Gly), represented upstream as `C1=CC=C(C=C1)COC(=O)N[C@@H](CCC(=O)N)C(=O)NCC(=O)O`.

Protein Preparation Wizard in Schrödinger 2023-4 assigned protonation states targeting pH 7.4. Cys64 was modeled as a neutral thiol, and the enzyme–ZQG complex remained noncovalent. There are 331 enzyme Cα atoms and 24 ZQG heavy atoms (42 ZQG atoms including hydrogens). These counts and the identity/connectivity of the reactive atoms are checked in the analysis scripts. Source atom IDs are topology-specific and must not be transferred between conditions without those checks.

## Solvent systems

OPLS_2005 was used with SPC water and orthorhombic periodic boxes. System composition is listed in the main README and checked in the preparation scripts. Volume-fraction labels are nominal preparation labels. The ethanol molecule-number fractions are 0.070407 and 0.223369 for the 20% and 50% mixtures, respectively.

The 20% condition was expanded before production by adding 1,408 ethanol and 18,590 water molecules while preserving the existing solute configuration. Recorded expansion settings were 0.500 g/cm³ initial packing density, vdW scale 0.50, 50 placement attempts per molecule, 20 attempts per density/scale setting and packing seed 2345. This initial packing density is not the equilibrated liquid density.

The 50% condition was packed with 5,584 ethanol and 19,416 waters. One water was replaced by chloride, giving 19,415 waters in the neutralized system. Replacement sites were required to be at least 8 Å from solute heavy atoms with at least 2.5 Å clearance from retained atoms. `prepare_1IU4_EtOH50.py` depends on the included `expand_1IU4_EtOH20.py` and `workflow_1IU4_EtOH20.py` helpers.

## Relaxation and production

The supplied 50% workflow uses five GUI-derived relaxation stages: 100 ps Brownian NVT at 10 K, 12 ps NVT at 10 K, 12 ps NPT at 10 K, 12 ps NPT at 298 K and 24 ps unrestrained NPT at 298 K. Early restraints select enzyme/ZQG heavy atoms and exclude solvent. A further 500 ps NPT equilibration precedes production. Its submission script checks composition, force field, density, volume drift and periodic-image separation. The expanded 20% system was also re-equilibrated before production; the actual relaxation MSJ is included.

| Parameter | Production setting |
|---|---|
| Engine | Desmond, Schrödinger Release 2023-4 |
| Force field / water | OPLS_2005 / SPC |
| Ensemble | NPT, unrestrained enzyme–substrate complex |
| Temperature / pressure | 298 K / 1.01325 bar, isotropic |
| Thermostat | Nosé–Hoover chain, 1.0 ps time constant |
| Barostat | Martyna–Tobias–Klein, 2.0 ps time constant |
| Multiple time steps | 0.002 / 0.002 / 0.006 ps |
| Short-range cutoff | 9 Å |
| Production time | 200,000 ps per condition |
| Recorded production seed | 2007 |
| Coordinate sampling | Approximately 200 ps; 1,002 saved frames |

All three engine logs contain a normal exit. Their last periodic `Chemical time` message is 199,999.2 ps; this logging cadence is distinct from the actual terminal saved frame near 200 ns. The CSV endpoint and frame-count checks are recorded separately. Original final configurations retain their `ORIG_CFG` section for Desmond energy-group analysis.

## Analysis definitions

| Quantity | Calculation | Interpretation |
|---|---|---|
| Protein Cα RMSD | Least-squares fit of 331 Cα atoms to the first saved frame | Global deviation within each trajectory |
| ZQG RMSD | RMSD of 24 ligand heavy atoms after the protein fit | Ligand movement relative to the enzyme |
| General active-site proximity | Minimum distance from Cys64 SG to any ligand heavy atom | Distinct from the reactive-carbon distance |
| Cα RMSF | 50–200 ns; three iterative fits to the mean configuration, then per-residue fluctuations | Local flexibility within a single trajectory |
| Reactive-carbon proximity | Minimum-image Cys64 SG–ZQG C28 distance; near ≤5 Å, far >5 Å; 4 Å also tabulated | A proximity state, not a catalytic-state assignment |
| Carbonyl approach angle | O22–C28–SG, vertex at C28; 5° histogram bins among near frames | Geometry of the noncovalent complex |
| Backbone H-bond proxy | Cys64 or Val65 N···O22 ≤3.5 Å and N–H···O22 ≥120° | Geometry-based contact, conditional on the selected frames |
| Residue contact | Any residue/ligand heavy-atom pair within 4.5 Å | Occupancy within the chosen state/window |
| Local solvent shell | Count water O and ethanol hydroxyl O within 4, 6 and 8 Å of Cys64 SG | Six-angstrom shell is the principal comparison |
| Local ethanol number fraction | Sum ethanol O counts / sum (ethanol O + water O counts) across selected frames | A local composition descriptor; not a free energy or preferential-interaction coefficient |
| Residue–ZQG energy | Desmond energygroup electrostatic + vdW pair terms | Direct force-field interaction energy; not MM/GBSA, binding ΔG or an activation barrier |

Periodic displacements use minimum-image distances. Protein coordinates are reconstructed/aligned where the relevant analysis requires it. The energy script verifies atom-group identities, nonoverlap, frame-to-reference time mapping and finite energy components.

## Sampling and comparisons

The dense analysis uses 1,002 frames over the full trajectory and 751 frames in the late window. Energy evaluation takes indices `0, 5, ..., 1000`: 201 frames per condition, approximately 1 ns apart. Actual saved timestamps determine window membership: the sampled early and late windows contain 51 and 150 frames. Nominal boundary labels must not be used to replace the actual timestamps.

The comparative energy heatmap includes 14 residues: Tyr62, Cys64, Val65, Tyr75, Arg208, Arg238, Val252, Asn253, Phe254, Asp255, His274, Asn276, His277 and Tyr278. It compares near frames in early water/20%/50% and late water/50% subsets. The late 20% near subset has only four sampled frames and is omitted from that heatmap. Small subsets are flagged in the tables.

Representative structures minimize the mean squared displacement of ligand heavy atoms from the state mean plus 0.25 times the corresponding protein Cα term. Subsets with fewer than 20 dense frames are excluded. The principal displayed frames are approximately 116.0 ns (water), 125.4 ns (20%) and 141.8 ns (50%). Exact timestamps and zero-based indices are preserved in the representative-frame tables. A selected frame is not an occupancy estimate.

## Scope of conclusions

There is one trajectory per solvent condition. Adjacent frames are correlated; frame SDs and counts are descriptive and do not quantify uncertainty across independent simulations. No significance test or replicate-based confidence interval is claimed. The neutral-cysteine classical force field does not model covalent catalysis or changes in protonation during a reaction. Favorable charged-residue pair energies can occur without direct short-range contacts and do not establish stronger binding. Neither the local ethanol shell nor the high 50%-ethanol near fraction establishes increased enzyme activity.
