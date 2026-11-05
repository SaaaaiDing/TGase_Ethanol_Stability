# Text extracted from the supplied methods document

Math objects are retained in the original DOCX; this text extraction may omit them.

Enzyme–substrate complex prediction and preparation

The amino acid sequence of microbial transglutaminase from Streptomyces mobaraensis was obtained from the Protein Data Bank entry 1IU4 [1], and the corresponding residue numbering was retained throughout the study. The enzyme–substrate complex was predicted using Boltz-2 [2]. The substrate, N-benzyloxycarbonyl-L-glutaminylglycine (Cbz-Gln-Gly; ZQG; CAS 6610-42-0), was represented by the isomeric SMILES string C1=CC=C(C=C1)COC(=O)N[C@@H](CCC(=O)N)C(=O)NCC(=O)O. The selected complex was prepared using the Protein Preparation Wizard in Maestro, Schrödinger Release 2023-4, with protonation-state assignment targeting pH 7.4 [3]. Explicit hydrogen atoms were included, and the prepared protonation states were retained during system construction. Cys64 was modeled as a neutral thiol, and the enzyme and ZQG were treated as a noncovalent complex. The same prepared enzyme–substrate model was used to construct all solvent conditions.

Construction of aqueous and ethanol–water systems

Three explicit-solvent systems were prepared: water, nominal 20% (v/v) ethanol–water, and nominal 50% (v/v) ethanol–water. The OPLS_2005 force field was used for the enzyme, substrate, and ethanol [4], together with the simple point charge (SPC) water model [5]. The aqueous system was constructed using the Desmond System Builder, whereas the mixed-solvent systems were assembled using the Schrödinger Disordered System Builder with the enzyme–ZQG complex treated as an immersed solute. Orthorhombic simulation cells and periodic boundary conditions in all three dimensions were used. The aqueous system contained 14,493 water molecules. The final 20% ethanol system contained 27,885 water and 2,112 ethanol molecules, and the final 50% ethanol system contained 19,415 water and 5,584 ethanol molecules. These compositions corresponded to ethanol solvent mole fractions of 0.07041 and 0.22337, respectively. Each system contained one enzyme molecule, one ZQG molecule, and one chloride counterion, giving total atom counts of 48,688, 107,872, and 113,710 for water, 20% ethanol, and 50% ethanol, respectively. 

For the 20% ethanol condition, the preliminary solvent box was enlarged before production simulation by adding 1,408 ethanol and 18,590 water molecules while retaining the existing molecular configuration. Steric packing during expansion used an initial density setting of 0.500 g cm⁻³, a van der Waals scale factor of 0.50, 50 placement attempts per molecule, 20 attempts per density/scale setting, and random seed 2345. Preservation of the solute geometry was verified by rigid-body alignment before and after packing. The 50% ethanol system was constructed by packing 5,584 ethanol and 19,416 water molecules around the prepared enzyme–substrate complex, followed by replacement of one water molecule with chloride. In the custom neutralization procedures, replacement-water sites were required to be at least 8 Å from solute heavy atoms and to have a minimum clearance of 2.5 Å from retained atoms. 

System relaxation and equilibration

All systems underwent relaxation before production MD, and the expanded 20% ethanol system was re-equilibrated after solvent addition. The 50% ethanol system underwent a GUI-derived, five-stage relaxation sequence comprising 100 ps of Brownian dynamics in the NVT ensemble at 10 K, 12 ps of NVT dynamics at 10 K, 12 ps of NPT dynamics at 10 K, 12 ps of NPT dynamics at 298 K, and 24 ps of unrestrained NPT dynamics at 298 K. Positional restraints during the initial stages were applied to enzyme and ZQG heavy atoms while excluding solvent molecules. This 160 ps relaxation sequence was followed by an additional 500 ps of NPT equilibration for the 50% ethanol system. 

Production molecular dynamics simulations

GPU-accelerated MD simulations were performed using Desmond within Schrödinger Release 2023-4 on the Roar Collab high-performance computing cluster at The Pennsylvania State University [6]. One 200 ns production trajectory was generated for each solvent condition, providing 600 ns of total production sampling. Simulations were performed in the NPT ensemble at 298 K and 1.01325 bar without positional restraints on the enzyme–substrate complex. Standard Desmond production settings were retained, including a Nosé–Hoover chain thermostat with a relaxation time of 1.0 ps [7] and an isotropically coupled Martyna–Tobias–Klein barostat with a relaxation time of 2.0 ps [8]. Multiple-time-step integration used time steps of 2, 2, and 6 fs for bonded, short-range nonbonded, and long-range interactions, respectively, and the short-range nonbonded cutoff was 9 Å. Coordinates were recorded at approximately 200 ps intervals, yielding 1,002 saved frames per trajectory; actual frame timestamps were used for analysis.

Trajectory processing and structural fluctuations

Trajectories were analyzed using custom Python scripts incorporating the Schrödinger Desmond topo, traj, and energygroup interfaces and NumPy. Molecular identities, atom selections, and topology-to-trajectory mappings were verified independently for each system. Periodic displacements were evaluated using the minimum-image convention, and protein coordinates were reconstructed across periodic boundaries where required. Protein Cα root-mean-square deviation (RMSD) was calculated after least-squares alignment using all 331 enzyme Cα atoms. ZQG heavy-atom RMSD was calculated after applying the protein alignment to its 24 heavy atoms, thereby retaining substrate displacement and reorientation relative to the enzyme. The minimum distance between Cys64 SG and any ZQG heavy atom was additionally monitored as a measure of general substrate proximity. Residue-level Cα root-mean-square fluctuations (RMSF) were calculated over 50–200 ns using 751 frames per condition. Coordinates were aligned through three iterations of fitting to the mean Cα configuration, and RMSF was calculated relative to each residue’s mean aligned position.

Reaction center geometry and hydrogen bonding

Reaction-center geometry was characterized by monitoring the distance between Cys64 Sγ, designated SG in the topology, and the glutamine side-chain carboxamide carbonyl carbon of ZQG, designated C28. This atom pair corresponds to the catalytic cysteine sulfur and the substrate acyl-donor carbon examined in previous modeling of microbial transglutaminase–ZQG recognition [9]. The identity of C28 was verified from its connectivity to carbonyl oxygen O22 and the primary amide nitrogen. The O22–C28–SG angle was calculated with C28 as the vertex. Frames were classified as proximal when SG–C28 ≤5 Å and distal when SG–C28 >5 Å; the fraction satisfying a 4 Å distance threshold was also calculated. Angular distributions over 50–200 ns were restricted to proximal frames and represented using normalized histograms with 5° bins. Hydrogen-bond geometry between ZQG O22 and the backbone NH groups of Cys64 or Val65 was assessed using a donor–acceptor distance ≤3.5 Å and a donor–hydrogen–acceptor angle ≥120°. Occupancies were calculated within the specified time window and proximity-state subset.

Residue contacts and local solvation

Enzyme–ZQG residue contacts were identified when any residue heavy atom was within 4.5 Å of any ZQG heavy atom. Contact occupancy was defined as the fraction of frames satisfying this criterion within the relevant analysis subset. Local solvation was quantified by counting water oxygen atoms and ethanol hydroxyl oxygen atoms within 4, 6, and 8 Å of Cys64 SG, with the 6 Å shell used for the principal comparison. The local ethanol oxygen number fraction was calculated as the summed ethanol oxygen counts divided by the summed ethanol-plus-water oxygen counts over the selected frames. This quantity was compared with the whole-system solvent number fraction, , and was distinguished from the nominal ethanol volume fraction used to designate each mixture. Analyses were summarized over 0–50 and 50–200 ns, with selected descriptors additionally evaluated in successive 50 ns windows.

Residue–substrate interaction energies

Residue–substrate interaction energies were recalculated using the Desmond energygroup interface and the corresponding simulation force field. Every fifth saved frame was analyzed, providing 201 sampled frames per trajectory at approximately 1 ns intervals. Electrostatic and van der Waals contributions were calculated separately and summed for interactions of ZQG with Tyr62, Cys64, Val65, Tyr75, Arg208, Arg238, Val252, Asn253, Phe254, Asp255, His274, Asn276, His277, and Tyr278. Energies were summarized by solvent condition, time window, and SG–C28 proximity state. The comparative energy heatmap included proximal frames; the late proximal subset for 20% ethanol was omitted because it contained only four sampled frames. Reported energies were direct nonbonded force-field interaction terms in kcal mol⁻¹.

Representative structures and descriptive analysis

Representative configurations were selected from the defined time-window/proximity subsets after protein alignment. Each candidate frame was scored using the mean squared displacement of ZQG heavy atoms from their state-specific mean coordinates plus 0.25 times the corresponding mean squared displacement of enzyme Cα atoms. The frame with the lowest score was retained, and subsets containing fewer than 20 frames were excluded from representative selection. Structures and selected atomic contacts were visualized in PyMOL. Numerical results were summarized using means, medians, quartiles, and frame occupancies as appropriate. Each solvent condition comprised one production trajectory, and comparisons were descriptive of the conformations sampled within those trajectories.

References

[1] Kashiwagi T, et al. Crystal structure of microbial transglutaminase from Streptoverticillium mobaraense. Journal of Biological Chemistry. 2002;277:44252–44260. https://doi.org/10.1074/jbc.M203933200

[2] Passaro S, et al. Boltz-2: Towards accurate and efficient binding affinity prediction. bioRxiv. 2025. https://doi.org/10.1101/2025.06.14.659707

[3] Sastry GM, Adzhigirey M, Day T, Annabhimoju R, Sherman W. Protein and ligand preparation: Parameters, protocols, and influence on virtual screening enrichments. Journal of Computer-Aided Molecular Design. 2013;27:221–234. https://doi.org/10.1007/s10822-013-9644-8

[4] Jorgensen WL, Maxwell DS, Tirado-Rives J. Development and testing of the OPLS all-atom force field on conformational energetics and properties of organic liquids. Journal of the American Chemical Society. 1996;118:11225–11236. https://doi.org/10.1021/ja9621760

[5] Berendsen HJC, Postma JPM, van Gunsteren WF, Hermans J. Interaction models for water in relation to protein hydration. In: Pullman B, editor. Intermolecular Forces. Dordrecht: D. Reidel; 1981. p. 331–342. https://doi.org/10.1007/978-94-015-7658-1_21

[6] Bowers KJ, et al. Scalable algorithms for molecular dynamics simulations on commodity clusters. Proceedings of the ACM/IEEE SC 2006 Conference. 2006. https://doi.org/10.1109/SC.2006.54

[7] Martyna GJ, Klein ML, Tuckerman M. Nosé–Hoover chains: The canonical ensemble via continuous dynamics. Journal of Chemical Physics. 1992;97:2635–2643. https://doi.org/10.1063/1.463940

[8] Martyna GJ, Tobias DJ, Klein ML. Constant pressure molecular dynamics algorithms. Journal of Chemical Physics. 1994;101:4177–4189. https://doi.org/10.1063/1.467468

[9] Tagami U, et al. Substrate specificity of microbial transglutaminase as revealed by three-dimensional docking simulation and mutagenesis. Protein Engineering, Design and Selection. 2009;22:747–752. https://doi.org/10.1093/protein/gzp061
