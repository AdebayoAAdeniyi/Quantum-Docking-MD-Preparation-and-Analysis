# Examples: real runs, with their input and output

Every folder here is the output of an actual run made while building this
repository, not an illustration. The programs used are the free ones, so these
examples can be reproduced on any machine.

| folder | what it is |
|---|---|
| `molecules/` | five 2-hydroxybenzophenones built from SMILES (xyz, pdb, mol2) and the SMILES table they came from; used as the test set throughout |
| `xtb-benzophenones/` | GFN2-xTB screening of those five in implicit water: energies, gaps, dipoles, IP, EA, chemical potential, hardness, electrophilicity, and the condensed Fukui indices of every atom |
| `kinetics-hcn-hnc/` | NWChem B3LYP/6-31G\* study of HCN → HNC: the three inputs, the extracted thermochemistry at four temperatures, the barriers, the Eyring rates with Wigner tunnelling, the Arrhenius fit and the energy profile figure |
| `amber-imatinib-water/` | the Amber protocol run end to end on imatinib in water: leap and MD inputs, the energy table and the equilibration figure (temperature ramp, total energy, density) |
| `docking-1iep/` | AutoDock Vina on Abl kinase with imatinib (PDB 1IEP): the box, all poses, the best per ligand with the redocking RMSD, and a consensus ranking example |
| `lammps-bead-spring/` | a coarse-grained polymer melt: the cooling curve and its glass transition fit, the stress-strain curve and modulus, the Green-Kubo viscosity, RDF and MSD |
| `metad-double-well/` | well-tempered metadynamics on a double well whose barrier is known exactly (5 units); the recovered free energy profile and its convergence series |

The numbers in the module READMEs come from these files.
