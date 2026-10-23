# Test fixtures

Small hand-made output fragments for program sections that have no public
reference output we could download. They copy the line layout documented in
the program manuals and tutorials; the numbers are invented and only test the
parsers, they are not results.

| file | layout taken from |
|---|---|
| `orca-neb-fragment.out` | ORCA 6.1 tutorial "Finding transition states with NEB-TS" (path summary table, `<= CI` / `<= TS` markers) |
| `orca-esd-fragment.out` | ORCA 6.0 manual, "Excited State Dynamics" (rate constant line and FC/HT line) |
| `nbo-neda-fragment.out` | water dimer NEDA summary printed in the NBO 6 tutorial, https://nbo6.chem.wisc.edu/tut_neda.htm (numbers as published there) |
| `aimall-bcp-fragment.sum` | three invented BCPs (covalent, metal–ligand, hydrogen bond) in the key = value layout that the original extraction script searched for in AIMAll .sum/.mgpviz files (`Type =`, `Rho =`, `DelSqRho =`, `Bond Ellipticity =`, `G =`, `V =`); no public AIMAll output was available |
| `simulate-cv.awk` | generator of simulated voltammograms (idealised shapes with analytically known E1/2 and onset) used by the electrochemistry tests; real potentiostat files were not available |
