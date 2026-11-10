# Score function weights for the hybridize stages

RosettaCM scores its three stages with three different functions: a coarse
centroid function for the large recombination moves, a second centroid function
for refinement, and a full atom function for the final relax. These six files
are those functions as used for the modelling in this module.

| file | stage | representation |
|---|---|---|
| `stage1.wts` | recombination | centroid, soluble |
| `stage2.wts` | centroid refinement | centroid, soluble |
| `stage3.wts` | full atom relax | full atom, soluble |
| `stage1_membrane.wts` | recombination | centroid, membrane |
| `stage2_membrane.wts` | centroid refinement | centroid, membrane |
| `stage3_rlx_membrane.wts` | full atom relax | full atom, membrane |

`rosetta-cm-hybridize.sh` picks the soluble set by default and the membrane set
with `-M`.

## What the membrane variants change

They add the membrane environment terms, and they reweight the generic ones to
suit a protein whose surface is lipid rather than water:

- `Menv`, `Menv_non_helix`, `Menv_termini`, `Menv_tm_proj` score residues
  against the predicted membrane environment, which is why these runs need the
  `.span` file;
- `Mpair` and `Mcbeta` replace the soluble `pair` and `cbeta` pair terms;
- `Mlipo` scores lipophilicity against depth in the bilayer;
- `vdw` rises from 0.2 to 3.0 and `rg` falls from 2.0 to 0.1. The radius of
  gyration term pulls a soluble protein into a globule, which is wrong for a
  transmembrane helix bundle, so it is nearly switched off and the clash term
  carries more of the work instead.

That `rg` change is the one to understand before using these. If you model a
membrane protein with the soluble weights, `rg` will compact the helices and
the result will look like a plausible globular protein and be wrong.

## Provenance

The soluble weights are the RosettaCM stage weights as published with the
method (Song, DiMaio, Wang, Kim, Miles, Brunette, Thompson and Baker,
*Structure* **2013**, 21, 1735). Rosetta ships equivalents in its demos. The
membrane variants are the adaptation used for the work in this module, built on
the Rosetta membrane framework of Yarov-Yarovoy, Schonbrun and Baker
(*Proteins* **2006**, 62, 1010).

They are included because they are the parameters the reported models were
produced with, and a protocol is not reproducible without them. They are score
term names and numbers rather than program code. If you would rather take them
from your own Rosetta installation, the soluble set is in Rosetta's demos and
`-w` points the hybridize script at any directory holding the six names above.
