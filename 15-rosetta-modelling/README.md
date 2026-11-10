# 15 · Protein structure modelling with Rosetta

Building a protein structure when there is no structure to start from: a
comparative model from several partial templates, assembly into an oligomer,
and refinement under experimental distance restraints.

This is where the systems that modules 08 to 12 simulate come from. The HCV E1E2
glycoprotein work behind this module had no experimental structure of the
complex, so the model had to be built, oligomerised, glycosylated and embedded
in a membrane before any molecular dynamics could begin.

## The pipeline

```
sequence + templates + topology prediction
        │
        ├─ 1. rosetta-cm-setup.sh     alignments to Grishin, topology to span
        ├─ 2. rosetta-cm-thread.sh    partial_thread each template
        ├─ 3. rosetta-cm-hybridize.sh recombine the threaded models
        └─ 4. rosetta-relax.sh        relax, then rank
                │
                ▼ one monomer
        ├─ rosetta-symm-dock.sh       Cn symmetry, for a homo-oligomer
        └─ rosetta-local-dock.sh      interface by interface, for a hetero-assembly
                │
                ▼
          rosetta-constrained-relax.sh  impose measured contacts
                │
                ▼ glycosylation, membrane, then modules 08 to 12
```

```
15-rosetta-modelling/
├── prepare/
│   ├── rosetta-cm-setup.sh          stage 1: Grishin alignments, span file
│   ├── rosetta-cm-thread.sh         stage 2: threading
│   └── rosetta-chains.sh            chain IDs and per-chain renumbering
├── run/
│   ├── rosetta-cm-hybridize.sh      stage 3: RosettaScripts Hybridize
│   ├── rosetta-relax.sh             stage 4: dual-space relax and ranking
│   ├── rosetta-symm-dock.sh         cyclic symmetric assembly
│   ├── rosetta-local-dock.sh        local refine and dock, per interface
│   └── rosetta-constrained-relax.sh relax under distance restraints
├── extract/
│   └── rosetta-scores.sh            score files to ranked CSV
├── plot/
│   └── plot-score-distribution.sh   score histogram and funnel
└── wts/                             the three-stage score weights, soluble and membrane
```

## Workflow

A membrane protein modelled from three templates, then built into a homotrimer
and restrained to measured contacts:

```bash
# 1. alignments and membrane topology
15-rosetta-modelling/prepare/rosetta-cm-setup.sh -t e2 -a e2_templates.aln \
    -n E2_target -O e2.octopus

# 2. thread each template, naming the chain to use
15-rosetta-modelling/prepare/rosetta-cm-thread.sh -t e2 -n E2_target \
    6urh:C 2kzq:A 4mwf:C

# 3. recombine them
15-rosetta-modelling/run/rosetta-cm-hybridize.sh -t e2 -n 10000 -M -q
bin/qdm-submit 3_hybridize

# 4. relax and rank, enforcing the disulfides
15-rosetta-modelling/run/rosetta-relax.sh -t e2 -s 3_hybridize/output_structures \
    -d e2.disulfide -M -n 5 -q
bin/qdm-submit 4_relax

# 5. check the symmetry definition cheaply, then dock
15-rosetta-modelling/run/rosetta-symm-dock.sh -p e2_best.pdb -N 3 -r 40 -S
15-rosetta-modelling/run/rosetta-symm-dock.sh -p e2_best.pdb -N 3 -r 40 -n 100 -q

# 6. impose the measured contacts
15-rosetta-modelling/run/rosetta-constrained-relax.sh -p e2_trimer.pdb \
    -c contacts.cst -n 1000 -W 4 -q

# 7. scores and plots at any stage
15-rosetta-modelling/extract/rosetta-scores.sh -o scores -k total_score 4_relax
15-rosetta-modelling/plot/plot-score-distribution.sh -i scores/scores.csv
```

The complete sequence, with what to check between stages, is in
[`../docs/running.md`](../docs/running.md#i-building-a-structure-with-rosetta).
Every script also prints its own usage block with `-h`.

## Inputs this module does not produce

Four things come from web services, and they are prerequisites rather than
optional extras:

| input | where it comes from | needed by |
|---|---|---|
| `TARGET_3.frags`, `TARGET_9.frags` | the Robetta fragment server | stage 3 |
| `TARGET.octopus` | TOPCONS or OCTOPUS, for a membrane protein | stage 1 with `-M` |
| the template alignment | a multiple aligner such as Clustal Omega | stage 1 |
| template PDB files | the PDB, searched for homologues | stage 2 |

Sequences must be upper case throughout. Rosetta reads lower case as a
deletion marker, so a sequence that came back from an aligner in mixed case
silently produces a truncated model. `rosetta-cm-setup.sh` warns and
uppercases its own copy.

One note from experience, recorded in the original project notes: for a
straightforward homology model, the Robetta server with its RoseTTAFold option
and no templates at all often beats a hand-built comparative model. This
pipeline earns its keep when you need control over which templates contribute
where, when the target is a membrane protein, or when the model has to satisfy
measured contacts.

## Restraints

`rosetta-constrained-relax.sh` takes Rosetta's own `AtomPair` syntax, or a
short form that it expands:

```text
# RES:ATOM   RES:ATOM   lower  upper
163A:H       167B:H     1.5    6.5     # amide proton contact between chains
176A:ND2     537B:CG    2.6    3.4     # interhelical hydrogen bond
179A:NZ      537B:CG    2.5    4.0     # salt bridge
```

Residue numbers carry a chain letter because pose numbering runs continuously
across chains, so residue 163 of chain B is not residue 163 in pose numbering.
PDB numbering with the chain says what you mean.

Before writing restraints, find the atoms in the structure you will actually
relax, not in a reference:

```bash
grep ' LYS ' complex.pdb | awk '$6 == 179 && $3 == "NZ"'
```

## Corrections made to the original scripts

**Cluster paths and account details are gone.** Every script began by sourcing
a specific `/mnt/lustre/groups/.../rosetta_bin_linux_2021.16.61629_bundle`
install, and the PBS headers carried a project code, personal output paths and
`-M sss` as an email address. Rosetta's location, build suffix and setup lines
now come from `config/cluster.local.conf`, which is gitignored, and job
submission goes through `bin/qdm-submit` like every other module.

**Model ranking was reading the wrong column.** Models were picked with
`cat *.fasc | sort -nk2` and `sort -nk5 | head -1 | awk '{print $NF}'`. Three
things break that, and all three were reachable:

- a Rosetta score file has a header row, and a header sorts numerically as
  zero, so with positive interface scores the "best model" read off the top is
  the literal word `description`, and the next command tries to copy
  `description.pdb`;
- column 2 is `total_score` and column 5 is `I_sc` only for one particular set
  of score terms. Drop one term and the same command silently ranks on
  whatever moved into that position. On a two-model test this picks a model
  that is worse by 3.5 energy units;
- newer Rosetta writes `.fasc` as one JSON object per line, where field 2 is
  the decoy name for every row, so `sort -nk2` does nothing at all and the
  "best" model is whichever was written first.

`extract/rosetta-scores.sh` looks terms up by name from each file's own header,
drops header rows, reads both the classic and the JSON layouts, and reports
when score files with different columns have been pooled.

**Interfaces were ranked on total score in places.** A pose can score better
overall while making a worse interface, which is the opposite of what a
docking run is being asked. The docking stages rank on `I_sc` and say so.

**Restraints that name a missing atom were silently dropped.** Rosetta does not
stop for this: the run completes, the scores look normal, and nothing was
enforced. The constrained relax now checks every atom in the restraint file
against the structure first and refuses to start otherwise. This matters most
for amide hydrogen restraints, because a structure straight from a crystal or
from some preparation steps has no hydrogens, so the restraint can never be
satisfied.

**`docking_flexible_symm` only scored.** Despite its name it called
`score.default.linuxgccrelease`, never a docking protocol, so it evaluated the
symmetric assembly without sampling it. It also ended with a bare line reading
`score_local_dock.sc`, which the shell tried to execute as a command, and then
sorted that file, which that script never produced; both were copied from the
other docking script. `rosetta-symm-dock.sh` runs `SymDock` by default and
scores only with `-S`, which is worth doing first because it is cheap and
checks the symmetry definition.

**Two poses is not a search.** Local docking ran with `-nstruct 2` under a 3
degree and 8 Angstrom perturbation, then took the better of the two. The
default is 10 now, and the help says hundreds is normal.

**The relax stage changed score function mid-pipeline.** Hybridize used
`stage3_rlx_membrane.wts` and the following relax used
`membrane_highres_Menv_smooth.wts`, so rankings from the two stages were not
comparable. The weights are an explicit argument and the value used is printed.

**24 MPI ranks were requested and the serial binary was run**, so 23 cores sat
idle for jobs of thousands of structures. The run scripts use an MPI build when
`ROSETTA_MPI_SUFFIX` is configured and say plainly when they are falling back
to one core.

**Template and disulfide lists lived inside the scripts**, so every new target
meant editing the code. They are arguments now: templates on the command line
as `6urh:C`, disulfides in a file that is checked for being two residue numbers
per line.

**Template sequences were matched with `grep -i NAME`**, which matches
substrings, so a template called `2kz` would also pull the sequence of `2kzq`
and the wrong structure would be threaded with no error. Names are compared as
whole fields, and the script also checks that every aligned sequence has the
same length, which a valid alignment must.

**`mv S_*.pdb output_structures/` ran unconditionally.** A crashed run ended
with `mv: cannot stat 'S_*.pdb'`, which says nothing about the real cause. The
original project notes record this error and guess at memory, which is often
right. The run script now reports that no models were produced and names the
three usual causes.

**Chain renumbering went through tleap**, which rebuilds residues and can
change atom names and protonation on the way. `rosetta-chains.sh` edits the
residue number field directly and verifies the atom count is unchanged. It also
writes one output file rather than two similarly named ones, only one of which
was renumbered.

**`rm ss2b` and `rm tmp/*` without `-f`** made the first run of a script fail
on a file that did not exist yet. Intermediate files are handled properly.

## Tested

`tests/test-15-rosetta.sh`, 24 checks. Rosetta itself is not installed in the
test environment, so what is tested is everything up to the point of launching
it, which is where the defects above lived:

- Grishin alignment generation from a Clustal alignment, including that a
  template whose name is a substring of another gets its own sequence, that an
  alignment with unequal lengths is rejected, and that a missing target name
  gives a usable error;
- score parsing for both file layouts, the header row dropped, lookup by term
  name, ranking on `total_score` and on `I_sc`, and detection of pooled files
  with mismatched columns;
- restraint expansion from the short form, and refusal when a restraint names
  an atom that is not in the structure;
- chain relabelling and per-chain renumbering with the atom count preserved;
- rejection of a multi-chain subunit for symmetric docking, and of a partner
  string naming a chain that is not present.

The Rosetta runs themselves are not reproduced here. The published models were
built with the 2021.16 release on a PBS cluster over days of wall time.

## Requirements

Rosetta, free for non-commercial users, from
[rosettacommons.org/download](https://rosettacommons.org/download/). Set
`ROSETTA_BIN`, `ROSETTA_SUFFIX` and optionally `ROSETTA_MPI_SUFFIX` and
`ROSETTA_TOOLS` in `config/cluster.local.conf`; see
[`../install/README.md`](../install/README.md).

Two Rosetta scripts are used when available and are not redistributed here,
because they are other people's work: `clean_pdb.py` by the Baker and Meiler
labs, and `octopus2span.pl` by Bjorn Wallner and Vladimir Yarov-Yarovoy. Both
ship with Rosetta under `main/tools/`. See
[`../docs/third-party-material.md`](../docs/third-party-material.md#rosetta).
Neither is required: templates can be cleaned by hand and `-O` is only needed
for a membrane protein.

## Scope

A comparative model is a hypothesis. The score function ranks models against
each other and does not measure how close any of them is to the truth, so a
lowest-scoring model is a best guess and not a structure. Where this module's
output matters most is as a starting point for the simulations in modules 08
to 12, and the honest way to report it is with the templates, their coverage
and identity, the sampling depth, and whatever independent evidence the
restraints came from.

The glycosylation and membrane steps that followed in the original work used
the GLYCAM glycoprotein builder and CHARMM-GUI, both web services, and are not
scripted here. `docs/running.md` records where they sit in the sequence.
