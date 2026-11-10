# Test fixtures

Small hand-made output fragments for program sections that have no public
reference output to download. They copy the line layout documented in the
program manuals and tutorials. Except where noted, the numbers are invented
and only test the parsers; they are not results.

| file | layout taken from | present |
|---|---|---|
| `orca-neb-fragment.out` | ORCA 6.1 tutorial "Finding transition states with NEB-TS" (path summary table, `<= CI` / `<= TS` markers) | **no** |
| `orca-esd-fragment.out` | ORCA 6.0 manual, "Excited State Dynamics" (rate constant line and FC/HT line) | **no** |
| `nbo-neda-fragment.out` | water dimer NEDA summary in the NBO 6 tutorial, https://nbo6.chem.wisc.edu/tut_neda.htm (numbers as published there) | **no** |
| `aimall-bcp-fragment.sum` | three invented BCPs (covalent, metal-ligand, hydrogen bond) in the `key = value` layout the original extraction script searched for in AIMAll `.sum`/`.mgpviz` files (`Type =`, `Rho =`, `DelSqRho =`, `Bond Ellipticity =`, `G =`, `V =`); no public AIMAll output was available | yes |
| `simulate-cv.awk` | generator of simulated voltammograms with analytically known E1/2 and onset, used by the electrochemistry tests; real potentiostat files were not available | yes |

## Why three of them are missing, and how to restore them

The three `.out` fixtures never reached the public repository. The cause was a
blanket `*.out` rule in `.gitignore`, written to keep simulation output from
being committed, which also matched these fixtures. Git excluded them without
reporting anything, and the result was that `test-02-orca.sh` and
`test-03-gamess.sh` failed on a fresh clone with an awk error about a file that
does not exist. The `.sum` and `.awk` fixtures survived only because those
extensions are not in the ignore list.

Two changes were made so this cannot recur:

* `.gitignore` ends with `!tests/fixtures/**`, which must stay **below** the
  `*.out` and `*.log` rules. Git applies the last matching pattern, so a
  negation placed above them has no effect at all: `*.out` wins and the
  fixtures are dropped exactly as before. Check it with
  `git check-ignore -v tests/fixtures/orca-neb-fragment.out`, which should
  report the negation rather than `*.out`;
* `tests/testlib.sh` gained a `fixture` helper. A test that cannot find its
  fixture now prints `missing fixture tests/fixtures/NAME, checks that need it
  skipped` and carries on with the rest of its checks, instead of crashing. A
  missing fixture is a packaging problem and should not look like a parser
  failure.

**To restore them,** copy the three files back into this directory from a local
copy of the project and commit them. `git add` will now accept them. They are
most likely still in the packaged archive of this repository, since the
exclusion happened at commit time rather than when the files were written.

They were not regenerated from the parsers on purpose. A fixture written to
match the regular expression that reads it tests nothing: it would pass whether
or not the parser handles real program output. These fragments are only worth
having if they follow the layout the programs actually print, which means
copying them from the manuals and tutorials cited above.

## Adding a fixture

Keep it short, a few dozen lines at most, and covering only the section being
parsed. Record in the table above where the layout came from, and say whether
the numbers are real or invented. If they are real and published, cite the
source, because a test that checks a published value is worth more than one
that checks an invented one.
