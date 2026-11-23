# OEIS A386829, FORMULA 5, revisions 11–12 (2 November 2026)

The note proves Mathar's recurrence in A386829, FORMULA 5, revision 11
(approved as revision 12), for every
integer index n >= 2 and gives its algebraic ordinary generating function.
It also proves the displayed recurrence conjectures in A386830 for n >= 2
and A386831 for n >= 4. A two-weight algebraic family gives the finite-sum,
coefficient, exponential/reversion and hypergeometric formulas in related
entries. The coverage table treats twelve entries: A386829, A385669,
A386763, A386830, A386831, A385670, A385671, A386764, A386765, A386766,
A386769 and A383888. The last entry's stated asymptotic is also proved.

## Files

- `note.tex`: self-contained LaTeX source, including the bibliography.
- `note.pdf`: typeset mathematical note.
- `metadata.json`: authors, license, problem identifier and coverage.
- `verify.py`: exact integer coefficient checks using Python's standard library.
- `verification.json`: the checker's output and finite ranges.
- `verify_polynomials.py`: exact polynomial identities proving the added recurrences.
- `polynomial-identities.json`: their zero polynomial numerators and recurrence onsets.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, use a LaTeX installation containing the packages named in
`note.tex` and run `pdflatex note.tex` twice. The bibliography is embedded;
no separate `.bib` file or BibTeX invocation is needed.

## Reproduce the finite checks

```sh
python3 verify.py > verification.json
python3 verify_polynomials.py > polynomial-identities.json
```

The checker tests exact integer identities through n=80 for A386829 and
through n=30 for nine parameter choices, the shifted and reversion identities,
the hypergeometric formulas, and the added recurrences. Both scripts use
Python's standard library. The polynomial checker expands whole polynomials
with integer coefficients; it does not sample their values.
The note supplies the proofs for all indices.
The pinned Lean declaration
[`D5.S3.ArithSums.MatharA386829.result`](https://github.com/the-omega-institute/trureturing/blob/074df09ddbbef1843ddb383577bfaa827f753e36/D5/S3/ArithSums/MatharA386829.lean)
states the full literal recurrence and existence of a nonzero homogeneous
differential equation. The closed generating functions and sibling statements
are established by the written proofs in the note.

## Authors and license

Haobo Ma (ChronoAI Pte Ltd, Singapore), followed by Wenlin Zhang
(National University of Singapore; ChronoAI Pte Ltd, Singapore).
Creative Commons Attribution 4.0 International, CC BY 4.0 (`cc-by-4.0`).
