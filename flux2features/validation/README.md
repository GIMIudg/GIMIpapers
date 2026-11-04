# Validation of flux2feature

An independent check, written to be run by someone who did not write the tool,
that the Python and MATLAB implementations of `flux2feature` compute what they
declare and agree with each other. Nothing here changes the tool: a failing check
is a finding to document and discuss before anything is modified.

| file | what it does |
|---|---|
| `cases/` | small models as plain CSV (reactions with bounds and subsystem, stoichiometry as triplets, one flux vector), read identically by both languages: a three-reaction model and a two-compartment network written with BiGG (`atp_m`) and VMH (`atp[m]`) identifiers |
| `validate_python.py` | the unit tests, every case and two states of RBC-GEM, and properties any correct implementation must have (orientation, batches, reordering, scaling, the absence rule, refusal of bad input, identifier conventions) |
| `validate_matlab.m` | the same, in MATLAB |
| `cross_check.py` | Python against MATLAB, descriptor by descriptor, on every case |
| `results/` | written by the three scripts: the values per case and language, and one report each |

## Running it

From the `flux2features` folder:

```bash
python validation/validate_python.py      # Python 3.11+, NumPy; SciPy to read RBC-GEM
```

```matlab
run('validation/validate_matlab.m')       % base MATLAB, no toolbox
```

```bash
python validation/cross_check.py          # after both
```

Each prints PASS or FAIL per check and writes a report in `results/`
(`python_report.md`, `matlab_report.md`, `cross_check_report.md`) with the
versions used. Record the date, machine, and versions of Python, NumPy, SciPy
and MATLAB with the reports.

## What the properties mean

- **Orientation and batches.** One state as a vector, as a row or as a column,
  and a matrix of states, give the same values state by state.
- **Reordering.** Permuting the reactions and metabolites of the model (and the
  flux vector with them) changes no descriptor: nothing depends on position.
- **Scaling.** Doubling every flux doubles every descriptor in flux units and
  leaves ratios, shares, entropies and statistics of shape unchanged. Counts
  (which compare fluxes with a tolerance) and capacity use (which compares them
  with bounds) are excluded by design.
- **Absence.** Every declared descriptor is returned with `which="all"`; one the
  model cannot support is NaN with a reason, never 0.
- **Refusal.** Duplicate identifiers, non-finite fluxes and a vector of the wrong
  length raise an error instead of returning numbers.
- **Identifier conventions.** The compartment descriptors are the same whether a
  metabolite is written `atp_m` or `atp[m]`. Descriptors that name BiGG
  identifiers (`turnover_atp_c`, `rxn_g6pd`) are declared unavailable under
  other identifiers until a mapping is supplied, never silently wrong; the
  network descriptors agree once the currency list is given in the same
  namespace.
