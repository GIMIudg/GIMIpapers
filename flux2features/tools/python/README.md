# flux2feature (Python)

Interpretable, declared descriptors of metabolic flux vectors.

Give `flux2feature` a constraint-based model and one or more flux vectors (a
flux balance solution, a sampled state, a matrix of sampled states) and it
returns named biochemical coordinates instead of thousands of reaction fluxes:
pathway activity and share, metabolite turnover, cofactor allocation, NADPH and
glutathione budgets, conversion yields, boundary exchange, compartment shares.
Before computing anything it tells you which descriptors your model can
support and, for the others, why not.

It is the reusable part of the representation Φ_M(v) of the accompanying study
and holds 179 descriptors: 170 of the 407 declared there (the ones that do not need
the study's curated metadata) and 9 simpler variants of its own. The MATLAB twin in `../matlab` (COBRA Toolbox
layout) implements the same 179 descriptors, the same tiers and the same
reference values.

- One module, one dependency (NumPy). pandas is optional (`as_frame=True`).
- Three rules: every descriptor is **declared** before it is computed (name,
  units, block, requirement, reading); a descriptor whose requirement is missing
  is **NaN, never 0**; nothing is renamed, reordered or dropped silently.

## Install

```bash
pip install flux2feature                       # from PyPI, once released
# or straight from the repository
pip install "git+https://github.com/GIMIudg/GIMIpapers.git#subdirectory=flux2features/tools/python"
# or from a clone, with the test dependencies
pip install -e "tools/python[test]"
```

Python ≥ 3.9 and NumPy ≥ 1.22. Optional: pandas (`[frame]`), COBRApy (`[cobra]`),
pytest and SciPy (`[test]`).

## Minimal example

```python
import numpy as np
import flux2feature as f2f

model = f2f.Model(S, rxns, mets, lb, ub, subsystems)  # or f2f.Model.from_cobra(cobra_model)

rows = f2f.audit(model)                    # what is computable here, and why the rest is not
phi = f2f.flux2feature(model, v)           # the core of Phi: {descriptor: value}
phi = f2f.flux2feature(model, V)           # a matrix of states: one value per state
phi = f2f.flux2feature(model, V, which=("core", "module"))
phi = f2f.flux2feature(model, V, which="all", as_frame=True)   # every descriptor, a DataFrame
```

`examples/quickstart.py` runs this on a 15-reaction toy erythrocyte and prints:

```text
179 descriptors declared, 78 computable here
  core      36 of 74 computable
  ...
rxn_pgl is not computable here: PGL is not in this model

the core: 74 coordinates, 36 finite on this model
  path_pentose_phosphate_share        0.03409
  turnover_nadph_c                    0.4
  eff_nadph_per_glucose               0.2
  rxn_gpx                             0.2
  conv_lactate_per_glucose_carbon     nan
                                     (no glucose or lactate exchange, or no formulas with carbon)
  cmp_atp_made_in_mitochondrion       nan
```

| Example | What it shows | Needs |
|---|---|---|
| `examples/quickstart.py` | the audit, the core, many states, mapping an identifier | NumPy |
| `examples/rbc_gem_reference.py` | a genome-scale model (RBC-GEM) and the shared reference values | SciPy |
| `examples/ecoli_core_cobrapy.py` | another organism through COBRApy: pFBA and sampling of *E. coli* core | COBRApy |

## The tiers

Every descriptor carries the tier of the release decided by the audit of the
study (Supplementary Table S8). The default output is the core.

| Tier | In this tool | Meaning |
|---|---:|---|
| `core` | 74 | functional coordinates with a biochemical reading |
| `module` | 61 | descriptors whose use depends on the model, the medium or the kind of state |
| `generic` | 18 | statistics of the flux vector (norms, moments, concentration); features for learning, not functional coordinates |
| `alias` | 26 | the same number as, or always redundant with, a canonical descriptor; computed only on request |

`f2f.STATUS[name]` gives `(tier, canonical)`; `audit()` repeats both per descriptor.

The tool's core is the 72 descriptors of the paper's release core (141) that the tool
implements, plus 2 of its own; the other 69 release-core descriptors need the
study's full pathway map, more named metabolites and reactions, or curated roles
of exchanged metabolites, and are computed by the analysis engine.

## What a model needs

Families that depend on nothing named work on any model. The others need the
named part of the model, and are NaN with a reason when it is absent.

| Families (descriptors) | Needs |
|---|---|
| magnitude, sparsity, entropy, architecture, distribution, network (20) | the stoichiometric matrix |
| pathway (32) | subsystem annotations, or reaction identifiers for groups defined by identity |
| turnover (20), cofactor (3), efficiency (7) | the named metabolites |
| biochemical (11), oxidative defence (10), budget (2) | the named reactions and signatures |
| elemental (12), source diversity (5), amino acids (22), conversion (4) | metabolite formulas (and names, for amino acids) |
| compartment (11) | at least two intracellular compartments (`_c`, `_m`, …) |
| capacity (4) | finite bounds |
| exchange burden (4), boundary metabolites (9), coupling (3) | boundary reactions (`EX_` prefix, or `exchange=`) |

## Another reconstruction

The defaults name human erythrocyte metabolism (RBC-GEM, BiGG-style
identifiers). Map the labels to your model's identifiers; the coordinates keep
their names, so states of different models stay comparable:

```python
declaration = dict(
    named_reactions={"g6pd": "G6PDH2r", "ldh": "LDH_D"},
    pathways={"glycolysis": ("Glycolysis/Gluconeogenesis",)},
    tracked_exchanges={"lac__D_e": "lactate"},
)
f2f.audit(model, **declaration)
f2f.flux2feature(model, V, **declaration)
```

`pathways`, `pathway_reactions`, `named_reactions`, `tracked_exchanges` and
`cofactor_pairs` are merged with the defaults: a label you give replaces its
default identifiers and the others are kept, so mapping never removes a
coordinate. `turnover_metabolites` and `currency` replace the defaults. Pass
`exchange=` to `Model` if your boundary reactions do not start with `EX_`.
Metabolite identifiers are read in BiGG style (`atp_c`); convert `atp[c]` first.

## Reading the values

- The core is the release's set of interpretable functional coordinates for the
  study; it is not a claim that these are the most important descriptors of any
  system.
- Computed over many feasible states (a sample), the spread of a descriptor
  separates properties the constraints fix from properties they leave free
  (metabolic reserve). A range taken over sampled states is an empirical range
  of those states, not a confidence interval and not a physiological limit.
- Φ does not replace flux variability analysis, sampling diagnostics or pathway
  methods; it can serve as their input or as a layer to read their output.

## Tests

```bash
cd tools/python
pip install -e ".[test]"
pytest
```

The tests check every family on toy networks worked out by hand, the absence
rule, the tiers, and a sampled state of RBC-GEM (`tests/data/`) against the
reference values that the study's engine and the MATLAB tool also reproduce
(179 descriptors, to 1e-5 relative). `python tests/test_flux2feature.py` runs
them without pytest.

## Cite

Please cite the article and the archived release (`CITATION.cff`). RBC-GEM:
Haiman ZB, Key A, D'Alessandro A, Palsson BO (2025) *PLoS Comput Biol* 21(3):
e1012109. The file `tests/data/RBC-GEM_flux2feature.mat` is RBC-GEM reduced to
the fields the tool reads.

MIT license. Changes are listed in `CHANGELOG.md`. Version 1.0.0 is frozen with
the article; validation on further reconstructions (Recon3D among them) is
planned for 1.1 and will not change the 1.0 coordinates.

## Releasing (maintainers)

```bash
cd tools/python
python -m pip install build twine
python -m build                      # dist/flux2feature-1.0.0.tar.gz and .whl
python -m twine check dist/*
python -m twine upload dist/*        # PyPI account and token required
```

Tag the repository (`flux2feature-v1.0.0`), publish a GitHub release, and let
the Zenodo–GitHub integration archive it; add the DOI to `CITATION.cff`.
