# flux2feature: the tools

The representation Φ_M(v) of the study, packaged for other reconstructions, in
two languages that implement the same 179 descriptors, the same tiers and the
same values (version 1.0.0).

| Folder | What it is | Publish to |
|---|---|---|
| [`python/`](python/README.md) | the Python package `flux2feature` (NumPy only): `pyproject.toml`, tests, three examples | PyPI and a tagged GitHub release archived by Zenodo |
| [`matlab/`](matlab/README.md) | `flux2feature.m`, its two tests and a tutorial, laid out as in the COBRA Toolbox | `opencobra/cobratoolbox` (function and tests) and `opencobra/COBRA.tutorials` (tutorial) |
| `ci/` | a GitHub Actions workflow for the Python tests | `GIMIpapers/.github/workflows/` |

## One set of numbers

The analysis notebook (stage 28) checks the Python tool against the study's
engine on a sampled state of RBC-GEM and writes, into the test folder of each
tool, that state, its 179 descriptor values and RBC-GEM reduced to the fields
the tools read. The Python and the MATLAB tests both check those values to
1e-5 relative, so the engine and the two tools are held to one set of numbers.

| | Python | MATLAB |
|---|---|---|
| declared descriptors | 179 | 179 |
| core / module / generic / alias | 74 / 61 / 18 / 26 | 74 / 61 / 18 / 26 |
| default output | the core | the core |
| a requirement missing | NaN, with the reason in `audit()` | NaN, with the reason in `info` |
| identifiers of another model | labelled declarations, merged with the defaults | the same, as name-value pairs |
| reference check (RBC-GEM) | `tests/test_flux2feature.py` | `testFlux2feature.m` |
| another organism | `examples/ecoli_core_cobrapy.py` | `testFlux2featureEcoliCore.m`, tutorial |

The tool's core is the 72 descriptors of the paper's release core (141) that the tool
implements, plus 2 of its own; the other 69 release-core descriptors need the
study's full pathway map, more named metabolites and reactions, or curated roles
of exchanged metabolites, and are computed by the analysis engine.

The engine that reproduces the paper (all 407 descriptors) is in `../python/`;
the tools hold 179: the 170 registry descriptors that need no metadata curated
for the study, and 9 simpler variants of their own.
Version 1.0 is frozen with the article; validation on further reconstructions
(Recon3D among them) is planned for 1.1 and will not change the 1.0 coordinates.
