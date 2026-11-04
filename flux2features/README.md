# flux2features

Code for *An auditable representation of sampled metabolic states: compression,
parsimony and comparability across reconstructions* (Preciat, Sotelo Rodriguez
and Gómez Márquez).

A flux vector of a constraint-based model, whether a flux balance optimum, a
sampled state or fluxes estimated from data, is a list of thousands of reaction
rates. This repository maps any such vector $v$ of a reconstruction $M$ to
$\Phi_M(v)$, a vector of declared, interpretable descriptors (pathway activity,
metabolite turnover, cofactor allocation and budgets, production by compartment,
capacity use, exchange and the coupling between inside and boundary), each with
its units, its requirements and the rule for when the reconstruction lacks
them (structurally unavailable, never a zero). The mapping needs a flux vector,
not a sampler; the study evaluates it on states sampled from the human red blood
cell reconstruction RBC-GEM and from its thermodynamically flux-consistent
subset.

`flux2feature` is an independent tool, not part of the COBRA Toolbox.

This folder is part of [GIMIpapers](https://github.com/GIMIudg/GIMIpapers); releases are
archived at Zenodo ([doi:10.5281/zenodo.19339596](https://doi.org/10.5281/zenodo.19339596)).

## Contents

| path | what it is |
|---|---|
| `matlab/generateModelsAndSamples_flux2features.m` | builds both reconstructions, imposes the 3 × 3 design plus the reference, and samples every condition with three CHRR chains, checking every perturbed chain against its bound before saving it |
| `matlab/flux2feature.m` | the representation packaged for reuse, base MATLAB, no toolbox |
| `python/flux2feature.py` | the same tool in Python (NumPy only); the notebook checks it against its engine (stage 28) |
| `python/test_flux2feature.py`, `matlab/test_flux2feature.m` | the tests: every descriptor of a three-reaction model and of a two-compartment network worked by hand, the absence rule, and a real state of RBC-GEM against reference values |
| `matlab/reference_state.csv`, `matlab/flux2feature_reference.csv` | that state (reaction, flux) and its descriptors, written by the notebook after checking the tool against its engine to 1e-5 |
| `validation/` | an independent validation of both implementations: the same small models and states mapped in Python and MATLAB, properties any correct implementation must have, and a descriptor-by-descriptor cross-check (see `validation/README.md`) |
| `notebook/metabolic_flux_phenotyping.ipynb` | the whole analysis, in order: the integrity of every chain, the descriptors and their quality control, the response to the gradient, raw fluxes against the representation, parsimony, transfer between reconstructions, secondary validations, the audit, every figure (shown below the analysis it comes from) and supplementary table, and the build of the manuscript |
| `environment.yml`, `CITATION.cff` | the Python environment, and how to cite this code |
| `data/RBC-GEM.mat` | RBC-GEM v1.3.0 as used ([commit 43fd381](https://github.com/z-haiman/RBC-GEM/commit/43fd381d9eb8a8a4079629012362312960e1c009); MD5 `294801a126d7dd89977ee7c3061de4b5`). Cite Haiman ZB, Key A, D'Alessandro A, Palsson BO. *PLoS Comput Biol* 2025;21(3):e1012109 |

## Using the tool

```python
from flux2feature import Model, flux2feature, audit
model = Model(S, rxns, mets, lb, ub, subsystems)   # or Model.from_cobra(cobra_model)
audit(model)                 # every declared descriptor, computable or not, and why
phi = flux2feature(model, v) # v: one flux vector, or a states-by-reactions matrix
```

```matlab
info = flux2feature(model);          % the audit
f = flux2feature(model, v);          % v: nRxns x 1, or nRxns x nStates as sampleCbModel returns
```

Run the tests with `python python/test_flux2feature.py` and, in MATLAB,
`test_flux2feature` from `matlab/`.

## Running the analysis

The code and the data are kept apart. The project folder holds `data/`
(`RBC-GEM.mat` and `metadata/`), `results/` and `reports/paper/`; the sampling
archives are several gigabytes and are deposited separately.

1. **Sampling (MATLAB, COBRA Toolbox v3, Gurobi and/or MOSEK).** Run
   `matlab/generateModelsAndSamples_flux2features.m` and choose the project
   folder. It writes `results/raw/{models,fluxes,samples}` and
   `results/raw/sampling_log.txt` (solver used for the rounding, seed and bound
   check of every chain).
2. **Analysis (Python 3.11 or later; `conda env create -f environment.yml`).**
   Point the notebook at the project folder, either with the environment
   variable `FLUX2FEATURES_PROJECT=/path/to/project` or with a file
   `notebook/project_dir.txt` holding that path, and run it from top to bottom
   (about four hours on two cores). It reproduces every table, figure and number
   of the manuscript and builds the manuscript. The definition cells open
   collapsed, so the notebook reads as text, runs and results; set
   `RECOMPUTE = False` in the first code cell to re-display a previous run
   without recomputing it.

## Licence

MIT; see `LICENSE`.
