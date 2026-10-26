# flux2features

Code for *An auditable representation of sampled metabolic states: compression,
parsimony and comparability across reconstructions* (Preciat, Sotelo Rodriguez
and Gómez Márquez).

A sampled flux vector is a complete metabolic state, but a list of thousands of
reaction rates is not a description of it. This repository maps every sampled
state $v$ of a constraint-based model to $\Phi(v)$, a vector of declared,
interpretable descriptors (pathway activity, metabolite turnover, cofactor
allocation, capacity use, exchange and the coupling between inside and
boundary), each with its units, its requirements and the rule for when a
requirement is absent (a missing value, never a zero), and tests that
representation on the human red blood cell reconstruction RBC-GEM and on its
thermodynamically consistent subset.

This folder is part of [GIMIpapers](https://github.com/GIMIudg/GIMIpapers); releases are
archived at Zenodo ([record 19339596](https://zenodo.org/records/19339596)).

## Contents

| path | what it is |
|---|---|
| `matlab/generateModelsAndSamples_flux2features.m` | builds both reconstructions, imposes the 3 × 3 design plus the reference, and samples every condition with three CHRR chains, checking every perturbed chain against its bound before saving it |
| `matlab/flux2feature.m` | the representation on its own, base MATLAB, no toolbox |
| `matlab/test_flux2feature.m`, `matlab/flux2feature_reference.csv` | its test against reference values from the Python engine |
| `notebook/metabolic_flux_phenotyping.ipynb` | the whole analysis: export of the MATLAB archives, sampling diagnostics, descriptors, quality control, every analysis, figure and supplementary table, and the build of the manuscript |
| `data/RBC-GEM.mat` | RBC-GEM v1.3.0 as used ([commit 43fd381](https://github.com/z-haiman/RBC-GEM/commit/43fd381d9eb8a8a4079629012362312960e1c009); MD5 `294801a126d7dd89977ee7c3061de4b5`). Cite Haiman ZB, Key A, D'Alessandro A, Palsson BO. *PLoS Comput Biol* 2025;21(3):e1012109 |

## Running it

The code and the data are kept apart. The project folder holds `data/`
(`RBC-GEM.mat` and `metadata/`), `results/` and `reports/paper/`; the sampling
archives are several gigabytes and are deposited separately.

1. **Sampling (MATLAB, COBRA Toolbox v3, Gurobi and/or MOSEK).** Run
   `matlab/generateModelsAndSamples_flux2features.m` and choose the project
   folder. It writes `results/raw/{models,fluxes,samples}` and
   `results/raw/sampling_log.txt` (solver used for the rounding, seed and bound
   check of every chain).
2. **Analysis (Python 3.11 or later, NumPy, SciPy, pandas, scikit-learn,
   matplotlib; pandoc for the manuscript).** Point the notebook at the project
   folder, either with the environment variable
   `FLUX2FEATURES_PROJECT=/path/to/project` or with a file
   `notebook/project_dir.txt` holding that path, and run it from top to bottom.
   The last cell, `run_all()`, reproduces every table, figure and number of the
   manuscript; set `RUN_EVERYTHING = False` to load the definitions only and
   call any stage on its own.

## Licence

MIT; see `LICENSE`.
