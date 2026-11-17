# flux2features

Code for *From metabolic fluxes to functional features: an interpretable
representation of constraint-based model states* (Preciat, Sotelo Rodriguez
and Gómez Márquez).

A flux vector of a constraint-based model, whether a flux balance optimum, a
sampled state or fluxes estimated from data, is a list of thousands of reaction
rates. This repository maps any such vector $v$ of a reconstruction $M$ to
$\Phi_M(v)$, a vector of declared, interpretable descriptors (pathway activity,
metabolite turnover, cofactor budgets, sources and allocation, oxidative defence,
production by compartment, elemental and amino-acid exchange, carbon and
nitrogen sources, capacity use and the coupling between inside and boundary),
each with its units, its requirements and the rule for when the reconstruction
lacks them (structurally unavailable, never a zero). The mapping needs only a
flux vector; the study evaluates it on states sampled from the human red blood
cell reconstruction RBC-GEM and from its thermodynamically flux-consistent
subset.

Two implementations serve two purposes. The **analysis engine**
(`python/flux2features_engine.py`) implements all 407 declared descriptors and
reproduces the paper. The **standalone tool** `flux2feature` (Python and MATLAB) is
the light interface for applying the representation to another reconstruction:
a Python package (`tools/python`) and a MATLAB function laid out for the COBRA
Toolbox (`tools/matlab`), with the same 179 descriptors. Version 1.0 of the
registry and of the tools is frozen with the article.

This folder is part of [GIMIpapers](https://github.com/GIMIudg/GIMIpapers); releases are
archived at Zenodo ([doi:10.5281/zenodo.19339596](https://doi.org/10.5281/zenodo.19339596)).

## Contents

| path | what it is |
|---|---|
| `matlab/generateModelsAndSamples_flux2features.m` | builds both reconstructions and samples every condition with three CHRR chains, each from its own fixed seed: the 3 × 3 design plus the reference in both models, and the reference and G6PDH2 0.10 of Original with NaKt held at no less than 95% and 90% of its maximum. Checks every chain against the condition it is named after and writes a SHA-256 fingerprint of each; `mode = 'check'` samples one chain and compares it with `matlab/sampling_fingerprints.csv` |
| `tools/python/` | the standalone tool `flux2feature` as a Python package (NumPy only; `pyproject.toml`, tests, three examples, its own README, changelog and citation); the notebook checks it against its engine (stage 28) |
| `tools/matlab/` | the same tool in MATLAB, laid out as in the COBRA Toolbox: `src/analysis/flux2feature/flux2feature.m`, two tests in `test/verifiedTests/analysis/testFlux2feature/` (toys worked by hand, the absence rule, the tiers, RBC-GEM against the reference values; *E. coli* core with an LP solver) and a tutorial |
| `tools/README.md` | the two tools side by side: the same 179 descriptors, tiers and reference values |
| `python/test_engine_release.py` | the tests of the release: every declared descriptor in one tier, every alias with a canonical form, the tool's tiers equal to the release's, and the engine on the reference state of RBC-GEM (NaN, never 0, where the model cannot support a descriptor) |
| `python/flux2features_engine.py`, `python/make_engine_module.py` | the analysis engine, copied from the notebook's own cells; the notebook checks it against itself on all 407 descriptors |
| `python/descriptor_release.csv` | the release: the tier of every declared descriptor, its canonical form if it repeats another, its module and the evidence for the tier |
| `tools/python/tests/data/`, `tools/matlab/test/verifiedTests/analysis/testFlux2feature/` | the reference state of RBC-GEM (reaction, flux), its 179 descriptor values and RBC-GEM reduced to the fields the tools read, written by the notebook after checking the tool against its engine to 1e-5 |
| `notebook/metabolic_flux_phenotyping.ipynb` | the whole analysis, in order: the integrity of every chain, the descriptors and their quality control, the response to the gradient, raw fluxes against the representation, parsimony, transfer between reconstructions, secondary validations, the audit, every figure (shown below the analysis it comes from) and supplementary table, and the build of the manuscript |
| `environment.yml`, `CITATION.cff` | the Python environment, and how to cite this code |
| `data/RBC-GEM.mat` | RBC-GEM v1.3.0 as used ([commit 43fd381](https://github.com/z-haiman/RBC-GEM/commit/43fd381d9eb8a8a4079629012362312960e1c009); MD5 `294801a126d7dd89977ee7c3061de4b5`). Cite Haiman ZB, Key A, D'Alessandro A, Palsson BO. *PLoS Comput Biol* 2025;21(3):e1012109 |

## The descriptor sets, in one table

The same numbers as Table 3 of the paper (RBC-GEM: Original / Pruned).

| Set | Declared | Defined | Default output | What it is |
|---|---:|---:|---:|---|
| Registry | 407 | 396 / 383 | — | every declared descriptor, with definition and requirements |
| QC-retained (core of the analysis) | — | 240 / 235 | — | kept by quality control for the analyses of the paper |
| Release core | 141 | 130 / 129 | — | the default coordinates of the released Φ, chosen after every analysis |
| Analysis engine (`python/`) | 407 | 396 / 383 | 141 | implements the registry; returns the release core by default |
| Standalone tool (`tools/`, Python and MATLAB) | 179 | 168 / 162 | 74 | 170 registry descriptors and 9 of its own; its default is the 72 release-core descriptors it implements and 2 of its own |

The tool returns 74 descriptors by default, not 141, because it implements only
the descriptors that need no curated metadata; the other 69 release-core
descriptors read the study's full pathway map, more named metabolites and
reactions, or the curated roles of exchanged metabolites.

## The release: what Φ returns by default

The registry declares 407 descriptors. After the audit of the study
(Supplementary Note 1 and Supplementary Table S8 of the paper) each one is in one tier,
and none is deleted:

| tier | what it is | returned |
|---|---|---|
| `core` | non-redundant functional coordinates with a biochemical reading (pathways, turnover, budgets, efficiencies, coupling, compartments, NADPH sources) | by default |
| `module` | descriptors whose use depends on the medium and the exchanges, the constraints, sparse states such as optima, several sources of a cofactor or an imposed oxidant load | on request |
| `generic` | statistics of the flux vector (norms, moments, concentration): features for learning, not functional coordinates | on request |
| `alias` | the same number as, or always redundant with, a canonical descriptor in the study | on request, under its canonical form |

`python/descriptor_release.csv` gives the tier, the canonical form and the
evidence for every descriptor. The tiers come from one reconstruction and one
perturbation. This release (v1.0) is frozen with the article; the same audit in
Recon3D is planned as a later validation (v1.1), not as a condition of v1.0.

The core is the set of interpretable coordinates the study releases, not a
claim that these are the most important descriptors of every system. Quality
control (which descriptors are retained for the analysis) and the release
(which tier each retained or removed descriptor goes to) are distinct stages.

## Using the tool

Install and full documentation: [`tools/python`](tools/python/README.md) and
[`tools/matlab`](tools/matlab/README.md).

```python
from flux2feature import Model, flux2feature, audit
model = Model(S, rxns, mets, lb, ub, subsystems)   # or Model.from_cobra(cobrapy_model)
audit(model)                                 # every declared descriptor: computable or not, and why
phi = flux2feature(model, v)                 # v: one flux vector, or a states-by-reactions matrix;
                                             # the core of the release by default
flux2feature(model, v, which=("core", "module"))   # add the modules
flux2feature(model, v, which="all")                # every declared descriptor
```

On the reference state of RBC-GEM shipped with the code
(`tools/python/tests/data/reference_state.csv`), the call returns, among the core:

```text
path_pentose_phosphate_activity           7013.83
turnover_6pgc_c                           1465.28
bud_oxppp_cover_of_glutathione_reduction     1.48
eff_ppp_flux_split                         328.62
cmp_flux_in_mitochondrion                     nan   # RBC-GEM has one compartment
```

and `audit(model)` says why the last one is missing:

```text
{'name': 'cmp_atp_made_in_mitochondrion', 'computable': False,
 'reason': 'atp is in 1 intracellular compartment(s) (c); needs m and one other', 'tier': 'core', ...}
```

```matlab
info = flux2feature(model);                          % the audit
f = flux2feature(model, samples);                    % samples: nRxns x nStates, as sampleCbModel returns
f = flux2feature(model, samples, 'Which', 'all');    % every declared descriptor
```

Over many feasible states, the distribution of a descriptor separates
properties the constraints fix from properties they leave free (metabolic
reserve). A range taken over sampled states is an empirical range of those
states, not a confidence interval and not a physiological limit. Φ does not
replace flux variability analysis, sampling diagnostics or pathway methods; it
can serve as their input or as a layer to read their output.

### What a model needs

| family | what it reads from the model |
|---|---|
| magnitude, distribution, sparsity, entropy, architecture, coupling, capacity | nothing named: the flux vector, the bounds, and which reactions are exchanges (`EX_` or one metabolite) |
| pathway | the subsystem of every reaction, grouped into the declared pathway groups |
| turnover, cofactor, efficiency | the declared metabolites (`nadph_c`, `atp_c`, ...) by identifier |
| named reactions, budgets | the declared reactions (`G6PDH2`, `HEX1`, ...) by identifier |
| elemental, source diversity, conversion | metabolite formulas of the exchanged metabolites |
| amino acids | exchange reactions of amino acids, recognised by BiGG or VMH stem or by name |
| compartments | metabolites in at least two intracellular compartments (`_c`, `_m`, ...) |

A descriptor whose requirement is missing is NaN with its reason, never 0.

### Another reconstruction: mapping identifiers

The defaults describe human erythrocyte metabolism with RBC-GEM (BiGG-style)
identifiers. For another reconstruction:

1. Run `audit(model)` first: it lists, for every descriptor, what is missing.
2. Map identifiers, not formulas. Pass your own declarations:
   `pathways` (pathway group → subsystem names), `named_reactions` (role →
   reaction identifier, e.g. `{"g6pd": "G6PDH2r"}`), `tracked_exchanges` and
   `cofactor_pairs` are merged with the defaults, so mapping never removes a
   coordinate; `turnover_metabolites` replaces them. In MATLAB the same names
   are name–value pairs (`'Pathways'`, `'NamedReactions'`, ...). Metabolite
   identifiers are read in BiGG style (`atp_c`; convert `atp[c]` first).
3. For another namespace (for example VMH or MetaNetX), translate the model's
   identifiers or the declarations first, with a mapping you check by hand on
   the reactions the core reads; descriptors that need nothing named work
   unchanged.
4. Run `audit(model)` again and compare the computable core with the one you
   expect; a coordinate that stays NaN tells you what is still unmapped.

Run the tests with `pytest` in `tools/python` (the tool) and `pytest python/`
(the release and the engine), and, in MATLAB, `testFlux2feature` (no solver)
and `testFlux2featureEcoliCore` (an LP solver) from
`tools/matlab/test/verifiedTests/analysis/testFlux2feature/`.

## Running the analysis

The code and the data are kept apart. The project folder holds `data/`
(`RBC-GEM.mat` and `metadata/`), `results/` and `reports/paper/`. The sampling
archives (about 13 GB) are not distributed: the MATLAB script regenerates them,
every chain from its own fixed seed.

1. **Sampling (MATLAB, COBRA Toolbox v3, Gurobi and/or MOSEK, about 9 hours).**
   Run `matlab/generateModelsAndSamples_flux2features.m` and choose the project
   folder. It writes `results/raw/{models,fluxes,samples,samples_objective}`,
   `results/raw/sampling_log.txt` and `results/raw/sampling_fingerprints.csv`
   (seed, rounding solver, software versions and SHA-256 of every chain).
   To check first that your installation reproduces the published chains, set
   `mode = 'check'`: it samples one chain (about 7 minutes) into a temporary
   folder and compares it with `matlab/sampling_fingerprints.csv`. A different
   version of MATLAB, the COBRA Toolbox or the solver can change the arithmetic
   of the rounding, and with it the chains; the comparison says so.
2. **Analysis (Python 3.11 or later; `conda env create -f environment.yml`).**
   Point the notebook at the project folder, either with the environment
   variable `FLUX2FEATURES_PROJECT=/path/to/project` or with a file
   `notebook/project_dir.txt` holding that path, and run it from top to bottom
   (about ten hours on two cores, half of it in the audit of stage 21). It reproduces every table, figure and number
   of the manuscript and builds the manuscript. The definition cells open
   collapsed, so the notebook reads as text, runs and results; set
   `RECOMPUTE = False` in the first code cell to re-display a previous run
   without recomputing it.

## Licence

MIT; see `LICENSE`.
