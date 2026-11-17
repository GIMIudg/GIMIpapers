# flux2feature (MATLAB, COBRA Toolbox)

`flux2feature.m` maps flux vectors of a COBRA model (a flux balance solution,
the output of `sampleCbModel`, any matrix of states) to declared, interpretable
descriptors, and audits which descriptors the model can support. One file, base
MATLAB (R2016b or later), no toolbox; the tests and the tutorial use the COBRA
Toolbox. Version 1.0.0.

It implements the same 179 descriptors, the same tiers and the same values as
the Python package in `../python`, and its test checks the same reference values
of a sampled state of RBC-GEM, to 1e-5 relative, that the Python tests check.

The folders follow the layout of the COBRA Toolbox, so that they can be copied
into it as they are:

| Here | In the COBRA Toolbox |
|---|---|
| `src/analysis/flux2feature/flux2feature.m` | `cobratoolbox/src/analysis/flux2feature/` |
| `test/verifiedTests/analysis/testFlux2feature/` (two tests, the reference state and values, RBC-GEM reduced to the fields the tool reads, ~100 kB) | `cobratoolbox/test/verifiedTests/analysis/testFlux2feature/` |
| `tutorials/analysis/flux2feature/tutorial_flux2feature.m` | `COBRA.tutorials/analysis/flux2feature/` |

## Usage

```matlab
info = flux2feature(model);                         % the audit: what is computable, and why not
f = flux2feature(model, samples);                   % the core of Phi, one value per state
f = flux2feature(model, samples, 'Which', {'core', 'module'});
f = flux2feature(model, samples, 'Which', 'all');   % all 179 declared descriptors
[f, info, cmp] = flux2feature(model, samples, vReference);   % differences and log2 ratios
```

`samples` may be `nRxns x nStates` (as `sampleCbModel` returns it) or one flux
vector. A descriptor the model cannot support is `NaN`, never 0, and its reason
is in `info`.

The defaults name human erythrocyte metabolism (RBC-GEM). On another
reconstruction, map the labels to its identifiers; the labelled declarations
are merged with the defaults, so no coordinate is lost:

```matlab
model.mets = strrep(regexprep(model.mets, '\[(\w+)\]$', '_$1'), '-', '__');   % atp[c] -> atp_c
f = flux2feature(model, samples, ...
    'NamedReactions', {'g6pd', 'G6PDH2r'; 'ldh', 'LDH_D'}, ...
    'Pathways', {'pentose_phosphate', {'Pentose Phosphate Pathway'}});
```

## Tiers

| Tier | Descriptors | Meaning |
|---|---:|---|
| `core` | 74 | functional coordinates with a biochemical reading (the default output) |
| `module` | 61 | descriptors whose use depends on the model, the medium or the kind of state |
| `generic` | 18 | statistics of the flux vector; features for learning, not functional coordinates |
| `alias` | 26 | the same number as, or always redundant with, a canonical descriptor |

The tool's core is the 72 descriptors of the paper's release core (141) that the tool
implements, plus 2 of its own; the other 69 release-core descriptors need the
study's full pathway map, more named metabolites and reactions, or curated roles
of exchanged metabolites, and are computed by the analysis engine.

The core is the study's release of interpretable coordinates, not a claim that
they are the most important descriptors of every system. Over sampled states,
the spread of a descriptor separates properties the constraints fix from those
they leave free; a range over sampled states is an empirical range, not a
confidence interval or a physiological limit.

## Tests and tutorial

```matlab
initCobraToolbox(false)
cd test/verifiedTests/analysis/testFlux2feature
testFlux2feature            % no solver: toys worked by hand, the absence rule, the tiers,
                            % and RBC-GEM against the Python reference values (179 descriptors)
testFlux2featureEcoliCore   % needs an LP solver: E. coli core, identifiers mapped, values
                            % against the same quantities computed from the model
```

The tutorial (`tutorial_flux2feature.m`, MATLAB cell format) applies the tool to
*E. coli* core: the audit, a parsimonious optimum, the distribution of the
descriptors over states sampled with CHRR, and a perturbed state against the
reference. Open it in MATLAB and save it as a live script (`.mlx`) for
COBRA.tutorials.

## Submitting to the COBRA Toolbox

Contributions go through pull requests from a fork, usually with
[MATLAB.devTools](https://github.com/opencobra/MATLAB.devTools):

1. In MATLAB: `contribute('cobratoolbox')`, option 1, to start a branch (for
   example `add-flux2feature`) in your fork.
2. Copy `src/analysis/flux2feature/` and
   `test/verifiedTests/analysis/testFlux2feature/` into the same paths of
   `fork-cobratoolbox`.
3. Run both tests there, then `contribute('cobratoolbox')`, option 3, to push
   and open the pull request against `develop`.
4. The tutorial goes to `COBRA.tutorials` the same way
   (`contribute('COBRA.tutorials')`), as `analysis/flux2feature/tutorial_flux2feature.mlx`
   together with this `.m` file.

## Cite

Preciat G, Sotelo Rodriguez R, Gómez Márquez CE, *From metabolic fluxes to
functional features: an interpretable representation of constraint-based model
states* (and the archived release of the code). RBC-GEM: Haiman ZB, Key A,
D'Alessandro A, Palsson BO (2025) *PLoS Comput Biol* 21(3): e1012109. MIT
license.
