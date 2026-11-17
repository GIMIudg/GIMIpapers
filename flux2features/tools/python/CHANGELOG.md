# Changelog

All notable changes to `flux2feature`. Versions follow [Semantic Versioning](https://semver.org).

## 1.0.0 — first public release

- 179 declared descriptors in 22 families, each with its name, units, block
  (internal, external, combined), what it needs from the model and what it means.
- Every descriptor carries its tier of the release (core, module, generic, alias),
  from the audit of the accompanying study; the default output is the core (74
  descriptors).
- `audit(model)` reports, before anything is computed, which descriptors a model
  supports and why the others cannot be computed.
- A descriptor whose requirement is missing returns NaN, never 0.
- The same descriptors and tiers in the MATLAB twin (`flux2feature.m`, COBRA
  Toolbox layout), tested against the same reference values.
- Checked against the study's engine on a sampled state of RBC-GEM; the reference
  values ship with the tests of both languages.

This release is frozen with the article. Validation on further reconstructions
(Recon3D among them) is planned for 1.1 and will not change the 1.0 coordinates.
