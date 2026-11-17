"""flux2feature on another organism, through COBRApy: E. coli core.

    pip install cobra
    python examples/ecoli_core_cobrapy.py

The defaults of flux2feature name human erythrocyte metabolism. On another
reconstruction, map the labels to that model's identifiers; the coordinates
keep their names, and whatever the model does not have stays NaN, with the
reason in the audit. This example maps the BiGG identifiers of E. coli core,
solves a parsimonious FBA and draws a few hundred feasible states.
"""
import numpy as np
from cobra.flux_analysis import pfba
from cobra.io import load_model
from cobra.sampling import sample

import flux2feature as f2f

ecoli = load_model("textbook")                   # E. coli core, shipped with COBRApy
model = f2f.Model.from_cobra(ecoli)

# the same labels, E. coli core's identifiers (merged with the defaults, so
# every coordinate is kept and the ones E. coli core lacks stay NaN)
declaration = dict(
    pathways={"glycolysis": ("Glycolysis/Gluconeogenesis",),
              "pentose_phosphate": ("Pentose Phosphate Pathway",),
              "transport": ("Transport, Extracellular",)},
    named_reactions={"g6pd": "G6PDH2r", "ldh": "LDH_D"},
    tracked_exchanges={"lac__D_e": "lactate"},
)

rows = f2f.audit(model, **declaration)
core = [r for r in rows if r["tier"] == "core"]
print(f"{sum(r['computable'] for r in core)} of {len(core)} core descriptors computable "
      f"on E. coli core; the others, and why:")
for r in core:
    if not r["computable"]:
        print(f"  {r['name']}: {r['reason']}")

v = pfba(ecoli).fluxes.reindex([r.id for r in ecoli.reactions]).to_numpy()
phi = f2f.flux2feature(model, v, **declaration)
print("\npFBA optimum, a few coordinates:")
for name in ("path_glycolysis_share", "path_pentose_phosphate_share", "eff_ppp_per_glycolysis",
             "turnover_nadph_c", "cof_nadph_nadp_share", "conv_atp_per_carbon"):
    print(f"  {name:32s} {phi[name]: .4g}")

states = sample(ecoli, 200, seed=1)               # one state per row, reactions as columns
V = states[[r.id for r in ecoli.reactions]].to_numpy()
table = f2f.flux2feature(model, V, **declaration)
spread = {k: (np.nanpercentile(x, 5), np.nanpercentile(x, 95)) for k, x in table.items()
          if np.isfinite(x).any()}
print(f"\n{len(V)} sampled states; empirical 5-95% range of the sampled states "
      f"(not a confidence interval):")
for name in ("path_glycolysis_share", "path_pentose_phosphate_share"):
    low, high = spread[name]
    print(f"  {name:32s} {low:.3g} to {high:.3g}")
