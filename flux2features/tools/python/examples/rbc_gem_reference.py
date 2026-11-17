"""flux2feature on a genome-scale model: a sampled state of RBC-GEM.

    python examples/rbc_gem_reference.py            # needs SciPy to read the .mat file

Reads RBC-GEM (Haiman et al., 2025; reduced to the fields flux2feature uses,
in tests/data) and the first retained state of the healthy reference chain of
the accompanying study, maps the state to the core of Phi, and checks it
against the reference values the notebook (stage 28) and the MATLAB tool share.
Set FLUX2FEATURE_RBCGEM to the full RBC-GEM.mat to use it instead.
"""
import csv
import math
import os
from pathlib import Path

import numpy as np
from scipy.io import loadmat

import flux2feature as f2f

DATA = Path(__file__).resolve().parents[1] / "tests" / "data"
MODEL = Path(os.environ.get("FLUX2FEATURE_RBCGEM", DATA / "RBC-GEM_flux2feature.mat"))

# a COBRA model structure saved by MATLAB or by the COBRA Toolbox
raw = loadmat(MODEL, squeeze_me=True, struct_as_record=False)
cobra = next(value for key, value in raw.items() if not key.startswith("__"))
model = f2f.Model(cobra.S, list(cobra.rxns), list(cobra.mets), cobra.lb, cobra.ub,
                  list(cobra.subSystems), formulas=list(cobra.metFormulas),
                  met_names=list(cobra.metNames))
print(f"RBC-GEM: {len(model.mets)} metabolites, {len(model.rxns)} reactions, "
      f"{int(model.exchange.sum())} boundary reactions")

rows = f2f.audit(model)
print(f"{sum(r['computable'] for r in rows)} of {len(rows)} declared descriptors computable")
for r in rows:
    if r["tier"] == "core" and not r["computable"]:
        print(f"  core, not computable: {r['name']} ({r['reason']})")

with (DATA / "reference_state.csv").open() as handle:
    state = {row["reaction"]: float(row["flux"]) for row in csv.DictReader(handle)}
v = np.array([state[r] for r in model.rxns])

core = f2f.flux2feature(model, v)
print(f"\nthe core of Phi on this state ({len(core)} coordinates), a few of them:")
for name in ("path_glycolysis_share", "path_pentose_phosphate_share", "turnover_atp_c",
             "turnover_nadph_c", "cof_glutathione_share", "cpl_boundary_share"):
    print(f"  {name:42s} {core[name]: .5g}")

with (DATA / "flux2feature_reference.csv").open() as handle:
    reference = {row["descriptor"]: float(row["value"]) if row["value"] else math.nan
                 for row in csv.DictReader(handle)}
every = f2f.flux2feature(model, v, which="all")
worst = max(abs(every[k] - x) / max(1.0, abs(x)) for k, x in reference.items()
            if not (math.isnan(x) and math.isnan(every[k])))
print(f"\nall {len(reference)} descriptors match the shared reference values "
      f"(worst relative difference {worst:.1e})")

# Many states: rows (or columns, as COBRA samplers return them) of a matrix.
# With pandas installed, as_frame=True gives one row per state:
#   table = f2f.flux2feature(model, samples, which=("core", "module"), as_frame=True)
