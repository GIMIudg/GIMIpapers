"""flux2feature in five minutes, on a network small enough to read.

    python examples/quickstart.py

A toy erythrocyte: glucose enters, goes down glycolysis to lactate, and a tenth
of it enters the oxidative pentose phosphate pathway, whose NADPH keeps
glutathione reduced. The flux vector below is a steady state (S v = 0).
"""
import numpy as np

import flux2feature as f2f

# ---------------------------------------------------------------- the model
mets = ["glc__D_e", "glc__D_c", "g6p_c", "6pgc_c", "ru5p__D_c", "pyr_c", "lac__L_c",
        "lac__L_e", "nadph_c", "nadp_c", "co2_c", "co2_e", "atp_c", "adp_c",
        "gthrd_c", "gthox_c"]
reactions = {  # identifier: ({metabolite: coefficient}, subsystem, lb, ub)
    "EX_glc__D_e": ({"glc__D_e": -1}, "", -1, 0),
    "GLCt1": ({"glc__D_e": -1, "glc__D_c": 1}, "Transport, extracellular", 0, 10),
    "HEX1": ({"glc__D_c": -1, "atp_c": -1, "g6p_c": 1, "adp_c": 1},
             "Glycolysis / Gluconeogenesis", 0, 10),
    "GLYC": ({"g6p_c": -1, "adp_c": -2, "pyr_c": 2, "atp_c": 2},          # lumped
             "Glycolysis / Gluconeogenesis", 0, 10),
    "LDH_L": ({"pyr_c": -1, "lac__L_c": 1}, "Glycolysis / Gluconeogenesis", -10, 10),
    "L_LACt": ({"lac__L_c": -1, "lac__L_e": 1}, "Transport, extracellular", -10, 10),
    "EX_lac__L_e": ({"lac__L_e": -1}, "", 0, 10),
    "G6PDH2": ({"g6p_c": -1, "nadp_c": -1, "6pgc_c": 1, "nadph_c": 1},   # with PGL
               "Pentose phosphate pathway", 0, 10),
    "GND": ({"6pgc_c": -1, "nadp_c": -1, "ru5p__D_c": 1, "nadph_c": 1, "co2_c": 1},
            "Pentose phosphate pathway", 0, 10),
    "PPPret": ({"ru5p__D_c": -1, "g6p_c": 1}, "Pentose phosphate pathway", 0, 10),
    "GTHOr": ({"gthox_c": -1, "nadph_c": -1, "gthrd_c": 2, "nadp_c": 1},
              "Glutathione metabolism", 0, 10),
    "GTHP": ({"gthrd_c": -2, "gthox_c": 1}, "Glutathione metabolism", 0, 10),
    "CO2t": ({"co2_c": -1, "co2_e": 1}, "Transport, extracellular", -10, 10),
    "EX_co2_e": ({"co2_e": -1}, "", 0, 10),
    "ATPM": ({"atp_c": -1, "adp_c": 1}, "", 0, 10),
}
rxns = list(reactions)
S = np.zeros((len(mets), len(rxns)))
for j, (coefficients, *_ ) in enumerate(reactions.values()):
    for met, c in coefficients.items():
        S[mets.index(met), j] = c
model = f2f.Model(S, rxns, mets,
                  lb=[r[2] for r in reactions.values()],
                  ub=[r[3] for r in reactions.values()],
                  subsystems=[r[1] for r in reactions.values()])
# From COBRApy instead:  model = f2f.Model.from_cobra(cobra_model)

v = np.array([-1, 1, 1, 1, 2, 2, 2, 0.1, 0.1, 0.1, 0.2, 0.2, 0.1, 0.1, 1])
assert np.allclose(S @ v, 0), "the example state is a steady state"

# ---------------------------------------------------------------- 1. audit first
rows = f2f.audit(model)
print(f"{len(rows)} descriptors declared, {sum(r['computable'] for r in rows)} computable here")
for tier in f2f.TIERS:
    mine = [r for r in rows if r["tier"] == tier]
    print(f"  {tier:8s} {sum(r['computable'] for r in mine):3d} of {len(mine)} computable")
why = {r["name"]: r["reason"] for r in rows}
print("rxn_pgl is not computable here:", why["rxn_pgl"])

# ---------------------------------------------------------------- 2. the core of Phi
core = f2f.flux2feature(model, v)            # which="core" is the default
finite = [k for k, x in core.items() if np.isfinite(x)]
print(f"\nthe core: {len(core)} coordinates, {len(finite)} finite on this model")
for name in ("path_pentose_phosphate_share", "turnover_nadph_c", "eff_nadph_per_glucose",
             "rxn_gpx", "conv_lactate_per_glucose_carbon"):
    print(f"  {name:34s} {core[name]: .4g}")
# conv_* need metabolite formulas (Model(..., formulas=...)), which this toy
# does not give, and the compartment shares need a second compartment; both are
# NaN here, with the reason in the audit:
print(f"  {'':34s} ({why['conv_lactate_per_glucose_carbon']})")
print(f"  {'cmp_atp_made_in_mitochondrion':34s} {core['cmp_atp_made_in_mitochondrion']: .4g}")
assert np.isnan(core["cmp_atp_made_in_mitochondrion"]), "absent is NaN, never 0"

# ---------------------------------------------------------------- 3. more than the core
every = f2f.flux2feature(model, v, which="all")         # the 179 declared descriptors
named = f2f.flux2feature(model, v, which=["rxn_g6pd", "rxn_pgl", "int_L1"])
print(f"\nwhich='all': {len(every)} coordinates; by name: {named}")

V = np.vstack([v, 2 * v])                    # one state per row (or per column)
many = f2f.flux2feature(model, V, which=("core", "module"))
print(f"{len(many)} coordinates for {V.shape[0]} states; "
      f"turnover_nadph_c = {many['turnover_nadph_c']}")

# ---------------------------------------------------------------- 4. another naming scheme
# Had this model called G6PD "G6PDH2r" (as BiGG's E. coli core does), map the
# label to the identifier; the coordinate keeps its name, rxn_g6pd.
renamed = f2f.Model(S, [r if r != "G6PDH2" else "G6PDH2r" for r in rxns], mets,
                    subsystems=[r[1] for r in reactions.values()])
mapping = {"g6pd": "G6PDH2r"}                # merged with the defaults; nothing is lost
print("without the mapping:", f2f.flux2feature(renamed, v, which=["rxn_g6pd"]))
print("with the mapping:   ",
      f2f.flux2feature(renamed, v, which=["rxn_g6pd"], named_reactions=mapping))
