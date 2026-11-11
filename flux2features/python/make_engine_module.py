"""Write python/flux2features_engine.py from the analysis notebook.

The notebook is the one place the descriptor engine is written: the
declarations of the study (pathway groups, named reactions, tolerances), the
model reader, the reaction annotation, every descriptor family and the engine
that runs them. This script copies those cells, unchanged and in notebook
order, into a module that can be imported without the notebook, and adds a
small interface (``phi``, ``declared``). It never edits the engine: to change
a descriptor, change the notebook and run this script again.

    python python/make_engine_module.py

Stage 28 of the notebook checks that the module and the notebook's own engine
return the same 407 values on the reference state of RBC-GEM.
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTEBOOK = HERE.parent / "notebook" / "metabolic_flux_phenotyping.ipynb"
TARGET = HERE / "flux2features_engine.py"

#: A cell is part of the engine if it defines one of these. Found by content, not
#: by position, so that adding a cell to the notebook does not break the export.
ENGINE_MARKERS = (
    "TOL = Tolerances()",            # the declarations of the study
    "class GEM:",                    # the model reader
    "def annotate_reactions(",       # the reaction annotation
    "class FeatureSpec:",            # the registry and the shared helpers
    "class StructuralFamily",        # the generic families
    "class PathwayFamily",           # pathway, turnover, named reactions, network
    "class ExternalFamily",          # the boundary
    "class CofactorFamily",          # cofactors, capacity, coupling
    "class BudgetFamily",            # budgets
    "class CompartmentFamily",       # compartments
    "class ElementalFamily",         # elements, amino acids, sources, defence, allocation
    "class FeatureEngine:",          # the engine
)

HEADER = '''"""flux2features engine: every declared descriptor of Phi_M(v).

GENERATED from notebook/metabolic_flux_phenotyping.ipynb by
python/make_engine_module.py. Do not edit: change the notebook and regenerate.

This is the engine the study ran: it computes all {declared} declared
descriptors, each NaN when the reconstruction lacks what it requires, with
NumPy, SciPy and pandas. The standalone tool (flux2feature.py) is a smaller,
dependency-light subset of it with a MATLAB twin.

    import flux2features_engine as fe
    gem = fe.load_gem("RBC-GEM.mat")          # a COBRA model saved by MATLAB
    table = fe.phi(gem, V)                    # V: states x reactions
    fe.declared(gem)                          # every descriptor, defined or not, and why

The declarations of the study (PATHWAY_GROUPS, TURNOVER_METABOLITES,
KEY_REACTIONS, EXCHANGE_ADMISSIBILITY, ...) describe RBC-GEM. For another
reconstruction, replace them before calling phi: a descriptor whose
requirement the model does not meet is reported as NaN, never as zero.
"""
'''

INTERFACE = '''

# --------------------------------------------------------------------------
# Interface
# --------------------------------------------------------------------------

def engine(gem: "GEM") -> "FeatureEngine":
    """The descriptor engine for one reconstruction."""
    return FeatureEngine(gem, annotate_reactions(gem))


def declared(gem: "GEM") -> pd.DataFrame:
    """Every declared descriptor, with whether this reconstruction defines it and why not."""
    built = engine(gem)
    table = built.registry_frame()
    table["defined"] = ~table["feature"].isin(built.unavailable)
    table["why_not"] = table["feature"].map(built.unavailable).fillna("")
    return table


def phi(gem: "GEM", V, *, built: "FeatureEngine | None" = None) -> pd.DataFrame:
    """Phi_M(v) of every state in V (states x reactions): every declared column,
    NaN where the reconstruction cannot define the descriptor."""
    built = built or engine(gem)
    V = np.atleast_2d(np.asarray(V, dtype=float))
    if V.shape[1] != len(gem.rxns) and V.shape[0] == len(gem.rxns):
        V = V.T
    values = built.transform(V)
    columns = built.registry_frame()["feature"].tolist()
    return values.reindex(columns=columns)
'''


def main() -> None:
    cells = json.loads(NOTEBOOK.read_text())["cells"]
    chosen, found = [], set()
    for cell in cells:
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        hits = [m for m in ENGINE_MARKERS if m in source]
        if hits:
            chosen.append(source.rstrip() + "\n")
            found.update(hits)
    missing = [m for m in ENGINE_MARKERS if m not in found]
    if missing:
        raise SystemExit(f"cells not found in the notebook: {missing}")
    body = "\n\n".join(chosen)
    declared_count = "every"
    TARGET.write_text(HEADER.format(declared=declared_count) + "\nimport pandas as pd\n\n"
                      + body + INTERFACE)
    print(f"{len(chosen)} cells written to {TARGET}")


if __name__ == "__main__":
    main()
