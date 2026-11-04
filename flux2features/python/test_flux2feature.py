"""Tests of flux2feature.py: by hand, and against the reference state of RBC-GEM.

    python test_flux2feature.py          # or: pytest test_flux2feature.py

1. The self-test of ``flux2feature.py``: a three-reaction model and a toy network
   with a mitochondrion, every value worked out on paper, and the absence rule.
2. The reference check: the first retained state of the healthy reference chain
   of RBC-GEM (``matlab/reference_state.csv``) mapped by this file must give the
   values in ``matlab/flux2feature_reference.csv`` to 1e-5 relative. The notebook
   (stage 28) writes both files after checking the tool against its own engine on
   that state, and ``matlab/test_flux2feature.m`` reads the same two files, so
   the three implementations are held to one set of numbers.

The reference check needs SciPy to read the model; the tool itself needs NumPy only.
"""
import csv
import math
from pathlib import Path

import numpy as np

import flux2feature as f2f

HERE = Path(__file__).resolve().parent
MODEL = HERE.parent / "data" / "RBC-GEM.mat"
STATE = HERE.parent / "matlab" / "reference_state.csv"
REFERENCE = HERE.parent / "matlab" / "flux2feature_reference.csv"


def test_by_hand():
    f2f._self_test()


def _load_model(path: Path) -> "f2f.Model":
    from scipy.io import loadmat
    raw = loadmat(str(path), squeeze_me=True, struct_as_record=False)
    model = next(v for k, v in raw.items() if not k.startswith("__"))
    return f2f.Model(model.S, list(model.rxns), list(model.mets), model.lb, model.ub,
                     list(model.subSystems))


def test_reference_state():
    if not (MODEL.exists() and STATE.exists() and REFERENCE.exists()):
        print("reference files not found; skipped")
        return
    model = _load_model(MODEL)
    with STATE.open() as handle:
        rows = list(csv.DictReader(handle))
    assert [r["reaction"] for r in rows] == model.rxns, "state and model list different reactions"
    state = np.array([float(r["flux"]) for r in rows])
    mine = f2f.flux2feature(model, state, which="all")
    with REFERENCE.open() as handle:
        expected = {r["descriptor"]: float(r["value"]) for r in csv.DictReader(handle)}
    worst = max((abs(mine[k] - v) / max(1.0, abs(v)), k) for k, v in expected.items()
                if not (math.isnan(v) and math.isnan(mine[k])))
    assert worst[0] <= 1e-5, f"{worst[1]} differs by {worst[0]:.2e}"
    print(f"reference: {len(expected)} descriptors agree to 1e-5 (worst {worst[0]:.1e})")


if __name__ == "__main__":
    test_by_hand()
    test_reference_state()
