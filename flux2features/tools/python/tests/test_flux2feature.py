"""Tests of the flux2feature package: by hand, and against the reference state of RBC-GEM.

    pytest tests/                        # from tools/python, after: pip install -e .

1. The self-test of the package: a three-reaction model, a toy network
   with a mitochondrion and a toy exchanging glucose, glutamine, lactate and CO2
   with peroxide made and removed two ways, every value worked out on paper,
   the absence rule, and the tiers of the release (the default output is the
   core; generic statistics and aliases come only on request).
2. The reference check: the first retained state of the healthy reference chain
   of RBC-GEM (``tests/data/reference_state.csv``) mapped by the package must give the
   values in ``tests/data/flux2feature_reference.csv`` to 1e-5 relative. The notebook
   (stage 28) writes both files after checking the tool against its own engine on
   that state, and the MATLAB test (``testFlux2feature.m``) reads the same two files, so
   the three implementations are held to one set of numbers.

The reference check needs SciPy to read the model (``tests/data/RBC-GEM_flux2feature.mat``,
RBC-GEM reduced to the fields the tool reads); the tool itself needs NumPy only.
"""
import csv
import math
import os
from pathlib import Path

import numpy as np

import flux2feature as f2f

try:                       # the file also runs as a plain script, without pytest
    import pytest
except ImportError:        # pragma: no cover
    pytest = None

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
#: RBC-GEM reduced to the fields the tool reads (written by the notebook's
#: stage 28); FLUX2FEATURE_RBCGEM may point at the full model instead.
MODEL = Path(os.environ.get("FLUX2FEATURE_RBCGEM", DATA / "RBC-GEM_flux2feature.mat"))
STATE = DATA / "reference_state.csv"
REFERENCE = DATA / "flux2feature_reference.csv"


def test_by_hand():
    f2f._self_test()


def _load_model(path: Path) -> "f2f.Model":
    from scipy.io import loadmat
    raw = loadmat(str(path), squeeze_me=True, struct_as_record=False)
    model = next(v for k, v in raw.items() if not k.startswith("__"))
    return f2f.Model(model.S, list(model.rxns), list(model.mets), model.lb, model.ub,
                     list(model.subSystems),
                     formulas=list(getattr(model, "metFormulas", [""] * len(model.mets))),
                     met_names=list(getattr(model, "metNames", [""] * len(model.mets))))


def test_reference_state():
    if not (MODEL.exists() and STATE.exists() and REFERENCE.exists()):
        if pytest is not None and __name__ != "__main__":
            pytest.skip("reference files not found")
        print("reference files not found; skipped")
        return
    model = _load_model(MODEL)
    with STATE.open() as handle:
        rows = list(csv.DictReader(handle))
    assert [r["reaction"] for r in rows] == model.rxns, "state and model list different reactions"
    state = np.array([float(r["flux"]) for r in rows])
    mine = f2f.flux2feature(model, state, which="all")
    with REFERENCE.open() as handle:
        expected = {r["descriptor"]: float(r["value"]) if r["value"] else math.nan
                    for r in csv.DictReader(handle)}
    worst = max((abs(mine[k] - v) / max(1.0, abs(v)), k) for k, v in expected.items()
                if not (math.isnan(v) and math.isnan(mine[k])))
    assert worst[0] <= 1e-5, f"{worst[1]} differs by {worst[0]:.2e}"
    print(f"reference: {len(expected)} descriptors agree to 1e-5 (worst {worst[0]:.1e})")


def test_release_tiers():
    """Every declared descriptor has a tier, and every alias names a canonical form."""
    assert f2f.STATUS, "the package carries no descriptor status; run the notebook's stage 28"
    declared = {r["name"] for r in f2f.audit(_toy())}
    assert declared == set(f2f.STATUS), "the status block and the registry disagree"
    assert {t for t, _ in f2f.STATUS.values()} <= set(f2f.TIERS)
    for name, (tier, canonical) in f2f.STATUS.items():
        assert (tier == "alias") == bool(canonical), name
        if canonical:
            assert f2f.STATUS[canonical][0] != "alias", f"{name} points to another alias"
    counts = {t: sum(v[0] == t for v in f2f.STATUS.values()) for t in f2f.TIERS}
    print("tiers: " + ", ".join(f"{n} {t}" for t, n in counts.items()))


def _toy() -> "f2f.Model":
    S = np.array([[1.0, -1.0, 0.0], [0.0, 1.0, -1.0]])
    return f2f.Model(S, ["EX_a", "R1", "EX_b"], ["a_c", "b_c"], lb=[-10, 0, 0], ub=[0, 10, 10],
                     subsystems=["", "Glycolysis / Gluconeogenesis", ""])


if __name__ == "__main__":
    test_by_hand()
    test_release_tiers()
    test_reference_state()
