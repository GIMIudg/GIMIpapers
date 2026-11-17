"""Tests of the released engine and of the release of Phi.

    pytest test_engine_release.py        # from python/, or: python test_engine_release.py

1. The release table (``descriptor_release.csv``, written by the notebook's
   stage 28 from the audit of stage 27c) is complete and consistent: every
   declared descriptor appears once, in one of the four tiers; every alias names
   a canonical descriptor that exists and is not itself an alias; a module names
   its module.
2. The standalone tool (``tools/python``) carries the same tiers as the release for every
   descriptor the two share; an alias whose canonical form the tool does not
   implement stands in for it and takes its tier.
3. The released engine (``flux2features_engine.py``) maps the reference state of
   RBC-GEM: it returns every declared coordinate, NaN (never 0) for the
   compartment family that a model with one compartment cannot define, and the
   core, modules, generic statistics and aliases as the release says.

The engine needs NumPy, SciPy and pandas; the tests skip the third part if the
repository's model or reference state is absent.
"""
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))                                        # the released engine
sys.path.insert(0, str(HERE.parent / "tools" / "python" / "src"))    # the standalone tool
import flux2feature as f2f  # noqa: E402

RELEASE = HERE / "descriptor_release.csv"
MODEL = HERE.parent / "data" / "RBC-GEM.mat"
STATE = HERE.parent / "tools" / "python" / "tests" / "data" / "reference_state.csv"
TIERS = ("core", "module", "generic", "alias")


def _release() -> pd.DataFrame:
    return pd.read_csv(RELEASE, keep_default_na=False)


def test_release_is_complete_and_consistent():
    table = _release()
    assert table["feature"].is_unique
    assert set(table["release"]) <= set(TIERS)
    assert all((table["release"] == tier).any() for tier in TIERS)
    tier = table.set_index("feature")["release"]
    aliases = table[table["release"] == "alias"]
    assert (aliases["canonical"] != "").all(), "every alias names its canonical form"
    assert aliases["canonical"].isin(tier.index).all(), "every canonical form is declared"
    assert (tier.reindex(aliases["canonical"]) != "alias").all(), "a canonical form is not an alias"
    assert (table.loc[table["release"] != "alias", "canonical"] == "").all()
    modules = table[table["release"] == "module"]
    assert (modules["module"] != "").all(), "every module descriptor names its module"
    assert (table["evidence"] != "").all(), "every tier carries its evidence"


def test_tool_agrees_with_the_release():
    release = _release().set_index("feature")
    tier, canonical = release["release"], release["canonical"]
    shared = [name for name in f2f.STATUS if name in tier.index]
    assert len(shared) > 100

    def expected(name):
        # an alias whose canonical form the tool does not declare stands in for it
        if tier[name] == "alias" and canonical[name] not in f2f.STATUS:
            return tier[canonical[name]]
        return tier[name]
    disagree = [name for name in shared if f2f.STATUS[name][0] != expected(name)]
    assert not disagree, f"tiers differ between the tool and the release: {disagree[:10]}"


@pytest.mark.skipif(not (MODEL.exists() and STATE.exists()),
                    reason="the repository's model or reference state is absent")
def test_engine_on_the_reference_state():
    import flux2features_engine as fe
    gem = fe.load_gem(MODEL)
    with open(STATE, newline="") as handle:
        rows = list(csv.reader(handle))[1:]
    flux = dict((reaction, float(value)) for reaction, value in rows)
    v = np.array([flux[r] for r in gem.rxns], dtype=float)
    built = fe.engine(gem)
    every = fe.phi(gem, v, built=built)
    release = _release()
    assert list(every.columns) == fe.declared(gem)["feature"].tolist()
    assert len(every.columns) == len(release)
    compartments = [c for c in every.columns if c.startswith("cmp_")]
    assert compartments and every[compartments].isna().all(axis=None), \
        "a model with one compartment gives NaN, not 0, for the compartment family"
    for tier in TIERS:
        part = fe.phi(gem, v, built=built, tiers=tier)
        assert set(part.columns) == set(release.loc[release["release"] == tier, "feature"])
        assert np.allclose(part.to_numpy(), every[part.columns].to_numpy(), equal_nan=True)
    with pytest.raises(ValueError):
        fe.phi(gem, v, built=built, tiers="everything")


if __name__ == "__main__":
    test_release_is_complete_and_consistent()
    test_tool_agrees_with_the_release()
    if MODEL.exists() and STATE.exists():
        test_engine_on_the_reference_state()
    print("all release checks passed")
