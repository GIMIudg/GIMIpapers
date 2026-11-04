"""Independent validation of the Python flux2feature, written to be rerun by anyone.

    python validation/validate_python.py          (from the flux2features folder)

What it checks, each reported PASS or FAIL in validation/results/python_report.md:

  1. The unit tests of python/test_flux2feature.py (values worked by hand, the
     absence rule, and the reference state of RBC-GEM against reference values).
  2. Every case in validation/cases/ and two RBC-GEM states (the reference state
     and the zero vector) mapped with which="all"; the values are written to
     validation/results/python/<case>.csv for the cross-check with MATLAB.
  3. Properties any correct implementation must have, on the RBC-GEM state:
       - one state as a vector, as a row and as a column of a matrix give the same values;
       - a matrix of states gives, row by row, what each state gives alone;
       - reordering the reactions and metabolites of the model (and the flux
         vector with them) changes no descriptor;
       - doubling every flux doubles every descriptor in flux units and leaves
         every ratio, share, entropy and statistic of shape unchanged (counts and
         capacity use, which compare fluxes with a tolerance or with bounds, are
         excluded and listed);
       - every declared descriptor appears with which="all", the unavailable ones
         as NaN with a reason, never 0;
       - duplicate identifiers, non-finite fluxes and a vector of the wrong
         length are refused with an error.
  4. Identifier conventions: the two-compartment network written with BiGG
     identifiers (atp_m) and with VMH identifiers (atp[m]) gives the same values.

It changes nothing in the tool. A FAIL is a finding to document, not to fix here.
"""
import csv
import math
import platform
import sys
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "python"))
import flux2feature as f2f            # noqa: E402
import test_flux2feature as unit      # noqa: E402

OUT = HERE / "results"
CASES = HERE / "cases"
RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(condition), detail))
    print(f"  {'PASS' if condition else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))


def load_case(folder: Path):
    with (folder / "reactions.csv").open() as h:
        rows = list(csv.DictReader(h))
    rxns = [r["reaction"] for r in rows]
    with (folder / "stoichiometry.csv").open() as h:
        triplets = list(csv.DictReader(h))
    mets = list(dict.fromkeys(t["metabolite"] for t in triplets))
    S = np.zeros((len(mets), len(rxns)))
    for t in triplets:
        S[mets.index(t["metabolite"]), rxns.index(t["reaction"])] = float(t["coefficient"])
    model = f2f.Model(S, rxns, mets, [float(r["lb"]) for r in rows], [float(r["ub"]) for r in rows],
                      [r["subsystem"] for r in rows])
    with (folder / "flux.csv").open() as h:
        flux = {r["reaction"]: float(r["flux"]) for r in csv.DictReader(h)}
    return model, np.array([flux[r] for r in rxns])


def load_rbc():
    model = unit._load_model(unit.MODEL)
    with unit.STATE.open() as h:
        state = np.array([float(r["flux"]) for r in csv.DictReader(h)])
    return model, state


def write_values(name: str, model, v) -> dict:
    values = f2f.flux2feature(model, v, which="all")
    audit = {r["name"]: r for r in f2f.audit(model)}
    (OUT / "python").mkdir(parents=True, exist_ok=True)
    with (OUT / "python" / f"{name}.csv").open("w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["descriptor", "value", "computable"])
        for k, x in values.items():
            w.writerow([k, "NaN" if math.isnan(x) else repr(float(x)),
                        int(bool(audit[k]["computable"]))])
    return values


def same(a: dict, b: dict, rtol=1e-9) -> tuple[bool, str]:
    worst, where = 0.0, ""
    for k in a:
        x, y = a[k], b[k]
        if math.isnan(x) and math.isnan(y):
            continue
        if math.isnan(x) != math.isnan(y):
            return False, f"{k}: {x} vs {y}"
        d = abs(x - y) / max(1.0, abs(y))
        if d > worst:
            worst, where = d, k
    return worst <= rtol, f"worst {worst:.1e} at {where}" if where else "identical"


def refuses(call) -> bool:
    try:
        call()
    except (ValueError, KeyError):
        return True
    return False


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Python {platform.python_version()}, NumPy {np.__version__}, {platform.platform()}")

    print("\n1. unit tests")
    for test in (unit.test_by_hand, unit.test_reference_state):
        try:
            test()
            check(test.__name__, True)
        except Exception as error:                        # noqa: BLE001
            check(test.__name__, False, f"{type(error).__name__}: {error}")
            traceback.print_exc()

    print("\n2. cases written for the cross-check")
    for folder in sorted(p for p in CASES.iterdir() if p.is_dir()):
        model, v = load_case(folder)
        write_values(folder.name, model, v)
        check(f"case {folder.name} mapped", True)
    model, state = load_rbc()
    reference = write_values("rbc_reference_state", model, state)
    write_values("rbc_zero_vector", model, np.zeros_like(state))
    check("RBC-GEM reference and zero vector mapped", True)

    print("\n3. properties, on the RBC-GEM state")
    alone = f2f.flux2feature(model, state, which="all")
    row = {k: float(x[0]) for k, x in f2f.flux2feature(model, state[None, :], which="all").items()}
    column = {k: float(x[0]) for k, x in f2f.flux2feature(model, state[:, None], which="all").items()}
    check("vector, row and column give the same values",
          same(alone, row)[0] and same(alone, column)[0])
    batch = f2f.flux2feature(model, np.vstack([state, 0.5 * state, state]), which="all")
    second = f2f.flux2feature(model, 0.5 * state, which="all")
    check("a matrix gives what each state gives alone",
          same(alone, {k: float(x[0]) for k, x in batch.items()})[0]
          and same(second, {k: float(x[1]) for k, x in batch.items()})[0])

    rng = np.random.default_rng(0)
    pr, pm = rng.permutation(len(model.rxns)), rng.permutation(len(model.mets))
    shuffled = f2f.Model(model.S[np.ix_(pm, pr)], [model.rxns[i] for i in pr],
                         [model.mets[i] for i in pm], model.lb[pr], model.ub[pr],
                         [model.subsystems[i] for i in pr])
    ok, detail = same(alone, f2f.flux2feature(shuffled, state[pr], which="all"), rtol=1e-9)
    check("reordering reactions and metabolites changes nothing", ok, detail)

    declared = {r["name"]: r for r in f2f.audit(model)}
    doubled = f2f.flux2feature(model, 2 * state, which="all")
    excluded = {k for k, r in declared.items()
                if r["units"] == "count" or r["family"] == "capacity"}
    bad = []
    for k, r in declared.items():
        if k in excluded or not r["computable"] or math.isnan(alone[k]):
            continue
        expected = 2 * alone[k] if r["units"] == "flux" else alone[k]
        if abs(doubled[k] - expected) > 1e-9 * max(1.0, abs(expected)):
            bad.append(k)
    check("doubling fluxes scales flux descriptors and leaves the rest unchanged", not bad,
          f"excluded by design: {len(excluded)}; failing: {bad[:5]}" if bad
          else f"excluded by design: {len(excluded)}")

    unavailable = [k for k, r in declared.items() if not r["computable"]]
    check(f"all {len(declared)} declared descriptors returned with which='all'",
          set(alone) == set(declared))
    check(f"{len(unavailable)} unavailable descriptors are NaN with a reason, never 0",
          all(math.isnan(alone[k]) and declared[k]["reason"] for k in unavailable))
    check("every computable descriptor is finite on the reference state",
          sum(1 for k in declared if declared[k]["computable"] and math.isnan(alone[k])) == 0,
          f"{sum(1 for k in declared if declared[k]['computable'] and math.isnan(alone[k]))} NaN")

    check("duplicate reaction identifiers are refused",
          refuses(lambda: f2f.Model(np.hstack([model.S, model.S[:, :1]]),
                                    model.rxns + model.rxns[:1], model.mets)))
    check("duplicate metabolite identifiers are refused",
          refuses(lambda: f2f.Model(np.vstack([model.S, model.S[:1]]),
                                    model.rxns, model.mets + model.mets[:1])))
    bad_state = state.copy()
    bad_state[0] = np.nan
    check("non-finite fluxes are refused", refuses(lambda: f2f.flux2feature(model, bad_state)))
    check("a vector of the wrong length is refused",
          refuses(lambda: f2f.flux2feature(model, state[:-1])))

    print("\n4. identifier conventions")
    a = f2f.flux2feature(*load_case(CASES / "compartment_bigg"), which="all")
    vmh_model, vmh_flux = load_case(CASES / "compartment_vmh")
    b = f2f.flux2feature(vmh_model, vmh_flux, which="all")
    compartment = [k for k in a if k.startswith("cmp_")]
    ok, detail = same({k: a[k] for k in compartment}, {k: b[k] for k in compartment})
    check("compartment descriptors: atp_m (BiGG) and atp[m] (VMH) give the same values",
          ok, detail)
    # Descriptors that name BiGG identifiers (turnover_atp_c, rxn_...) cannot be resolved in
    # VMH identifiers without a mapping: they must be declared unavailable, never wrong. The
    # network descriptors exclude a currency list written in BiGG identifiers; given the same
    # list in VMH identifiers (currency=...), they must agree.
    import re
    currency = [re.sub(r"_([a-z])$", lambda m: f"[{m.group(1)}]", c) for c in f2f.CURRENCY]
    b = f2f.flux2feature(vmh_model, vmh_flux, which="all", currency=currency)
    vmh_audit = {r["name"]: r for r in f2f.audit(vmh_model, currency=currency)}
    differ = [k for k in a if not k.startswith("cmp_")
              and not (math.isnan(a[k]) and math.isnan(b[k]))
              and not math.isclose(a[k], b[k], rel_tol=1e-9, abs_tol=1e-12)]
    honest = all(math.isnan(b[k]) and not vmh_audit[k]["computable"] and vmh_audit[k]["reason"]
                 for k in differ)
    check("with the currency list mapped, the rest agree or are declared unavailable "
          "(named BiGG identifiers need a mapping)",
          honest, f"{len(differ)} unavailable without a mapping, e.g. {differ[:3]}")
    for name, expected in f2f.COMPARTMENT_TOY_EXPECTED.items():
        if not math.isclose(b[name], expected, rel_tol=1e-9):
            check(f"VMH network: {name}", False, f"{b[name]} != {expected}")
            break
    else:
        check("VMH network: the 11 compartment values worked by hand", True)

    passed = sum(ok for _, ok, _ in RESULTS)
    lines = ["# Python validation of flux2feature", "",
             f"Python {platform.python_version()}, NumPy {np.__version__}, {platform.platform()}", "",
             f"{passed} of {len(RESULTS)} checks passed.", "",
             "| check | result | detail |", "|---|---|---|"]
    lines += [f"| {n} | {'PASS' if ok else 'FAIL'} | {d} |" for n, ok, d in RESULTS]
    (OUT / "python_report.md").write_text("\n".join(lines) + "\n")
    print(f"\n{passed} of {len(RESULTS)} checks passed; report in {OUT / 'python_report.md'}")
    return 0 if passed == len(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
