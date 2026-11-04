"""Python against MATLAB, descriptor by descriptor, on the same inputs.

    python validation/cross_check.py      (after validate_python.py and validate_matlab.m)

For every case written by both validators (validation/results/python/<case>.csv
and validation/results/matlab/<case>.csv) it checks that the two implementations
declare the same descriptors, agree on which ones the model can compute, and give
the same values: NaN where the other gives NaN, and otherwise within 1e-9
relative (both read the same inputs; only the order of floating-point sums
differs). The report goes to validation/results/cross_check_report.md.
"""
import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
RTOL = 1e-9


def read(path: Path) -> dict:
    with path.open() as handle:
        return {r["descriptor"]: (float(r["value"]), r["computable"] in ("1", "True", "true"))
                for r in csv.DictReader(handle)}


def compare(python: dict, matlab: dict) -> tuple[bool, str]:
    if set(python) != set(matlab):
        only_p, only_m = sorted(set(python) - set(matlab)), sorted(set(matlab) - set(python))
        return False, f"only in Python: {only_p[:4]}; only in MATLAB: {only_m[:4]}"
    flags = [k for k in python if python[k][1] != matlab[k][1]]
    if flags:
        return False, f"computable differs: {flags[:4]}"
    worst, where = 0.0, ""
    for k, (x, _) in python.items():
        y = matlab[k][0]
        if math.isnan(x) or math.isnan(y):
            if math.isnan(x) != math.isnan(y):
                return False, f"{k}: Python {x}, MATLAB {y}"
            continue
        d = abs(x - y) / max(1.0, abs(x))
        if d > worst:
            worst, where = d, k
    ok = worst <= RTOL
    return ok, f"{len(python)} descriptors; worst relative difference {worst:.1e}" + (
        f" ({where})" if where else "")


def main() -> int:
    cases = sorted(p.stem for p in (RESULTS / "python").glob("*.csv"))
    rows = []
    for case in cases:
        matlab = RESULTS / "matlab" / f"{case}.csv"
        if not matlab.exists():
            rows.append((case, False, "no MATLAB result: run validate_matlab.m first"))
            continue
        ok, detail = compare(read(RESULTS / "python" / f"{case}.csv"), read(matlab))
        rows.append((case, ok, detail))
    for case, ok, detail in rows:
        print(f"  {'PASS' if ok else 'FAIL'}  {case}  ({detail})")
    passed = sum(ok for _, ok, _ in rows)
    lines = ["# Python against MATLAB", "", f"{passed} of {len(rows)} cases agree "
             f"(relative tolerance {RTOL:g}).", "", "| case | result | detail |", "|---|---|---|"]
    lines += [f"| {c} | {'PASS' if ok else 'FAIL'} | {d} |" for c, ok, d in rows]
    (RESULTS / "cross_check_report.md").write_text("\n".join(lines) + "\n")
    print(f"{passed} of {len(rows)} cases agree; report in {RESULTS / 'cross_check_report.md'}")
    return 0 if rows and passed == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
