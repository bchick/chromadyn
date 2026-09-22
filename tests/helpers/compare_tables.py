#!/usr/bin/env python3
"""Compare two TSVs as data, not as bytes.

    python tests/helpers/compare_tables.py golden.tsv observed.tsv [rtol]

Headers and non-numeric cells must match exactly. Numeric cells must agree to
a relative tolerance (default 1e-6). A byte-for-byte diff is wrong for tables
of floats: the same seeded run on a different CPU or BLAS build differs in the
eleventh significant digit, which is how tier 2 came to pass locally and fail
on every CI runner while the results were identical.
"""
import csv
import math
import sys


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except ValueError:
        return None


def main(a, b, rtol=1e-6, atol=1e-9):
    ra, rb = rows(a), rows(b)
    if len(ra) != len(rb):
        print(f"row count differs: {len(ra)} against {len(rb)}", file=sys.stderr)
        return 1
    bad = 0
    for i, (x, y) in enumerate(zip(ra, rb), start=1):
        if len(x) != len(y):
            print(f"line {i}: column count differs", file=sys.stderr); bad += 1; continue
        for j, (u, v) in enumerate(zip(x, y), start=1):
            fu, fv = num(u), num(v)
            same = (u == v) if (fu is None or fv is None or i == 1) else \
                math.isclose(fu, fv, rel_tol=rtol, abs_tol=atol)
            if not same:
                bad += 1
                if bad <= 10:
                    print(f"line {i} col {j}: {u!r} != {v!r}", file=sys.stderr)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2], *(float(t) for t in sys.argv[3:4])))
