"""Score the curated real-reaction dataset and reproduce the manuscript tables.

Reads ``curated_reactions.csv`` (one row per curated kinetic--thermodynamic
comparison, with the kinetic estimate ``keq_kin`` and the thermodynamic
comparator ``keq_thermo`` already reduced to apparent equilibrium constants,
together with full provenance) and reports, for each record,

    x   = keq_kin / keq_thermo
    C   = J(x) = 1/2 (x + 1/x) - 1          (reciprocal recognition cost)
    fold= max(x, 1/x)
    fold band                                (fixed reporting band)

These reproduce the phosphoglucose isomerase rows, the clean scored set (TPI
and the racemase controls), and the fumarase kinetic Haldane scores reported in
the manuscript. Standard library only.

Usage
-----
    python score_curated.py                       # print summary
    python score_curated.py -o scored.csv         # also write a scored CSV
"""
from __future__ import annotations

import argparse
import csv
import os

import haldane_consistency as hc

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_INPUT = os.path.join(HERE, "curated_reactions.csv")


def score_rows(path: str) -> list[dict]:
    rows: list[dict] = []
    with open(path, newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            keq_kin = float(raw["keq_kin"])
            keq_thermo = float(raw["keq_thermo"])
            x = hc.consistency_ratio(keq_kin, keq_thermo)
            row = dict(raw)
            row["x"] = x
            row["C_haldane"] = hc.reciprocal_cost(x)
            row["fold_error"] = hc.fold_error(x)
            row["consistency_class"] = hc.classify(x)
            rows.append(row)
    return rows


def print_summary(rows: list[dict]) -> None:
    header = f"{'reaction_id':24s} {'Kkin':>7s} {'Kthermo':>8s} {'x':>6s} {'fold':>5s} {'C':>7s}  fold band"
    print(header)
    print("-" * len(header))
    for r in rows:
        print(
            f"{r['reaction_id']:24s} {float(r['keq_kin']):7.3f} "
            f"{float(r['keq_thermo']):8.3f} {r['x']:6.3f} {r['fold_error']:5.2f} "
            f"{r['C_haldane']:7.4f}  {r['consistency_class']}"
        )


def write_scored(path: str, rows: list[dict]) -> None:
    extra = ["x", "C_haldane", "fold_error", "consistency_class"]
    base = [c for c in rows[0] if c not in extra] if rows else []
    columns = base + extra
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for r in rows:
            writer.writerow({c: r.get(c, "") for c in columns})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", default=DEFAULT_INPUT,
                        help="curated reactions CSV (default: curated_reactions.csv)")
    parser.add_argument("-o", "--output", default=None,
                        help="optional path for the scored CSV")
    args = parser.parse_args()

    rows = score_rows(args.input)
    print_summary(rows)
    if args.output:
        write_scored(args.output, rows)
        print(f"\nScored CSV written to {args.output}")


if __name__ == "__main__":
    main()
