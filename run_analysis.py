#!/usr/bin/env python3
"""Command-line driver for the Haldane-consistency analysis.

Reads a curated dataset (CSV), computes the reciprocal Haldane-consistency
score for every record, writes an augmented CSV, prints a short summary, and
optionally regenerates the figures.

Examples
--------
Score the bundled synthetic example and write results next to it::

    python run_analysis.py --input example_dataset.csv --output scored_example.csv

Also (re)generate Figures 1, 3 and 4 into ``figures/``::

    python run_analysis.py -i example_dataset.csv -o scored_example.csv --figures
"""

from __future__ import annotations

import argparse
import statistics
import sys
from collections import Counter
from typing import List

from haldane_consistency import (
    CLASS_LABELS,
    EnzymeRecord,
    ScoreResult,
    read_records,
    score_record,
    write_results,
)


def _score_all(records: List[EnzymeRecord]) -> List[ScoreResult]:
    results: List[ScoreResult] = []
    for rec in records:
        try:
            results.append(score_record(rec))
        except ValueError as exc:  # skip but report unscorable rows
            print(f"  [skip] {rec.record_id or '<no id>'}: {exc}", file=sys.stderr)
    return results


def _print_summary(results: List[ScoreResult]) -> None:
    if not results:
        print("No records scored.")
        return

    costs = [r.C_haldane for r in results]
    counts = Counter(r.consistency_class for r in results)

    print(f"\nScored {len(results)} record(s).")
    print(f"  C_Haldane: min={min(costs):.4g}  median={statistics.median(costs):.4g}  "
          f"max={max(costs):.4g}  mean={statistics.fmean(costs):.4g}")
    print("  Consistency classes:")
    for label in CLASS_LABELS:
        print(f"    {label:<24s} {counts.get(label, 0)}")

    worst = sorted(results, key=lambda r: r.C_haldane, reverse=True)[:5]
    print("  Largest inconsistencies:")
    for r in worst:
        print(f"    {r.record_id:<10s} x={r.ratio_x:8.3g}  "
              f"C_Haldane={r.C_haldane:8.4g}  [{r.consistency_class}]")


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compute reciprocal Haldane-consistency scores for enzyme records."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="input CSV with kinetic and thermodynamic constants")
    parser.add_argument("-o", "--output", default="scored_dataset.csv",
                        help="output CSV (default: scored_dataset.csv)")
    parser.add_argument("--figures", action="store_true",
                        help="also (re)generate Figures 1, 3 and 4")
    parser.add_argument("--figures-dir", default="figures",
                        help="directory for generated figures (default: figures)")
    args = parser.parse_args(argv)

    records = read_records(args.input)
    if not records:
        print(f"No records found in {args.input!r}.", file=sys.stderr)
        return 1

    results = _score_all(records)
    write_results(args.output, records, results)
    print(f"Wrote scored dataset to {args.output!r}.")
    _print_summary(results)

    if args.figures:
        import os

        import figures  # imported lazily so matplotlib is optional

        os.makedirs(args.figures_dir, exist_ok=True)
        f1 = figures.figure_cost_shape(os.path.join(args.figures_dir, "figure1_cost_shape.png"))
        f3 = figures.figure_consistency_scatter(
            results, os.path.join(args.figures_dir, "figure3_consistency_scatter.png"))
        f4 = figures.figure_score_distribution(
            results, os.path.join(args.figures_dir, "figure4_score_distribution.png"))
        print(f"Wrote figures: {f1}, {f3}, {f4}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
