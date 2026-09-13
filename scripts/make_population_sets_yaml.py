#!/usr/bin/env python3
"""Generate assets/population_sets.yaml for Phase 1.

Phase 1 scope (binding project ruling): only the 'all' population is
implemented. This script reads reconciliation_table.csv and writes a
single-key YAML file listing every strain with scope == INCLUDE.

Usage:
    python3 scripts/make_population_sets_yaml.py \
        --reconciliation /path/to/reconciliation_table.csv \
        --out assets/population_sets.yaml
"""
import argparse
import csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reconciliation", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    strains = []
    with open(args.reconciliation, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row["scope"] == "INCLUDE":
                strains.append(row["strain"])

    strains.sort()

    with open(args.out, "w") as out:
        out.write("all:\n")
        for s in strains:
            out.write(f"  - {s}\n")

    print(f"Wrote {len(strains)} strains to {args.out}")


if __name__ == "__main__":
    main()
