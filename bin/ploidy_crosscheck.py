#!/usr/bin/env python3
"""Cross-check computationally inferred ploidy against metadata.txt's
existing ploidy column, producing ploidy_review.csv for manual sign-off.

Agreement is NOT auto-finalized - a human still produces the actual
ploidy_overrides.csv the pipeline reads for calling (see the design spec's
"Ploidy inference" section). This script only flags what needs attention.
"""
import argparse
import csv


def crosscheck(inferred_rows, metadata_ploidy_by_strain):
    results = []
    for row in inferred_rows:
        strain = row["strain"]
        inferred = row["inferred_ploidy"]
        meta = metadata_ploidy_by_strain.get(strain, "")
        meta_norm = "haploid" if meta == "haploid_from_hybrid" else meta
        if inferred == "unknown":
            status = "UNKNOWN_INFERENCE"
        elif not meta_norm:
            status = "NO_METADATA"
        elif meta_norm == inferred:
            status = "AGREE"
        else:
            status = "DISAGREE"
        results.append({
            "strain": strain,
            "inferred_ploidy": inferred,
            "metadata_ploidy": meta,
            "status": status,
        })
    return results


def load_inferred_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def load_metadata_ploidy(path):
    by_strain = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            strain = row["strain"].strip()
            ploidy = row.get("ploidy", "").strip()
            if strain not in by_strain or ploidy:
                by_strain[strain] = ploidy
    return by_strain


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inferred", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    inferred_rows = load_inferred_csv(args.inferred)
    metadata_ploidy = load_metadata_ploidy(args.metadata)
    results = crosscheck(inferred_rows, metadata_ploidy)

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["strain", "inferred_ploidy", "metadata_ploidy", "status"])
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()
