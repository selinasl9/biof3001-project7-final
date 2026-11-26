#!/usr/bin/env python3
"""
build_contig_annotation.py

Build one row per contig in clusters 0-27, the table every figure and result joins on.

Inputs:
    --info      per-sequence cluster table from get_cluster_representatives.py
                (contig, cluster, length)
    --tree      the IQ-TREE .treefile the figures are drawn from
    --refs      the reference alignment, read for its headers only
    --hits      DIAMOND best-hits-with-metadata tables, in priority order
    --metadata  run and sample metadata tables (both halves of the dataset)

Outputs:
    contig_annotation.tsv   on stdout, one row per contig, 19 columns

tree_distance is measured in whichever tree --tree names, so it is comparable only
within one table.
"""

import argparse
import csv
import re
import sys

sys.path.insert(0, __file__.rsplit("/", 2)[0])
from common import newick
from common.mosquito import NOT_A_SPECIES, genus_and_group

# Set the output columns, the per-table run column names + the text rules
COLUMNS = ["contig", "cluster", "length_bp", "run", "bioproject", "country", "province",
           "lat_lon", "host_species", "host_genus", "host_group", "collection_date",
           "diamond_species", "pident", "qcovhsp", "nearest_reference",
           "nearest_reference_name", "tree_distance", "singleton_cluster"]

# A run column is called something different in each of the three metadata tables
RUN_COLUMNS = ["acc", "sample", "run", "Run"]

# Everything after one of these words is strain or assembly detail, not the organism
DESCRIPTION_TAIL = re.compile(
    r"\s+(?:strain|isolate|from|polyprotein|genomic|genome|mRNA|complete|partial)\b")


def read_tsv(path):
    """Yield the rows of a tab-separated file as dicts."""
    with open(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def read_cluster_info(path):
    """Return (contig, cluster, length) for each sequence, skipping the header."""
    rows = []
    for row in read_tsv(path):
        values = list(row.values())
        rows.append((values[0], values[1], values[2]))
    return rows


def read_best_hits(paths, contigs):
    """Return {contig: DIAMOND best hit}, taking the first table that has one."""
    # Match without _rev, which the reverse-complemented contigs carry only in the tree
    without_rev = {contig.replace("_rev", ""): contig for contig in contigs}

    hits = {}
    for path in paths:
        for row in read_tsv(path):
            query = row.get("qseqid", "")
            contig = without_rev.get(query)
            if contig and contig not in hits:
                hits[contig] = row
    return hits


def read_run_metadata(paths):
    """Return {run: host species, province, coordinates}, merged over the metadata tables."""
    runs = {}
    for path in paths:
        for row in read_tsv(path):
            run = next((row[column].strip() for column in RUN_COLUMNS
                        if row.get(column)), "")
            if not run:
                continue

            species = (row.get("host_biosample") or row.get("organism")
                       or row.get("host_organism") or "").strip()
            # Prefer the organism column where host_biosample only says "metagenome"
            if species in NOT_A_SPECIES and row.get("organism") not in NOT_A_SPECIES:
                species = row["organism"].strip()

            entry = runs.setdefault(run, {})
            if species and not entry.get("species"):
                entry["species"] = species
            for column in ("province", "lat_lon"):
                if row.get(column):
                    entry.setdefault(column, row[column].strip())
    return runs


def host_of(run, runs):
    """Return (species, genus, group) for a run, blank or "unknown" when nothing is recorded."""
    species = runs.get(run, {}).get("species", "")
    genus, group = genus_and_group(species)
    return ("" if genus == "unknown" else species), genus, group


def read_reference_names(path):
    """Return {accession: organism name} from the reference alignment headers."""
    names = {}
    for line in open(path):
        if not line.startswith(">"):
            continue

        accession, *rest = line[1:].split()
        name = " ".join(rest).split(",")[0]
        name = DESCRIPTION_TAIL.split(name)[0].strip()
        # Drop a repeated "flavivirus", as "Aedes flavivirus flavivirus" reads badly
        names[accession] = re.sub(r"\s+flavivirus$", "", name) if name.endswith(
            "virus flavivirus") else name
    return names


def nearest_reference(contig, tips, references):
    """Return (accession, distance) for the reference genome closest to this contig."""
    tip = tips.get(contig)
    if tip is None or not references:
        return "", ""

    closest, distance = newick.nearest(tip, references)
    return closest["name"].replace("_R_", ""), f"{distance:.4f}"


def runs_per_cluster(info):
    """Return {cluster: set of runs its contigs came from}, for the singleton rule."""
    clusters = {}
    for contig, cluster, _ in info:
        clusters.setdefault(cluster, set()).add(contig.split("_")[0])
    return clusters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--info", required=True,
                        help="per-sequence cluster table (contig, cluster, length)")
    parser.add_argument("--tree", required=True, help="IQ-TREE .treefile")
    parser.add_argument("--refs", required=True, help="reference alignment, for the names")
    parser.add_argument("--hits", nargs="+", required=True,
                        help="best-hits-with-metadata tables, in priority order")
    parser.add_argument("--metadata", nargs="+", required=True,
                        help="run and sample metadata tables")
    args = parser.parse_args()

    # Step 1: read the contigs to annotate, in the order the cluster table lists them
    info = read_cluster_info(args.info)
    contigs = {contig for contig, _, _ in info}

    # Step 2: read the DIAMOND best hit for each contig, for the virus columns
    hits = read_best_hits(args.hits, contigs)

    # Step 3: read the host species, province + coordinates for each run
    runs = read_run_metadata(args.metadata)

    # Step 4: read the organism names, so the nearest reference reads as a virus name
    reference_names = read_reference_names(args.refs)

    # Step 5: read the tree + the reference tips to measure each contig against
    tips = newick.tip_dict(newick.parse(args.tree))
    # Take the tips named with an accession, which are the references
    references = [tip for name, tip in tips.items() if name.startswith("NC_")]

    # Step 6: work out which runs each cluster was found in, for the singleton column
    clusters = runs_per_cluster(info)

    # Step 7: join it all onto one row per contig + write it out
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(COLUMNS)
    for contig, cluster, length in info:
        hit = hits.get(contig, {})
        metadata = runs.get(contig.split("_")[0], {})
        species, genus, group = host_of(contig.split("_")[0], runs)
        accession, distance = nearest_reference(contig, tips, references)

        writer.writerow([
            contig, cluster, length, contig.split("_")[0],
            hit.get("bioproject", ""), hit.get("country", ""),
            metadata.get("province", ""), metadata.get("lat_lon", ""),
            species, genus, group,
            hit.get("collection_date_sam", "").strip("[]"),
            hit.get("species", ""), hit.get("pident", ""), hit.get("qcovhsp", ""),
            accession, reference_names.get(accession, ""), distance,
            "yes" if len(clusters[cluster]) == 1 else "no",
        ])

    # Step 8: return summary in the terminal, so a dropped join shows up straight away
    print(f"{len(info)} contigs, {len(clusters)} clusters, "
          f"{sum(1 for c, _, _ in info if c in hits)} with a DIAMOND hit", file=sys.stderr)


if __name__ == "__main__":
    main()
