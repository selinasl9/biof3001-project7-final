#!/usr/bin/env python3
"""
make_result_tables.py

Write every number quoted in the results as a table, so each one can be checked.

Inputs:
    VIRUS_TREE    the .treefile the figures are drawn from
    HOST_TREE     the host COI .treefile, tip names like Culex_pipiens|KX...
    ANNOTATION    contig_annotation.tsv

Outputs, all into OUTDIR:
    association_tests.tsv     Cramer's V for each pair of variables, with a p value
    mantel_geography.tsv      genetic vs geographic distance, overall and per lineage
    lineage_sampling.tsv      how widely each lineage was sampled
    long_distance_pairs.tsv   near-identical contigs more than 2,000 km apart
    cophylogeny.tsv           host COI distance vs virus distance
    host_by_lineage.tsv       contigs per host species per lineage

Both trees are read for branch lengths only, so support values do not matter here.
"""

import argparse
import csv
import math
import random
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, __file__.rsplit("/", 2)[0])
from common import newick
from common.geography import parse_lat_lon
from common.palette import relabel_lineage

# Set the test settings + the thresholds the results section quotes
# The Cramer's V p values are quoted to four decimal places, which needs more shuffles,
# and a subset under MIN_SUBSET has too few pairs for a Mantel test to say anything
PERMUTATIONS = 999
PERMUTATIONS_CRAMER = 9999
RANDOM_SEED = 1

MIN_SUBSET = 6
MIN_LOCATIONS = 3
MIN_COUNTRIES = 2
SAME_VIRUS_DISTANCE = 0.05
LONG_DISTANCE_KM = 2000

EARTH_RADIUS_KM = 6371.0

ASSOCIATION_PAIRS = [("host_species", "cluster"), ("host_genus", "cluster"),
                     ("host_species", "lineage"), ("host_genus", "lineage"),
                     ("province", "cluster"), ("host_species", "province"),
                     ("country", "host_genus"), ("country", "cluster")]


def pearson(x, y):
    """Return the Pearson correlation of two equal-length lists, 0 if either is flat."""
    n = len(x)
    mean_x, mean_y = sum(x) / n, sum(y) / n
    sd_x = math.sqrt(sum((a - mean_x) ** 2 for a in x))
    sd_y = math.sqrt(sum((b - mean_y) ** 2 for b in y))
    if not sd_x or not sd_y:
        return 0.0
    return sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y)) / (sd_x * sd_y)


def haversine_km(a, b):
    """Return the great-circle distance in km between two (latitude, longitude) pairs."""
    lat1, lon1 = map(math.radians, a)
    lat2, lon2 = map(math.radians, b)
    h = (math.sin((lat2 - lat1) / 2) ** 2
         + math.cos(lat1) * math.cos(lat2) * math.sin((lon2 - lon1) / 2) ** 2)
    return 2 * EARTH_RADIUS_KM * math.asin(math.sqrt(h))


def cramers_v(pairs):
    """Return Cramer's V for a contingency table given as (a, b) pairs, with chi-square."""
    n = len(pairs)
    observed = Counter(pairs)
    rows_a, rows_b = Counter(a for a, _ in pairs), Counter(b for _, b in pairs)

    chi2 = 0.0
    for a, count_a in rows_a.items():
        for b, count_b in rows_b.items():
            expected = count_a * count_b / n
            chi2 += (observed.get((a, b), 0) - expected) ** 2 / expected

    smaller = min(len(rows_a), len(rows_b))
    v = math.sqrt(chi2 / (n * (smaller - 1))) if smaller > 1 else 0.0
    return {"n": n, "levels_a": len(rows_a), "levels_b": len(rows_b),
            "chi2": chi2, "df": (len(rows_a) - 1) * (len(rows_b) - 1), "V": v}


def permutation_p(statistic, observed, shuffle, permutations):
    """Return the share of shuffles reaching the observed statistic, +1 to keep p > 0."""
    at_least = 0
    for _ in range(permutations):
        if statistic(shuffle()) >= observed:
            at_least += 1
    return (at_least + 1) / (permutations + 1)


def read_tree(path, clean):
    """Return {tip name: node} for a tree, names passed through `clean` first."""
    root = newick.parse(path)
    return {clean(tip["name"]): tip for tip in newick.tips(root)}


def virus_tip_name(name):
    """Strip MAFFT's _R_ prefix, which is all a virus tip name carries."""
    return name.replace("_R_", "")


def host_tip_name(name):
    """Turn a host tip "Culex_pipiens|KX..." into the annotation's "Culex pipiens"."""
    return name.split("|")[0].replace("_", " ")


def read_contigs(path):
    """Return the annotation rows that are contigs, each with its lineage worked out."""
    with open(path) as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row.get("contig")]
    for row in rows:
        row["lineage"] = relabel_lineage(row["nearest_reference_name"], row["cluster"])
    return rows


def write_tsv(path, header, records):
    """Write a header row and its records to a tab-separated file."""
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(records)


def write_association_tests(rows, outdir):
    """Write Cramer's V for each pair of variables, with a permutation p value beside it."""
    records = []
    for a, b in ASSOCIATION_PAIRS:
        pairs = [(row[a], row[b]) for row in rows if row.get(a) and row.get(b)]
        if not pairs:
            continue

        result = cramers_v(pairs)
        first = [a for a, _ in pairs]
        second = [b for _, b in pairs]

        # Shuffle one variable, which breaks the link but keeps both sets of row totals
        def shuffle():
            random.shuffle(second)
            return list(zip(first, second))

        p = permutation_p(lambda pairs: cramers_v(pairs)["V"], result["V"],
                          shuffle, PERMUTATIONS_CRAMER)
        records.append([a, b, result["n"], result["levels_a"], result["levels_b"],
                        round(result["chi2"], 1), result["df"], round(result["V"], 3),
                        p, PERMUTATIONS_CRAMER])

    write_tsv(f"{outdir}/association_tests.tsv",
              ["variable_a", "variable_b", "n", "levels_a", "levels_b", "chi2", "df",
               "cramers_V", "p_perm", "permutations"], records)


def mantel(subset, distance):
    """Run a Mantel test of geographic against genetic distance for one set of contigs."""
    geographic, genetic = [], []
    for i, a in enumerate(subset):
        for b in subset[i + 1:]:
            geographic.append(haversine_km(parse_lat_lon(a["lat_lon"]),
                                           parse_lat_lon(b["lat_lon"])))
            genetic.append(distance(a["contig"], b["contig"]))

    observed = pearson(geographic, genetic)
    order = list(range(len(subset)))

    # Shuffle which contig sits at which place, keeping both distance matrices intact
    def shuffle():
        random.shuffle(order)
        return [distance(subset[order[i]]["contig"], subset[order[j]]["contig"])
                for i in range(len(subset)) for j in range(i + 1, len(subset))]

    p = permutation_p(lambda permuted: pearson(geographic, permuted), observed,
                      shuffle, PERMUTATIONS)
    return [len(subset), len(geographic), round(observed, 3), round(p, 3), PERMUTATIONS]


def write_mantel_geography(rows, distance, in_tree, outdir):
    """Write the Mantel test overall, within China + within each lineage big enough."""
    located = [row for row in rows
               if parse_lat_lon(row["lat_lon"]) and row["contig"] in in_tree]

    records = [["all contigs with coordinates", *mantel(located, distance)]]
    china = [row for row in located if row["country"] == "China"]
    records.append(["China only", *mantel(china, distance)])

    by_lineage = defaultdict(list)
    for row in located:
        by_lineage[row["lineage"]].append(row)
    for lineage, subset in sorted(by_lineage.items(), key=lambda item: -len(item[1])):
        if len(subset) >= MIN_SUBSET:
            records.append([f"within {lineage}", *mantel(subset, distance)])

    write_tsv(f"{outdir}/mantel_geography.tsv",
              ["subset", "n_contigs", "n_pairs", "mantel_r", "p_value", "permutations"],
              records)
    return located


def write_lineage_sampling(rows, located, distance, outdir):
    """Write how many places and countries each lineage was found in, and how far apart."""
    by_lineage = defaultdict(list)
    for row in rows:
        by_lineage[row["lineage"]].append(row)
    located_by_lineage = defaultdict(list)
    for row in located:
        located_by_lineage[row["lineage"]].append(row)

    records = []
    for lineage, subset in sorted(by_lineage.items(), key=lambda item: -len(item[1])):
        countries = sorted({row["country"] for row in subset})
        places = {row["province"] or row["country"] for row in subset}

        with_coordinates = located_by_lineage[lineage]
        spans, distances = [], []
        for i, a in enumerate(with_coordinates):
            for b in with_coordinates[i + 1:]:
                spans.append(haversine_km(parse_lat_lon(a["lat_lon"]),
                                          parse_lat_lon(b["lat_lon"])))
                distances.append(distance(a["contig"], b["contig"]))

        identities = [float(row["pident"]) for row in subset if row.get("pident")]
        records.append([
            lineage, len(subset), len(countries), len(places), "; ".join(countries),
            round(max(spans), 0) if spans else 0,
            round(median(distances), 4) if distances else "",
            min(identities) if identities else "",
            max(identities) if identities else "",
            "yes" if len(places) >= MIN_LOCATIONS else "no",
            "yes" if len(countries) >= MIN_COUNTRIES else "no",
        ])

    write_tsv(f"{outdir}/lineage_sampling.tsv",
              ["lineage", "n_contigs", "n_countries", "n_locations", "countries",
               "max_span_km", "median_genetic_distance", "pident_min", "pident_max",
               f"passes_{MIN_LOCATIONS}_locations", f"passes_{MIN_COUNTRIES}_countries"],
              records)


def write_long_distance_pairs(located, distance, outdir):
    """Write the pairs of near-identical contigs found a long way apart, furthest first."""
    found = []
    for i, a in enumerate(located):
        for b in located[i + 1:]:
            genetic = distance(a["contig"], b["contig"])
            kilometres = haversine_km(parse_lat_lon(a["lat_lon"]),
                                      parse_lat_lon(b["lat_lon"]))
            if genetic < SAME_VIRUS_DISTANCE and kilometres > LONG_DISTANCE_KM:
                found.append((kilometres, genetic, a, b))

    records = [[round(kilometres), round(genetic, 4), a["lineage"],
                a["contig"], a["province"] or a["country"], a["host_species"],
                b["contig"], b["province"] or b["country"], b["host_species"]]
               for kilometres, genetic, a, b in sorted(found, key=lambda item: -item[0])]

    write_tsv(f"{outdir}/long_distance_pairs.tsv",
              ["km_apart", "genetic_distance", "lineage", "contig_a", "place_a", "host_a",
               "contig_b", "place_b", "host_b"], records)


def median(values):
    """Return the middle value of a list, taking the upper of two middles."""
    return sorted(values)[len(values) // 2]


def write_cophylogeny(rows, virus_distance, host_distance, in_virus_tree,
                      host_tips, outdir):
    """Write host COI distance against virus distance, over every pair of contigs."""
    usable = [row for row in rows if row["contig"] in in_virus_tree
              and (row["host_species"] or "") in host_tips]

    host, virus, same_species = [], [], []
    for i, a in enumerate(usable):
        for b in usable[i + 1:]:
            host.append(0.0 if a["host_species"] == b["host_species"]
                        else host_distance(a["host_species"], b["host_species"]))
            virus.append(virus_distance(a["contig"], b["contig"]))
            same_species.append(a["host_species"] == b["host_species"])

    observed = pearson(host, virus)
    order = list(range(len(usable)))

    # Shuffle which host species each contig came from, keeping the virus distances
    def shuffle():
        random.shuffle(order)
        permuted = []
        for i in range(len(usable)):
            for j in range(i + 1, len(usable)):
                a, b = usable[order[i]]["host_species"], usable[order[j]]["host_species"]
                permuted.append(0.0 if a == b else host_distance(a, b))
        return permuted

    p = permutation_p(lambda permuted: pearson(permuted, virus), observed,
                      shuffle, PERMUTATIONS)

    same = [v for v, flag in zip(virus, same_species) if flag]
    different = [v for v, flag in zip(virus, same_species) if not flag]
    write_tsv(f"{outdir}/cophylogeny.tsv", ["measure", "value", "note"], [
        ["contigs used", len(usable), "have both a virus tip and a host in the COI tree"],
        ["pairs", len(host), ""],
        ["mantel_r", round(observed, 3),
         "host COI distance vs virus patristic distance"],
        ["p_value", round(p, 3), f"{PERMUTATIONS} permutations"],
        ["median virus distance, same host species", round(median(same), 4),
         f"n={len(same)} pairs"],
        ["median virus distance, different host species", round(median(different), 4),
         f"n={len(different)} pairs"],
        ["caveat", "", "r is inflated: same-species pairs are 0 in both matrices"],
    ])


def write_host_by_lineage(rows, outdir):
    """Write the contigs per host species per lineage, the table behind host_vs_virus.pdf."""
    carriage = defaultdict(Counter)
    for row in rows:
        if row["host_species"]:
            carriage[row["host_species"]][row["lineage"]] += 1

    lineages = sorted({row["lineage"] for row in rows})
    records = [[species, *[carriage[species].get(lineage, 0) for lineage in lineages],
                sum(carriage[species].values())] for species in sorted(carriage)]

    write_tsv(f"{outdir}/host_by_lineage.tsv",
              ["host_species", *lineages, "total"], records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("virus_tree")
    parser.add_argument("host_tree")
    parser.add_argument("annotation")
    parser.add_argument("outdir", nargs="?", default="results/tables")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    # Fix the seed, so the permutation p values come out the same on a re-run
    random.seed(RANDOM_SEED)

    # Step 1: read both trees, as a lookup from tip name to node
    virus_tips = read_tree(args.virus_tree, virus_tip_name)
    host_tips = read_tree(args.host_tree, host_tip_name)

    def virus_distance(a, b):
        return newick.patristic_distance(virus_tips[a], virus_tips[b])

    def host_distance(a, b):
        return newick.patristic_distance(host_tips[a], host_tips[b])

    # Step 2: read the contigs, with the relabelled lineage each belongs to
    rows = read_contigs(args.annotation)
    print(f"{len(rows)} contigs, {len(virus_tips)} virus tips, {len(host_tips)} host tips")

    # Step 3: measure how strongly each pair of variables is associated
    write_association_tests(rows, args.outdir)

    # Step 4: test genetic against geographic distance, which also gives the located
    # contigs, reused by the next two tables
    located = write_mantel_geography(rows, virus_distance, virus_tips, args.outdir)

    # Step 5: measure how widely each lineage was sampled + find the long-distance pairs
    write_lineage_sampling(rows, located, virus_distance, args.outdir)
    write_long_distance_pairs(located, virus_distance, args.outdir)

    # Step 6: test the host tree against the virus tree
    write_cophylogeny(rows, virus_distance, host_distance, virus_tips,
                      host_tips, args.outdir)

    # Step 7: write the host x lineage counts behind the host figure
    write_host_by_lineage(rows, args.outdir)

    print("wrote into", args.outdir)
    for name in ["association_tests", "mantel_geography", "lineage_sampling",
                 "long_distance_pairs", "cophylogeny", "host_by_lineage"]:
        print("  ", name + ".tsv")


if __name__ == "__main__":
    main()
