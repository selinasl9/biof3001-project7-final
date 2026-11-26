#!/usr/bin/env python3
"""
get_cluster_representatives.py

Pull each long-enough CD-HIT cluster's representative, or all its members, from the source FASTA.

Inputs:
    clstr           CD-HIT .clstr file
    source_fasta    the FASTA the clustering was run on

Outputs:
    out_fasta       the chosen sequences, longest first
    <stem>_info.tsv sequence, cluster, length
    <stem>_ids.txt  one sequence name per line
"""

import argparse
import re
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from common.fasta_io import read_fasta, sequence_id, write_fasta

# Set the .clstr member-line pattern + the minimum representative length
# Member line e.g. "0  10771nt, >SRR22751784_4127... *", non-representatives end in "at ..."
MEMBER_LINE = re.compile(r"^\d+\s+(\d+)nt,\s+>(.+?)\.\.\.\s+(\*|at )")

MIN_LENGTH = 5000


def read_clusters(path):
    """Return {cluster: [(name, length, is_representative), ...]} from a .clstr file."""
    clusters = {}
    cluster = None

    for line in open(path):
        if line.startswith(">Cluster"):
            cluster = line.split()[1]
            clusters[cluster] = []
            continue

        match = MEMBER_LINE.match(line)
        if match:
            clusters[cluster].append((match.group(2), int(match.group(1)),
                                      match.group(3) == "*"))

    return clusters


def select_sequences(clusters, min_length, all_members):
    """Return {name: (cluster, length)} for the sequences wanted from clusters that pass."""
    wanted = {}

    for cluster, members in clusters.items():
        representative = next((m for m in members if m[2]), None)

        # Skip clusters whose representative is too short to be a complete genome
        if representative is None or representative[1] < min_length:
            continue

        for name, length, is_representative in members:
            if all_members or is_representative:
                wanted[name] = (cluster, length)

    return wanted


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("clstr")
    parser.add_argument("source_fasta")
    parser.add_argument("out_fasta")
    parser.add_argument("--min-length", type=int, default=MIN_LENGTH)
    parser.add_argument("--all-members", action="store_true",
                        help="every contig of those clusters, not just the representative")
    args = parser.parse_args()

    # Step 1: read the clusters + pick the names we want
    clusters = read_clusters(args.clstr)
    wanted = select_sequences(clusters, args.min_length, args.all_members)

    kept = sorted({int(c) for c, _ in wanted.values()})
    what = "contigs" if args.all_members else "representatives"
    print(f"clusters with a representative >= {args.min_length} bp: {len(kept)} "
          f"(cluster {kept[0]} to {kept[-1]})")
    print(f"{what} wanted: {len(wanted)}")

    # Step 2: look the names up in the source FASTA
    found = []
    for header, seq in read_fasta(args.source_fasta):
        sid = sequence_id(header)
        if sid in wanted:
            found.append((sid, seq))

    missing = sorted(set(wanted) - {sid for sid, _ in found})
    if missing:
        print(f"WARNING: {len(missing)} not in {args.source_fasta}: "
              f"{', '.join(missing[:5])}")

    # Step 3: sort longest first, so the alignment starts with the fullest genomes
    found.sort(key=lambda record: -len(record[1]))

    # Step 4: write the FASTA
    out = Path(args.out_fasta)
    out.parent.mkdir(parents=True, exist_ok=True)
    stem = out.name.split(".")[0]
    info = out.with_name(stem + "_info.tsv")
    ids = out.with_name(stem + "_ids.txt")

    write_fasta(out, found)
    print(f"wrote {len(found)} sequences to {out}")

    # Step 5: write the ids + info tables for the tree annotation and later commands
    ids.write_text("".join(f"{sid}\n" for sid, _ in found))
    print(f"wrote {ids} (one name per line)")

    rows = ["sequence\tcluster\tlength\n"]
    rows += [f"{sid}\t{wanted[sid][0]}\t{len(seq)}\n" for sid, seq in found]
    info.write_text("".join(rows))
    print(f"wrote {info} (sequence, cluster, length)")


if __name__ == "__main__":
    main()
