#!/usr/bin/env python3
"""
fasta_tools.py

Run the six FASTA jobs this project needs, as subcommands sharing one reader.

    gaps          how gappy an alignment is, printed, nothing written
    ungap         strip gaps so sequences can go into the BLAST web form
    dedup         de-duplicate the sequences downloaded from BLAST
    drop          remove named sequences and the gap-only columns they leave
    trim-gappy    drop columns that too few references cover
    split-length  split a FASTA into long / moderate / short by sequence length

Inputs:
    one or more FASTA files, per subcommand

Outputs:
    a FASTA, per subcommand, plus a table for dedup and split-length
"""

import argparse
import re
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from common.fasta_io import read_fasta, sequence_id, write_fasta

# Set the length window, the gap characters, the header patterns + the occupancy cutoffs
# MIN_GENOME to MAX_GENOME is the complete-genome window for this family, in bases
MIN_GENOME = 5000
MAX_GENOME = 11500

GAP_CHARACTERS = "-."

# Short sequences carrying SEGMENT_WORD are segments of a segmented genome (Jingmen group)
SEGMENT_WORD = "segment"

# NCBI downloads come as gi|123|ref|NC_001.1| or lcl|Query_1, so dig the accession out
ACCESSION = re.compile(r"([A-Z]{1,2}_?[0-9]{5,8}(?:\.[0-9]+)?)")

OCCUPANCY_CUTOFFS = (0.1, 0.3, 0.5, 0.8)


def ungap_sequence(seq):
    """Return one sequence with the alignment gaps taken out."""
    for character in GAP_CHARACTERS:
        seq = seq.replace(character, "")
    return seq.upper()


def accession(header):
    """Return the accession in a header, or the first token if there isn't one."""
    first = header.split()[0]
    match = ACCESSION.search(first)
    return match.group(1) if match else first


def virus_name(header):
    """Return everything up to and including the first word ending in "virus"."""
    words = header.split()[1:]
    for i, word in enumerate(words):
        if word.lower().rstrip(",").endswith("virus"):
            return " ".join(words[:i + 1]).lower().rstrip(",")
    return " ".join(words[:3]).lower()


def cmd_gaps(args):
    """Print how gappy each alignment is, to choose a trimming cutoff by."""
    for path in args.alignment:
        # Step 1: read the alignment
        records = list(read_fasta(path))
        names = [name.split()[0] for name, _ in records]
        seqs = [seq for _, seq in records]
        width = len(seqs[0])

        # Step 2: count gaps per sequence + how many sequences cover each column
        per_sequence = [seq.count("-") / width for seq in seqs]
        occupancy = [sum(1 for seq in seqs if seq[i] != "-") / len(seqs)
                     for i in range(width)]

        # Step 3: return the counts in the terminal, a line per candidate cutoff
        print(f"\n{path}")
        print(f"  {len(seqs)} sequences x {width} columns")
        print(f"  mean gaps per sequence: {100 * sum(per_sequence) / len(per_sequence):.1f}%")
        for cutoff in OCCUPANCY_CUTOFFS:
            n = sum(1 for o in occupancy if o < cutoff)
            print(f"  columns filled by < {int(100 * cutoff):>2}% of sequences: "
                  f"{n:>6}  ({100 * n / width:.1f}%)")

        # Step 4: name the gappiest sequences, usually the ones to drop
        worst = sorted(zip(per_sequence, names), reverse=True)[:5]
        print("  gappiest sequences: "
              + ", ".join(f"{name} {100 * gap:.0f}%" for gap, name in worst))


def cmd_ungap(args):
    """Strip the gaps out so the sequences can go into the BLAST web form."""
    # Step 1: strip the gaps, holding back anything too short to be worth a BLAST
    kept, short = [], []
    for header, seq in read_fasta(args.alignment):
        header = sequence_id(header) if header.startswith("_R_") else header
        bases = ungap_sequence(seq)

        if len(bases) < args.min_length:
            short.append((bases and header.split()[0], len(bases)))
            continue
        kept.append((header, bases))

    # Step 2: write the ungapped sequences + report the length range
    write_fasta(args.out, kept)
    lengths = sorted(len(seq) for _, seq in kept)
    print(f"{len(kept)} sequences, {lengths[0]}-{lengths[-1]} bp, {sum(lengths)} bp in total")
    print("wrote", args.out)

    if short:
        print(f"left out {len(short)} sequences under {args.min_length} bp: "
              + ", ".join(f"{name} {length}" for name, length in short))

    # Step 3: write the sequences in batches if --chunk, as web BLAST takes one job per query
    if args.chunk:
        stem = args.out.rsplit(".", 1)[0]
        for start in range(0, len(kept), args.chunk):
            part = f"{stem}_part{start // args.chunk + 1}.fa"
            write_fasta(part, kept[start:start + args.chunk])
            print(f"wrote {part} ({len(kept[start:start + args.chunk])} sequences)")


def cmd_dedup(args):
    """De-duplicate BLAST downloads, with a verdict for every record."""
    # Step 1: read the accessions + sequences we already have, so they are not kept twice
    have_accessions, have_sequences = set(), set()
    for path in args.exclude:
        for header, seq in read_fasta(path):
            have_accessions.add(accession(header))
            have_sequences.add(ungap_sequence(seq))
    print(f"{len(have_accessions)} sequences we already have, from {len(args.exclude)} files")

    # Step 2: give one verdict per record, in the order the rules are checked
    kept, rows = [], []
    seen_accessions, seen_sequences = set(), set()
    counts = dict.fromkeys(["kept", "same accession", "same sequence", "already have it",
                            "too short", "too long", "query sequence"], 0)

    for path in args.fasta:
        for header, raw in read_fasta(path):
            seq = ungap_sequence(raw)
            acc = accession(header)
            if not seq:
                continue

            # The aligned download includes the query itself as Query_1 / lcl|Query_1
            if re.search(r"query", header.split()[0], re.I):
                verdict = "query sequence"
            elif acc in have_accessions or seq in have_sequences:
                verdict = "already have it"
            elif acc in seen_accessions:
                verdict = "same accession"
            elif seq in seen_sequences:
                verdict = "same sequence"
            elif len(seq) < args.min_length:
                verdict = "too short"
            elif len(seq) > args.max_length:
                verdict = "too long"
            else:
                verdict = "kept"
                seen_accessions.add(acc)
                seen_sequences.add(seq)
                kept.append((header, seq))

            counts[verdict] += 1
            rows.append((acc, len(seq), verdict, header))

    # Step 3: return the verdict counts in the terminal, then write what is left
    for verdict, n in counts.items():
        print(f"  {verdict}: {n}")

    write_fasta(args.out, kept)
    lengths = sorted(len(seq) for _, seq in kept)
    print(f"{len(kept)} sequences left, {lengths[0]}-{lengths[-1]} bp" if kept else "nothing left")
    print("wrote", args.out)

    # Step 4: write the per-record report, so a surprising count can be chased down
    if args.report:
        with open(args.report, "w") as out:
            out.write("accession\tlength_bp\tverdict\theader\n")
            for row in rows:
                out.write("\t".join(str(field) for field in row) + "\n")
        print("wrote", args.report)


def cmd_drop(args):
    """Remove named sequences, then the columns that are gaps in everything left."""
    # Step 1: read the alignment + check every name to drop is actually in it
    records = list(read_fasta(args.alignment))
    names = [sequence_id(header) for header, _ in records]

    missing = [name for name in args.drop if name not in names]
    if missing:
        raise SystemExit(f"not found in the alignment: {', '.join(missing)}")

    # Step 2: drop the named sequences
    kept = [record for record, name in zip(records, names) if name not in args.drop]
    print(f"removed {len(records) - len(kept)} sequences, {len(kept)} left")

    # Step 3: delete the gap-only columns, as MEGA does, if every sequence is the same length
    width = len(kept[0][1])
    if len({len(seq) for _, seq in kept}) == 1:
        columns = [i for i in range(width) if any(seq[i] != "-" for _, seq in kept)]
        print(f"columns: {width} -> {len(columns)} "
              f"({width - len(columns)} gap-only columns removed)")
        kept = [(header, "".join(seq[i] for i in columns)) for header, seq in kept]
    else:
        print("sequences have different lengths, so this is not an alignment, "
              "no columns removed")

    # Step 4: write it out
    write_fasta(args.out, kept)
    print("wrote", args.out)


def cmd_trim_gappy(args):
    """Keep the columns enough reference sequences cover, drop the rest."""
    # Step 1: read the alignment
    records = list(read_fasta(args.alignment))
    if not records:
        sys.exit(f"ERROR: no sequences in {args.alignment}")
    width = len(records[0][1])

    # Step 2: mark which sequences are references, as the columns are judged on those only
    reference_names = {header.split()[0] for header, _ in read_fasta(args.references)}
    is_reference = [sequence_id(header) in reference_names for header, _ in records]

    n_references = sum(is_reference)
    if n_references == 0:
        sys.exit("ERROR: no reference names matched the alignment")

    # Step 3: count reference coverage per column, as our contigs are far too short to judge on
    counts = [0] * width
    for keep, (_, seq) in zip(is_reference, records):
        if not keep:
            continue
        for i, character in enumerate(seq):
            if character not in GAP_CHARACTERS:
                counts[i] += 1

    # Step 4: keep the columns that pass the occupancy cutoff
    cutoff = args.min_occupancy * n_references
    columns = [i for i, count in enumerate(counts) if count >= cutoff]
    if not columns:
        sys.exit("ERROR: no columns passed the cutoff, try a lower --min-occupancy")

    # Step 5: trim every sequence to those columns
    out_records, empty = [], 0
    for header, seq in records:
        trimmed = "".join(seq[i] for i in columns)

        # A contig that sat entirely in the trimmed-away columns carries nothing
        if set(trimmed) <= set(GAP_CHARACTERS):
            empty += 1
            if args.drop_empty:
                continue
        out_records.append((header, trimmed))

    # Step 6: write it + return what the trim cost in the terminal
    write_fasta(args.out, out_records)

    print(f"{n_references} references, {len(records)} sequences total")
    print(f"{width} columns -> {len(columns)} kept "
          f"({len(columns) / width * 100:.1f}%) at {args.min_occupancy:.0%} occupancy")
    print(f"{empty} sequences are all gaps after trimming"
          f"{', left out' if args.drop_empty else ', kept, use --drop-empty to leave them out'}")
    print(f"saved {args.out}")


def cmd_split_length(args):
    """Split a FASTA into long / moderate / short, segments told apart from fragments."""
    # Step 1: read the sequences, ungapped, so length means bases
    records = [(header, ungap_sequence(seq)) for header, seq in read_fasta(args.fasta)]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = args.prefix or Path(args.fasta).name.split(".")[0]

    # Step 2: note the segmented viruses, as some segments are named after their gene
    segmented = {virus_name(header) for header, seq in records
                 if len(seq) < args.min_length and SEGMENT_WORD in header.lower()}

    # Step 3: sort each sequence into its length group
    groups = {"long": [], "moderate": [], "short_segmented": [], "short_nonsegmented": []}
    lengths = []

    for header, seq in records:
        if len(seq) > args.max_length:
            group = "long"
        elif len(seq) < args.min_length:
            group = ("short_segmented"
                     if SEGMENT_WORD in header.lower() or virus_name(header) in segmented
                     else "short_nonsegmented")
        else:
            group = "moderate"

        groups[group].append((header, seq))
        lengths.append((header.split()[0], len(seq), group))

    if segmented:
        print(f"segmented viruses found: {', '.join(sorted(segmented))}")

    # Step 4: write one FASTA per group + a table of every sequence's length and group
    for group, group_records in groups.items():
        out = outdir / f"{stem}_{group}.fa"
        write_fasta(out, group_records)
        print(f"{group:<19} {len(group_records):>5}  {out}")

    lengths_path = outdir / f"{stem}_lengths.tsv"
    with open(lengths_path, "w") as out:
        out.write("sequence\tlength\tgroup\n")
        for sid, length, group in lengths:
            out.write(f"{sid}\t{length}\t{group}\n")

    # Step 5: list the short non-segmented ones, the ones worth looking at by eye
    if groups["short_nonsegmented"]:
        print("short, not from a segmented virus:")
        for header, seq in groups["short_nonsegmented"]:
            print(f"  {len(seq):>6} bp  {header[:90]}")

    # Step 6: return summary in the terminal
    print(f"total               {len(lengths):>5}  "
          f"(moderate = {args.min_length}-{args.max_length} bp, the near-complete genomes)")
    print(f"lengths table       {lengths_path}")


def build_parser():
    """Return the argument parser, with one subparser per subcommand."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    subparsers = parser.add_subparsers(dest="command", required=True)

    gaps = subparsers.add_parser("gaps", help="report how gappy an alignment is")
    gaps.add_argument("alignment", nargs="+")
    gaps.set_defaults(run=cmd_gaps)

    ungap = subparsers.add_parser("ungap", help="strip alignment gaps")
    ungap.add_argument("alignment")
    ungap.add_argument("-o", "--out", required=True)
    ungap.add_argument("--chunk", type=int, default=0,
                       help="also write the sequences in batches of this many")
    ungap.add_argument("--min-length", type=int, default=0)
    ungap.set_defaults(run=cmd_ungap)

    dedup = subparsers.add_parser("dedup", help="de-duplicate BLAST downloads")
    dedup.add_argument("fasta", nargs="+")
    dedup.add_argument("-o", "--out", required=True)
    dedup.add_argument("--exclude", nargs="*", default=[],
                       help="FASTA files whose sequences we already have")
    dedup.add_argument("--min-length", type=int, default=MIN_GENOME)
    dedup.add_argument("--max-length", type=int, default=MAX_GENOME)
    dedup.add_argument("--report", help="write a verdict per record")
    dedup.set_defaults(run=cmd_dedup)

    drop = subparsers.add_parser("drop", help="remove named sequences")
    drop.add_argument("alignment")
    drop.add_argument("out")
    drop.add_argument("--drop", nargs="+", required=True,
                      help="sequence names, the first word of the header")
    drop.set_defaults(run=cmd_drop)

    trim = subparsers.add_parser("trim-gappy", help="drop columns few references cover")
    trim.add_argument("alignment", help="references + contigs, aligned")
    trim.add_argument("references", help="the reference FASTA")
    trim.add_argument("-o", "--out", required=True)
    trim.add_argument("--min-occupancy", type=float, default=0.5)
    trim.add_argument("--drop-empty", action="store_true",
                      help="leave out sequences with no bases left")
    trim.set_defaults(run=cmd_trim_gappy)

    split = subparsers.add_parser("split-length", help="split by sequence length")
    split.add_argument("fasta")
    split.add_argument("outdir")
    split.add_argument("--prefix", help="name the output files")
    split.add_argument("--min-length", type=int, default=MIN_GENOME)
    split.add_argument("--max-length", type=int, default=MAX_GENOME)
    split.set_defaults(run=cmd_split_length)

    return parser


if __name__ == "__main__":
    arguments = build_parser().parse_args()
    arguments.run(arguments)
