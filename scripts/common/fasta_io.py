#!/usr/bin/env python3
"""
fasta_io.py

Read and write FASTA, in one place instead of a copy in each script.

Used by: fasta_tools.py, get_cluster_representatives.py, coi/trim_coi_alignment.py,
coi/fetch_host_coi_references.py
"""


def read_fasta(path):
    """Yield (header, sequence) for each record, header is everything after '>'."""
    # Read line by line rather than load the file, as the reference alignment is 39,278 columns
    name, chunks = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks)
            name, chunks = line[1:].rstrip("\n"), []
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def read_fasta_dict(path):
    """Return the whole file as {header: sequence}, for callers that need random access."""
    return dict(read_fasta(path))


def write_fasta(path, records, width=60):
    """Write (header, sequence) pairs, wrapped at width characters per line."""
    with open(path, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for start in range(0, len(seq), width):
                out.write(seq[start:start + width] + "\n")


def sequence_id(header):
    """Return the first token of a header, with MAFFT's reverse-complement marker removed."""
    sid = header.split()[0]
    # Strip _R_ at the front only, as a plain replace would cut it from the middle of a name
    return sid[3:] if sid.startswith("_R_") else sid
