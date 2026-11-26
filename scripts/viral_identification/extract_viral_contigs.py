#!/usr/bin/env python3
"""
extract_viral_contigs.py

Write the contigs the host-cleaned hits point at into one FASTA, ready for clustering.

A hit with qstart > qend was found on the minus strand, so that contig is reverse
complemented and its name gets _rev.

Inputs:
    {mode}_host_cleaned_hits.tsv                     both halves of the dataset
    data/raw/logan_contigs/, results/viral_identification/PRJNA1003643/logan_contigs/  whatever is on disk
    Logan on S3 for the runs that are not

Outputs:
    results/viral_identification/viral_contigs/{mode}_viral_contigs.fa        the contigs
    results/viral_identification/viral_contigs/{mode}_viral_contigs_info.tsv  one row per contig written
    results/viral_identification/viral_contigs/{mode}_missing_contigs.tsv     asked for but not found
    results/viral_identification/viral_contigs/{mode}_counts.tsv              the run's own summary
"""

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

# Set the input tables, where contigs are looked for + the output folder
# Host-cleaned hit tables, all bioprojects combined
HITS = ["results/viral_identification/host_filter/{mode}_host_cleaned_hits.tsv",
        "results/viral_identification/PRJNA1003643/host_filter/CN_PRJNA1003643_{mode}_host_cleaned_hits.tsv"]

# Where contig files may already be on disk, before downloading from Logan
CONTIG_DIRS = ["data/raw/logan_contigs", "results/viral_identification/PRJNA1003643/logan_contigs"]
LOGAN_URL = "https://s3.amazonaws.com/logan-pub/c/{run}/{run}.contigs.fa.zst"

OUTDIR = Path("results/viral_identification/viral_contigs")

# Complement of each base, for reverse complementing
COMP = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")


def load_hits(paths):
    """Read both hit tables as one, one row per contig, with its FASTA name + strand."""
    # Step 1: read both hit tables + combine into one table
    frames = []
    for p in map(Path, paths):
        if not p.exists():
            sys.exit(f"ERROR: missing input file: {p}")
        df = pd.read_csv(p, sep="\t", low_memory=False, dtype={"staxids": str})
        print(f"  {p}: {len(df):,} hits, {df['sample'].nunique():,} runs")
        frames.append(df)
    hits = pd.concat(frames, ignore_index=True)

    # Step 2: add the columns a table might not carry, so the rest can assume them
    for col in ("bioproject", "country", "species", "sscinames"):
        if col not in hits:
            hits[col] = ""

    # Take the species name, then the subject name, then "unclassified virus"
    # Blanks have to become NA first or fillna puts an empty string in the column
    hits["species"] = (hits["species"].replace("", pd.NA)
                       .fillna(hits["sscinames"].replace("", pd.NA))
                       .fillna("unclassified virus"))

    # Step 3: tidy the ids, as they come from different pipelines
    hits["sample"] = hits["sample"].astype(str).str.strip()
    hits["qseqid"] = hits["qseqid"].astype(str).str.strip()

    # Mark the minus-strand hits, where qstart > qend
    hits["reverse"] = pd.to_numeric(hits["qstart"]) > pd.to_numeric(hits["qend"])

    # Step 4: keep one row per contig, warning if a contig was hit on both strands
    strands = hits.groupby(["sample", "qseqid"])["reverse"].nunique()
    if (strands > 1).any():
        print(f"  WARNING: {int((strands > 1).sum())} contigs hit on both strands, "
              f"keeping the first row")
    hits = hits.drop_duplicates(["sample", "qseqid"]).reset_index(drop=True)

    # Step 5: build the FASTA name, _rev marking a reverse complemented contig
    hits["header"] = hits["qseqid"] + hits["reverse"].map({True: "_rev", False: ""})

    # Step 6: stop if two contigs share a name
    dup = hits["header"].duplicated().sum()
    if dup:
        sys.exit(f"ERROR: {dup} duplicated FASTA headers, would be ambiguous for CD-HIT")

    return hits[["header", "sample", "qseqid", "bioproject", "country",
                 "species", "reverse"]]


def local_file(run):
    """Return the run's contig file if it is on disk, so only the missing ones download."""
    for d in CONTIG_DIRS:
        for name in (f"{run}.contigs.fa", f"{run}.contigs.fa.zst"):
            p = Path(d) / name
            if p.exists() and p.stat().st_size:
                return p
    return None


def download(run):
    """Fetch the run's contig file from Logan, writing to .part until it is complete."""
    dest = Path(CONTIG_DIRS[0]) / f"{run}.contigs.fa.zst"
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    rc = subprocess.run(["wget", "-q", "-c", LOGAN_URL.format(run=run), "-O", str(tmp)],
                        stdin=subprocess.DEVNULL).returncode
    if rc != 0 or not tmp.exists() or tmp.stat().st_size == 0:
        tmp.unlink(missing_ok=True)
        return None
    tmp.rename(dest)
    return dest


def read_wanted(path, wanted):
    """Read the wanted contigs out of one Logan file, in one pass."""
    # Step 1: decompress as it is read, as these files are not worth unpacking to disk
    proc = None
    if path.name.endswith(".zst"):
        proc = subprocess.Popen(["zstd", "-q", "-dc", str(path)],
                                stdout=subprocess.PIPE, text=True)
        handle = proc.stdout
    else:
        handle = open(path)

    # Step 2: walk the FASTA, collecting only the contigs asked for
    found = {}
    current = None      # Name of the contig being read, None while skipping
    chunks = []          # Its sequence lines so far

    try:
        for line in handle:
            if line[0] == ">":
                if current is not None:
                    found[current] = "".join(chunks)
                    if len(found) == len(wanted):
                        current = None
                        break      # Everything found, skip the rest of the file
                cid = line[1:].split(None, 1)[0]
                current = cid if cid in wanted else None
                chunks = []
            elif current is not None:
                chunks.append(line.strip())

        if current is not None:
            found[current] = "".join(chunks)
    finally:
        handle.close()
        if proc is not None:
            proc.kill()
            proc.wait()
    return found


def process(run, wanted, keep):
    """Find or download a run's contig file, read what is wanted, drop the download."""
    path, fetched = local_file(run), False
    if path is None:
        path, fetched = download(run), True
    if path is None:
        return run, {}, "no Logan contig file"
    try:
        return run, read_wanted(path, wanted), ""
    finally:
        if fetched and not keep:
            path.unlink(missing_ok=True)


def wrap(seq, width=60):
    """Wrap a sequence in 60-character lines, as every other FASTA in the project is."""
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser(
        description="Host-cleaned contigs out of the Logan files into one FASTA, "
                    "reverse complemented where the DIAMOND hit has qstart > qend "
                    "(those headers get _rev).")
    ap.add_argument("--mode", default="strict", help="strict, loose or evalue")
    ap.add_argument("--hits", nargs="+", help="hit tables instead of the defaults")
    ap.add_argument("--outdir", default=str(OUTDIR))
    ap.add_argument("--jobs", type=int, default=4, help="runs at a time")
    ap.add_argument("--keep-downloads", action="store_true",
                    help="keep the downloaded .zst files")
    args = ap.parse_args()

    # Step 1: read the hit tables for this threshold
    paths = args.hits or [t.format(mode=args.mode) for t in HITS]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    hits = load_hits(paths)

    # Step 2: group the contig names by run, so each file is downloaded + read once
    wanted = {run: set(g) for run, g in hits.groupby("sample")["qseqid"]}
    print(f"\n{len(hits):,} contigs to extract from {len(wanted):,} runs\n")

    fasta = outdir / f"{args.mode}_viral_contigs.fa"
    lengths = {}        # FASTA name to length, also the record of what was written
    missing = []        # Contigs asked for but not found
    no_file = []        # Runs with no contig file at all
    done = 0

    # Step 3: read several runs at a time, writing each one's contigs as it comes back
    with open(fasta, "w") as out, ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(process, run, ids, args.keep_downloads)
                   for run, ids in wanted.items()]

        for fut in as_completed(futures):
            run, seqs, problem = fut.result()
            if problem:
                no_file.append(run)

            # Write this run's contigs, reverse complementing the backwards ones
            rows = hits[hits["sample"] == run]
            for header, qseqid, rev in zip(rows["header"], rows["qseqid"],
                                           rows["reverse"]):
                seq = seqs.get(qseqid)
                if not seq:
                    missing.append((run, qseqid, problem or "not in contig file"))
                    continue
                if rev:
                    seq = seq.translate(COMP)[::-1]
                out.write(f">{header}\n{wrap(seq)}\n")
                lengths[header] = len(seq)

            done += 1
            if done % 100 == 0 or done == len(wanted):
                print(f"  {done:,}/{len(wanted):,} runs, {len(lengths):,} contigs written")

    # Step 4: write one row per contig, tracing a sequence back to its run + species
    info = hits[hits["header"].isin(lengths)].copy()
    n_reverse = int(info["reverse"].sum())
    info["strand"] = info["reverse"].map({True: "minus (reverse complemented)",
                                          False: "plus"})
    info["length"] = info["header"].map(lengths)
    info = info[["header", "sample", "bioproject", "country", "species",
                 "strand", "length"]].sort_values(["sample", "header"])
    info.to_csv(outdir / f"{args.mode}_viral_contigs_info.tsv", sep="\t", index=False)

    # Step 5: list the contigs that were not found (if any)
    if missing:
        pd.DataFrame(missing, columns=["sample", "qseqid", "reason"]).to_csv(
            outdir / f"{args.mode}_missing_contigs.tsv", sep="\t", index=False)

    # Step 6: write the counts for the methods section
    counts = pd.DataFrame([{
        "mode": args.mode,
        "contigs_requested": len(hits),
        "contigs_written": len(info),
        "contigs_missing": len(missing),
        "reverse_complemented": n_reverse,
        "runs_with_contigs": info["sample"].nunique(),
        "bioprojects": info["bioproject"].nunique(),
        "virus_species": info["species"].nunique(),
        "total_bp": int(info["length"].sum()),
        "median_contig_bp": int(info["length"].median()) if len(info) else 0,
    }])
    counts.to_csv(outdir / f"{args.mode}_counts.tsv", sep="\t", index=False)

    print(f"\n{args.mode.upper()}")
    for col, val in counts.iloc[0].items():
        if col != "mode":
            print(f"  {col.replace('_', ' '):24s} {val:,}")
    if no_file:
        print(f"\n  {len(no_file)} runs had no Logan contigs: "
              f"{', '.join(no_file[:5])}")
    print(f"\nSaved {fasta}")


if __name__ == "__main__":
    main()
