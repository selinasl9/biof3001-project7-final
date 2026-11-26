#!/usr/bin/env bash
#
# run_diamond_blastx.sh
#
# DIAMOND blastx the Logan contigs against the flavivirus protein database, run by run.
#
# Inputs:
#   PROJECT    the BioProject to screen, e.g. PRJNA1003643
#   CC         a country code for the file names, e.g. CN
#   runs/logan_matches_rna_stats.tsv   the run table, for this project's runs
#   db/arbovirus_refseq                the DIAMOND database, with taxonomy
#
# Outputs:
#   results/PROJECT/diamond/RUN.blastx.tsv   one hit table per run
#   results/PROJECT/summary/                 best hits, virus matrices, timing
#
set -u

PROJECT="${1:?give a BioProject, e.g. PRJNA1003643}"
CC="${2:?give a country code for file names, e.g. CN}"
PREFIX="${CC}_${PROJECT}"

# Set the job sizes, the input tables + the output folders
# JOBS x THREADS must stay under the CPU count, check nproc
JOBS=4
THREADS=8

RUN_TABLE="runs/logan_matches_rna_stats.tsv"
DB="db/arbovirus_refseq"
TAXDUMP="db/taxonomy"
OUTDIR="results/${PROJECT}"
CONTIG_DIR="data/logan_contigs"

mkdir -p "$OUTDIR/diamond" "$OUTDIR/timing" "$OUTDIR/summary" "$CONTIG_DIR"

# Step 1: build the run list for this BioProject
RUN_LIST="$OUTDIR/runs.txt"
awk -F'\t' -v p="$PROJECT" '$2==p {print $1}' "$RUN_TABLE" > "$RUN_LIST"
echo "$(wc -l < "$RUN_LIST") runs in $PROJECT, $JOBS at a time x $THREADS threads"
date

# Step 2: download the contigs + DIAMOND one run, called in parallel below
run_one() {
    acc="$1"
    ZST="${CONTIG_DIR}/${acc}.contigs.fa.zst"
    FA="${CONTIG_DIR}/${acc}.contigs.fa"
    OUT="${OUTDIR}/diamond/${acc}.blastx.tsv"
    TIMEFILE="${OUTDIR}/timing/${acc}.time.txt"

    [ -s "$OUT" ] && { echo "ALREADY DONE: $acc"; return; }

    if ! wget -q -c "https://s3.amazonaws.com/logan-pub/c/${acc}/${acc}.contigs.fa.zst" -O "$ZST" < /dev/null; then
        echo "MISSING: $acc (no Logan contigs)"
        rm -f "$ZST"
        return
    fi
    zstd -q -dc "$ZST" > "$FA"

    # Keep the pilot's settings, writing to .tmp so a crash is retried next time
    if /usr/bin/time -v -o "$TIMEFILE" \
        diamond blastx --query "$FA" --db "$DB" --out "$OUT.tmp" \
        --threads "$THREADS" --sensitive --evalue 1e-5 --max-target-seqs 5 \
        --header simple --outfmt 6 \
        qseqid sseqid pident length qlen slen qstart qend sstart send \
        evalue bitscore qcovhsp scovhsp staxids sscinames stitle \
        --quiet < /dev/null
    then
        mv "$OUT.tmp" "$OUT"
        echo "DONE: $acc  $(date +%H:%M)  [$(ls "$OUTDIR/diamond" | grep -c '\.blastx\.tsv$') finished]"
    else
        echo "FAILED: $acc (see $TIMEFILE)"
        rm -f "$OUT.tmp"
    fi
    rm -f "$FA" "$ZST"
}
export -f run_one
export OUTDIR CONTIG_DIR DB THREADS

xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {} < "$RUN_LIST"

echo "All runs processed"
date

# Step 3: write the summary tables
python3 - "$OUTDIR" "$RUN_TABLE" "$TAXDUMP" "$PREFIX" <<'PYEOF'
import glob, os, re, sys
import pandas as pd

outdir, run_table, taxdump, prefix = sys.argv[1:5]
summ = os.path.join(outdir, "summary")

meta = (pd.read_csv(run_table, sep="\t", dtype=str).fillna("")
          [["acc", "bioproject", "geo_loc_name_country_calc", "seqstats_contigs_nbseq"]]
          .drop_duplicates("acc")
          .rename(columns={"acc": "sample", "geo_loc_name_country_calc": "country",
                           "seqstats_contigs_nbseq": "n_contigs"}))
meta["n_contigs"] = meta["n_contigs"].str.removesuffix(".0")

# Take the best hit per contig, highest bitscore, ties broken by longer alignment
files = sorted(glob.glob(os.path.join(outdir, "diamond", "*.blastx.tsv")))
if not files:
    sys.exit(f"No DIAMOND output in {outdir}/diamond, nothing to summarise")
hits = pd.concat([pd.read_csv(f, sep="\t", dtype={"staxids": str})
                    .assign(sample=os.path.basename(f).split(".")[0]) for f in files],
                 ignore_index=True)
best = (hits.sort_values(["bitscore", "length"], ascending=False, kind="stable")
            .drop_duplicates("qseqid")
            .sort_values(["sample", "qseqid"]))

# Read the virus family and genus from staxids, with the db's own taxdump
parent, rank = {}, {}
with open(os.path.join(taxdump, "nodes.dmp")) as f:
    for line in f:
        t, p, r = line.rstrip("\t|\n").split("\t|\t")[:3]
        parent[t], rank[t] = p, r
taxids = set(best["staxids"].fillna("").str.split(";").str[0]) - {"", "0"}
anc = {}
for t in taxids:
    found, node = {}, t
    while node in parent and len(found) < 2:
        if rank[node] in ("family", "genus"):
            found[rank[node]] = node
        if parent[node] == node:
            break
        node = parent[node]
    anc[t] = found
needed = {a for f in anc.values() for a in f.values()}
names = {}
with open(os.path.join(taxdump, "names.dmp")) as f:
    for line in f:
        parts = line.rstrip("\t|\n").split("\t|\t")
        if parts[0] in needed and parts[3] == "scientific name":
            names[parts[0]] = parts[1]
first = best["staxids"].fillna("").str.split(";").str[0]
for r in ("family", "genus"):
    best["virus_" + r] = first.map(lambda t: names.get(anc.get(t, {}).get(r), ""))

best = best.merge(meta[["sample", "bioproject", "country"]], on="sample", how="left")
cols = list(hits.columns)
i = cols.index("sscinames") + 1
cols = cols[:i] + ["virus_family", "virus_genus"] + cols[i:]
cols.insert(cols.index("sample") + 1, "bioproject")
cols.insert(cols.index("bioproject") + 1, "country")
best[cols].to_csv(os.path.join(summ, f"{prefix}_best_hits.tsv"), sep="\t", index=False)

counts = pd.crosstab(best["sample"], best["sscinames"])
counts.to_csv(os.path.join(summ, f"{prefix}_virus_contig_counts.tsv"), sep="\t")
(counts > 0).astype(int).to_csv(os.path.join(summ, f"{prefix}_virus_presence_absence.tsv"), sep="\t")

# Build one timing row per run, from the time -v logs
fields = {"Elapsed (wall clock) time (h:mm:ss or m:ss)": "wall_s",
          "User time (seconds)": "user_s", "System time (seconds)": "system_s",
          "Percent of CPU this job got": "cpu_percent",
          "Maximum resident set size (kbytes)": "max_rss_kb", "Exit status": "exit_status"}
rows = []
for path in sorted(glob.glob(os.path.join(outdir, "timing", "*.time.txt"))):
    row = {"sample": os.path.basename(path).removesuffix(".time.txt")}
    for line in open(path):
        key, _, val = line.strip().rpartition(": ")
        if key == "Command being timed":
            m = re.search(r"--threads (\d+)", val)
            row["threads"] = m.group(1) if m else ""
        elif key in fields:
            row[fields[key]] = val.rstrip("%")
    if "wall_s" in row:
        secs = 0.0
        for part in row["wall_s"].split(":"):
            secs = secs * 60 + float(part)
        row["wall_s"] = round(secs, 2)
    rows.append(row)
timing = pd.DataFrame(rows).merge(meta, on="sample", how="left")
n_hits = hits.groupby("sample").size()
timing["n_hits"] = timing["sample"].map(n_hits).fillna(0).astype(int)
timing = timing.reindex(columns=["sample", "bioproject", "country", "n_contigs", "n_hits", "threads",
                 "wall_s", "user_s", "system_s", "cpu_percent", "max_rss_kb", "exit_status"])
timing.to_csv(os.path.join(summ, f"{prefix}_timing_summary.tsv"), sep="\t", index=False)

print(f"Summaries in {summ}: {len(files)} runs with output, {best['sample'].nunique()} with hits, "
      f"{len(best)} contigs with a best hit, {counts.shape[1]} viruses")
print(f"Total DIAMOND wall time {timing['wall_s'].sum() / 3600:.1f} h")
PYEOF
