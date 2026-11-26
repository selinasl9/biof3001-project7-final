#!/usr/bin/env bash
#
# assemble_runs_without_contigs.sh
#
# Assemble the runs Logan has no contigs for, so they are screened rather than dropped.
#
# Contigs are written in the Logan layout, so the DIAMOND step cannot tell the two apart.
#
# Inputs:
#   RUN_LIST   optional file of runs, one per line, default the last 8 of the run table
#   runs/logan_missing_runs.tsv   the runs with no Logan contigs
#
# Outputs:
#   results/logan_missing/assembly/RUN.contigs.fa   assembled contigs, Logan layout
#   results/logan_missing/diamond/RUN.blastx.tsv    one hit table per run
#   results/logan_missing/summary/                  best hits, virus matrices, timing
#
set -u

# Set the job sizes, the input table + the working folders
# Check nproc and free -g before raising these
THREADS=16
MEM_GB=200

RUN_TABLE="runs/logan_missing_runs.tsv"
OUTDIR="results/logan_missing"
WORK="data/logan_missing"      # Reads + SPAdes working dirs, large and temporary
PREFIX="MISSING"

# Take whichever DIAMOND db and taxdump exist in this clone
DB=""; TAXDUMP=""
for d in reference/arbo_tax.dmnd db/arbovirus_refseq.dmnd; do [ -s "$d" ] && { DB="$d"; break; }; done
for t in reference/taxonomy db/taxonomy; do [ -s "$t/nodes.dmp" ] && { TAXDUMP="$t"; break; }; done
[ -n "$DB" ] || { echo "No DIAMOND db found (looked for reference/arbo_tax.dmnd, db/arbovirus_refseq.dmnd); set DB= in the script"; exit 1; }
[ -n "$TAXDUMP" ] || { echo "No taxdump found (looked for reference/taxonomy, db/taxonomy); set TAXDUMP= in the script"; exit 1; }

# Stop with a clear message if any tool is not installed
missing=""
for tool in prefetch fasterq-dump fastp spades.py diamond python3; do
    command -v "$tool" > /dev/null || missing="$missing $tool"
done
if [ -n "$missing" ]; then
    echo "Missing tools:$missing"
    echo "Install into the env with:  conda install -c conda-forge -c bioconda sra-tools fastp spades"
    exit 1
fi

mkdir -p "$OUTDIR"/{fastp,assembly,diamond,timing,summary} "$WORK"

RUN_LIST="$OUTDIR/runs.txt"
if [ $# -ge 1 ]; then
    grep -v '^\s*$' "$1" > "$RUN_LIST"
else
    tail -n +2 "$RUN_TABLE" | tail -n 8 | cut -f1 > "$RUN_LIST"
fi

STATS="$OUTDIR/assembly_stats.tsv"
[ -s "$STATS" ] || printf "sample\tlayout\tspades_mode\treads_before_fastp\treads_after_fastp\tn_contigs\ttotal_bp\tmax_len\tn50\n" > "$STATS"

echo "$(wc -l < "$RUN_LIST") runs: $(tr '\n' ' ' < "$RUN_LIST")"
echo "DB=$DB  TAXDUMP=$TAXDUMP  THREADS=$THREADS  MEM=${MEM_GB}G"
echo "fastp $(fastp --version 2>&1 | tail -1) | $(spades.py --version 2>&1 | tail -1) | $(diamond version 2>&1 | tail -1)"
date

# Write one assembly_stats row from the fastp json + the renamed contigs
write_stats() {
    acc="$1"; layout="$2"; mode="$3"
    python3 - "$acc" "$layout" "$mode" "$OUTDIR/fastp/$acc.fastp.json" "$OUTDIR/assembly/$acc.contigs.fa" >> "$STATS" <<'PYEOF'
import json, sys
acc, layout, mode, js, fa = sys.argv[1:6]
s = json.load(open(js))["summary"]
lens, n = [], 0
for line in open(fa):
    if line.startswith(">"):
        lens.append(0)
    else:
        lens[-1] += len(line.strip())
lens.sort(reverse=True)
total, run, n50 = sum(lens), 0, 0
for l in lens:
    run += l
    if run >= total / 2:
        n50 = l
        break
print("\t".join(map(str, [acc, layout, mode, s["before_filtering"]["total_reads"],
      s["after_filtering"]["total_reads"], len(lens), total, lens[0] if lens else 0, n50])))
PYEOF
}

for acc in $(cat "$RUN_LIST"); do
    echo "==================== $acc  $(date +%H:%M)"
    CONTIGS="$OUTDIR/assembly/$acc.contigs.fa"
    OUT="$OUTDIR/diamond/$acc.blastx.tsv"
    R="$WORK/$acc"
    mkdir -p "$R"

    # Step 1: download the reads, trim them with fastp + assemble with SPAdes
    if [ ! -s "$CONTIGS" ]; then
        if ! ls "$R"/${acc}*.fastq > /dev/null 2>&1; then
            echo "Downloading reads"
            prefetch --max-size u -O "$WORK" "$acc" < /dev/null \
              && fasterq-dump --split-3 --threads "$THREADS" -O "$R" -t "$R/tmp" "$WORK/$acc/$acc.sra" < /dev/null \
              || fasterq-dump --split-3 --threads "$THREADS" -O "$R" -t "$R/tmp" "$acc" < /dev/null \
              || { echo "FAILED download: $acc"; continue; }
            rm -rf "$R/tmp" "$WORK/$acc/$acc.sra"
        fi

        # Pick the layout, as --split-3 gives _1/_2 for paired data and one file for single
        if [ -s "$R/${acc}_1.fastq" ] && [ -s "$R/${acc}_2.fastq" ]; then
            layout=PAIRED; mode=meta
            echo "Paired-end: running fastp"
            fastp -i "$R/${acc}_1.fastq" -I "$R/${acc}_2.fastq" \
                  -o "$R/${acc}_1.trim.fq.gz" -O "$R/${acc}_2.trim.fq.gz" \
                  -q 30 -l 50 --trim_poly_g --detect_adapter_for_pe \
                  -w "$THREADS" \
                  -j "$OUTDIR/fastp/$acc.fastp.json" -h "$OUTDIR/fastp/$acc.fastp.html" \
                  2> "$OUTDIR/fastp/$acc.fastp.log" || { echo "FAILED fastp: $acc"; continue; }
            echo "Running SPAdes --meta"
            spades.py --meta -1 "$R/${acc}_1.trim.fq.gz" -2 "$R/${acc}_2.trim.fq.gz" \
                      -o "$R/spades" -t "$THREADS" -m "$MEM_GB" > /dev/null 2>&1
        elif [ -s "$R/$acc.fastq" ] || [ -s "$R/${acc}_1.fastq" ]; then
            layout=SINGLE; mode=default_single_end
            SE="$R/$acc.fastq"; [ -s "$SE" ] || SE="$R/${acc}_1.fastq"
            echo "Single-end: running fastp"
            fastp -i "$SE" -o "$R/$acc.trim.fq.gz" \
                  -q 30 -l 50 --trim_poly_g -w "$THREADS" \
                  -j "$OUTDIR/fastp/$acc.fastp.json" -h "$OUTDIR/fastp/$acc.fastp.html" \
                  2> "$OUTDIR/fastp/$acc.fastp.log" || { echo "FAILED fastp: $acc"; continue; }
            echo "NOTE: metaSPAdes does not support single-end reads; running SPAdes without --meta"
            spades.py -s "$R/$acc.trim.fq.gz" \
                      -o "$R/spades" -t "$THREADS" -m "$MEM_GB" > /dev/null 2>&1
        else
            echo "FAILED: no FASTQ found for $acc in $R"; continue
        fi

        cp "$R/spades/spades.log" "$OUTDIR/assembly/$acc.spades.log" 2>/dev/null
        if [ ! -s "$R/spades/contigs.fasta" ]; then
            echo "FAILED SPAdes: $acc (see $OUTDIR/assembly/$acc.spades.log)"
            continue
        fi
        # Rename NODE_12_length_500_cov_3.2 to ACC_12, like Logan contig ids
        awk -v a="$acc" '/^>/ {split(substr($0,2),f,"_"); print ">" a "_" f[2]; next} {print}' \
            "$R/spades/contigs.fasta" > "$CONTIGS.tmp" && mv "$CONTIGS.tmp" "$CONTIGS"
        write_stats "$acc" "$layout" "$mode"
        echo "Assembled: $(grep -c '^>' "$CONTIGS") contigs"
        rm -rf "$R" "$WORK/$acc"
    else
        echo "Contigs already assembled"
    fi

    # Step 2: DIAMOND the contigs, same settings as the Logan ones
    if [ -s "$OUT" ]; then
        echo "DIAMOND already done"
        continue
    fi
    if /usr/bin/time -v -o "$OUTDIR/timing/$acc.time.txt" \
        diamond blastx --query "$CONTIGS" --db "$DB" --out "$OUT.tmp" \
        --threads "$THREADS" --sensitive --evalue 1e-5 --max-target-seqs 5 \
        --header simple --outfmt 6 \
        qseqid sseqid pident length qlen slen qstart qend sstart send \
        evalue bitscore qcovhsp scovhsp staxids sscinames stitle \
        --quiet < /dev/null
    then
        mv "$OUT.tmp" "$OUT"
        echo "DONE: $acc  $(( $(wc -l < "$OUT") - 1 )) hits"
    else
        echo "FAILED DIAMOND: $acc"
        rm -f "$OUT.tmp"
    fi
done

echo "All runs processed"
date

# Step 3: write the summary tables, same format as run_diamond_blastx.sh
python3 - "$OUTDIR" "$RUN_TABLE" "$TAXDUMP" "$PREFIX" "$STATS" <<'PYEOF'
import glob, os, re, sys
import pandas as pd

outdir, run_table, taxdump, prefix, stats = sys.argv[1:6]
summ = os.path.join(outdir, "summary")

meta = (pd.read_csv(run_table, sep="\t", dtype=str).fillna("")
          [["acc", "bioproject", "geo_loc_name_country_calc"]]
          .drop_duplicates("acc")
          .rename(columns={"acc": "sample", "geo_loc_name_country_calc": "country"}))
asm = pd.read_csv(stats, sep="\t", dtype=str).drop_duplicates("sample", keep="last")
meta = meta.merge(asm[["sample", "n_contigs"]], on="sample", how="left")

files = sorted(glob.glob(os.path.join(outdir, "diamond", "*.blastx.tsv")))
if not files:
    sys.exit(f"No DIAMOND output in {outdir}/diamond, nothing to summarise")
hits = pd.concat([pd.read_csv(f, sep="\t", dtype={"staxids": str})
                    .assign(sample=os.path.basename(f).split(".")[0]) for f in files],
                 ignore_index=True)
best = (hits.sort_values(["bitscore", "length"], ascending=False, kind="stable")
            .drop_duplicates("qseqid")
            .sort_values(["sample", "qseqid"]))

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
PYEOF
