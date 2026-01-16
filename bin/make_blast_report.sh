#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 SAMPLE_NAME" >&2
    exit 1
fi

SampleName="$1"

blast_file="${SampleName}_blast.tsv"
report_file="${SampleName}_report_blast.tsv"
best_hits_file="${SampleName}_report_blast_best.tsv"

if [[ ! -f "$blast_file" ]]; then
    echo "Error: BLAST file '$blast_file' not found" >&2
    exit 1
fi

# Always write header first
echo -e "queryid\tsubject_id\talignment length\tquery_coverage\t%identity\tevalue\tstaxids\tsscinames\tscomnames\tstitle" > "$report_file"

# Append all non-comment hit lines (if any)
awk 'substr($1,1,1) != "#" && NF > 0 { print }' "$blast_file" | sort | uniq >> "$report_file"

# Add one "none" row for each query that had no hits
awk '
    /^# Query:/       { q = $3 }
    /^# 0 hits found/ { nohit[q] = 1 }
    END {
        for (q in nohit)
            print q "\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone";
    }
' "$blast_file" >> "$report_file"

# If the report still has only the header (no hits and no per-query "none" rows),
# add a single global "none" row
if [[ "$(wc -l < "$report_file")" -le 1 ]]; then
    echo -e "none\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone\tnone" >> "$report_file"
fi

# Keep only the best hit per query (highest coverage, then highest identity)
# 1) write header
head -n 1 "$report_file" > "$best_hits_file"
# 2) sort by query, coverage (col 4 desc), identity (col 5 desc)
# 3) keep first row per query
#    (we assume input is already sorted by those keys)
tail -n +2 "$report_file" \
    | sort -t$'\t' -k1,1 -k4,4rn -k5,5rn \
    | awk -F'\t' 'BEGIN{OFS=FS} { if ($1 != prev) { print; prev = $1 } }' >> "$best_hits_file"
