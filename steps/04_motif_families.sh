#!/usr/bin/env bash
# Step 04: cross-run motif families.
#   work/motifs/tomtom_all_vs_all.tsv, families.tsv, family_reps.meme, family_rep_features.tsv
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
M="$WORK/motifs"
[[ -s "$M/all_discovered.meme" ]] || { log "no de novo motifs; run step 02 first"; exit 1; }
if [[ ! -s "$M/tomtom_all_vs_all.tsv" ]]; then
    log "Tomtom all-vs-all on $(grep -c '^MOTIF' "$M/all_discovered.meme") discovered motifs"
    tomtom -no-ssc -text -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 1 \
        "$M/all_discovered.meme" "$M/all_discovered.meme" > "$M/tomtom_all_vs_all.tsv"
fi
Rscript "$ROOT/R/04_cluster_families.R" "$CONFIG_YAML"
python3 "$ROOT/bin/write_family_reps.py" "$WORK"
