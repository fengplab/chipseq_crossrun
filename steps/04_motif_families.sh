#!/usr/bin/env bash
# Step 04: cross-run motif families.
#   work/motifs/tomtom_all_vs_all.tsv, families.tsv, family_reps.meme, family_rep_features.tsv
# Two stages: (1) Tomtom PWM similarity (R/04_cluster_families.R); (2) optional, motifs.family_cooccupancy_min:
# Tomtom families whose sites sit on the same positions of the peak sequences are merged (bin/cooccupancy.py).
#   work/motifs/families.tomtom.tsv        families after stage 1 only (kept when stage 2 runs)
#   work/motifs/cooccupancy/               stage-2 FIMO scans and family_cooccupancy_pairs.tsv
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
M="$WORK/motifs"
[[ -s "$M/all_discovered.meme" ]] || { log "no de novo motifs; run step 02 first"; exit 1; }
if [[ ! -s "$M/tomtom_all_vs_all.tsv" || "$M/all_discovered.meme" -nt "$M/tomtom_all_vs_all.tsv" ]]; then
    log "Tomtom all-vs-all on $(grep -c '^MOTIF' "$M/all_discovered.meme") discovered motifs"
    $MEME_PREFIX tomtom -no-ssc -text -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 1 \
        "$M/all_discovered.meme" "$M/all_discovered.meme" > "$M/tomtom_all_vs_all.tsv"
fi
Rscript "$ROOT/R/04_cluster_families.R" "$CONFIG_YAML"
python3 "$ROOT/bin/write_family_reps.py" "$WORK"

if [[ -n "${FAMILY_COOCC_MIN:-}" ]]; then
    C="$M/cooccupancy"; mkdir -p "$C"
    cp "$M/families.tsv" "$M/families.tomtom.tsv"
    cp "$M/family_reps.meme" "$C/tomtom_reps.meme"
    pgc=$(pgc_flag fimo)
    fimo_files=()
    while read -r s _; do
        fa="$WORK/seqs/$s.fa"
        [[ -s "$fa" ]] || { log "ERROR: $fa missing; step 03 writes it and the co-occupancy stage needs it"; exit 1; }
        # same scan as step 05 (FIMO --text: start/stop 1-based closed; cooccupancy.py converts to 0-based half-open)
        $MEME_PREFIX fimo --text --skip-matched-sequence $pgc --verbosity 1 --thresh "$FIMO_THRESH" \
            "$C/tomtom_reps.meme" "$fa" > "$C/$s.fimo.tsv.tmp"
        mv "$C/$s.fimo.tsv.tmp" "$C/$s.fimo.tsv"
        fimo_files+=("$C/$s.fimo.tsv")
    done < <(sample_col sample_id)
    [[ ${#fimo_files[@]} -gt 0 ]] || { log "ERROR: no samples to scan for the co-occupancy stage"; exit 1; }
    log "co-occupancy stage: merging Tomtom families with site-co-occupancy similarity >= $FAMILY_COOCC_MIN (containment rule: ${FAMILY_COOCC_CONTAINMENT}, 0 = off; floor ${FAMILY_COOCC_FLOOR})"
    python3 "$ROOT/bin/cooccupancy.py" merge --families "$M/families.tomtom.tsv" --fimo "${fimo_files[@]}" \
        --threshold "$FAMILY_COOCC_MIN" --overlap "$FAMILY_COOCC_OVERLAP" \
        --containment "$FAMILY_COOCC_CONTAINMENT" --floor "$FAMILY_COOCC_FLOOR" \
        --pairs-out "$C/family_cooccupancy_pairs.tsv" --out "$M/families.tsv"
    python3 "$ROOT/bin/write_family_reps.py" "$WORK"
fi
