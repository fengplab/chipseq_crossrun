#!/usr/bin/env bash
# Step 05: scan every sample with every motif family, and compare with known motifs.
#   work/scan/<s>.peaks.summary.tsv / .background.summary.tsv   prevalence in peaks vs genomic background
#   work/scan/<s>.peaks.spacing.tsv                              inter-site spacing histogram
#   work/scan/<s>.peaks.sites.bed.gz                             genomic motif sites
#   work/scan/<s>.site_labels.tsv                                motif sites falling in repeat/censat labels
#   work/known/tomtom_reps_vs_known.tsv, tomtom_discovered_vs_known.tsv
#   work/known/sea/<s>.sea.tsv                                   SEA of every known motif per sample
#   work/known/known_motif_index.tsv
#   work/motifs/reps_vs_repeat_consensus.tsv                     (if repeat_consensus_fasta given)
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$WORK/scan" "$WORK/known/sea"
REPS="$WORK/motifs/family_reps.meme"
FEAT="$WORK/motifs/family_rep_features.tsv"
[[ -s "$REPS" ]] || { log "run step 04 first"; exit 1; }

scan_sample() {
    set -euo pipefail
    local s="$1" q="$WORK/seqs/$1" o="$WORK/scan/$1"
    if [[ ! -s "$o.background.summary.tsv" ]]; then
        for set in peaks background; do
            local fa="$q.fa" bed="$q.scan.bed"
            [[ $set == background ]] && { fa="$q.bg.fa"; bed="$q.bg.bed"; }
            fimo --text --skip-matched-sequence --verbosity 1 --thresh "$FIMO_THRESH" "$REPS" "$fa" \
                | python3 "$ROOT/bin/summarize_fimo.py" --sample "$s" --set "$set" --scan-bed "$bed" \
                    --families "$FEAT" --spacing-max "$SPACING_MAX" --out "$o.$set"
        done
    fi
    # motif sites (midpoint) in repeat / censat labels
    printf "sample_id\tfamily_id\tfeature_set\tlabel\tn_sites\n" > "$o.site_labels.tsv"
    for fs in repeatmasker censat; do
        local cb="$WORK/annot/$fs.class.bed"
        [[ -s "$cb" ]] || continue
        zcat "$o.peaks.sites.bed.gz" \
            | awk 'BEGIN{FS=OFS="\t"} { m=int(($2+$3)/2); print $1, m, m+1, $4 }' | sort -k1,1 -k2,2n \
            | bedtools intersect -a stdin -b "$cb" -wa -wb -sorted -g "$(genome_file)" \
            | awk -v s="$s" -v fs="$fs" 'BEGIN{FS=OFS="\t"}
                { k=$1 SUBSEP $2 SUBSEP $4 SUBSEP $8; if (!(k in seen)) { seen[k]=1; n[$4 SUBSEP $8]++ }
                  a=$1 SUBSEP $2 SUBSEP $4; if (!(a in any)) { any[a]=1; na[$4]++ } }
                END { for (k in n) { split(k,p,SUBSEP); print s, p[1], fs, p[2], n[k] }
                      for (f in na) print s, f, fs, "ANY", na[f] }' >> "$o.site_labels.tsv"
    done
    # SEA with every known motif (same databases for every sample).
    # --text writes sea.tsv to stdout and skips the HTML report, so a broken/missing MEME HTML template
    # ("FATAL: Template does not contain data section") cannot prevent the results from being written.
    if [[ -n "$KNOWN_DBS" && ! -s "$WORK/known/sea/$s.sea.tsv" ]]; then
        local m_args="" ctrl=""
        for d in $KNOWN_DBS; do m_args="$m_args --m $d"; done
        [[ "$SEA_CONTROL" == "genomic" ]] && ctrl="--n $q.bg.fa"
        if ! sea --text --p "$q.fa" $ctrl $m_args --verbosity 1 > "$WORK/known/sea/$s.sea.tsv.tmp" 2> "$WORK/known/sea/$s.sea.log"; then
            log "ERROR: SEA failed for $s; see $WORK/known/sea/$s.sea.log"
            tail -5 "$WORK/known/sea/$s.sea.log" >&2
            return 1
        fi
        if ! grep -q "^RANK" "$WORK/known/sea/$s.sea.tsv.tmp"; then
            log "ERROR: SEA output for $s has no results header; see $WORK/known/sea/$s.sea.log"
            return 1
        fi
        mv "$WORK/known/sea/$s.sea.tsv.tmp" "$WORK/known/sea/$s.sea.tsv"
    fi
    log "scanned $s"
}
export -f scan_sample
export ROOT REPS FEAT
genome_file >/dev/null
cut -f1 <(sample_col peaks_file) | xargs -P "$THREADS" -I{} bash -c 'scan_sample "$@"' _ {}

if [[ -n "$KNOWN_DBS" ]]; then
    # every sample must have SEA results, otherwise the known-motif tables would be silently incomplete
    missing_sea=$(cut -f1 <(sample_col peaks_file) | while read -r s; do [[ -s "$WORK/known/sea/$s.sea.tsv" ]] || echo "$s"; done)
    if [[ -n "$missing_sea" ]]; then
        log "ERROR: SEA results missing for: $(echo $missing_sea | tr ' ' ',')"
        exit 1
    fi
    python3 "$ROOT/bin/meme_index.py" $KNOWN_DBS > "$WORK/known/known_motif_index.tsv"
    # Tomtom vs known databases is slow; reuse earlier results when the databases have not changed
    stamp="$WORK/known/tomtom_dbs.txt"
    tt1="$WORK/known/tomtom_reps_vs_known.tsv"; tt2="$WORK/known/tomtom_discovered_vs_known.tsv"
    if [[ -s "$tt1" && -s "$tt2" && ! -f "$stamp" ]]; then
        echo "$KNOWN_DBS" > "$stamp"   # results from an earlier version of the bundle: adopt them
    fi
    if [[ -s "$tt1" && -s "$tt2" && "$(cat "$stamp")" == "$KNOWN_DBS" \
          && ! "$WORK/motifs/family_reps.meme" -nt "$tt1" && ! "$WORK/motifs/all_discovered.meme" -nt "$tt2" ]]; then
        log "Tomtom vs known databases: reusing existing results (delete $WORK/known/tomtom_*.tsv to force a rerun)"
    else
        log "Tomtom: family representatives and all discovered motifs vs known databases"
        tomtom -no-ssc -text -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 1 \
            "$REPS" $KNOWN_DBS > "$tt1.tmp" && mv "$tt1.tmp" "$tt1"
        tomtom -no-ssc -text -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 1 \
            "$WORK/motifs/all_discovered.meme" $KNOWN_DBS > "$tt2.tmp" && mv "$tt2.tmp" "$tt2"
        echo "$KNOWN_DBS" > "$stamp"
    fi
else
    log "no known_motif_dbs configured: known-motif comparisons skipped"
fi

if [[ -n "$REPEAT_CONSENSUS_FASTA" ]]; then
    log "scanning repeat consensus sequences with family representatives"
    fimo --text --skip-matched-sequence --verbosity 1 --thresh "$FIMO_THRESH" "$REPS" "$REPEAT_CONSENSUS_FASTA" \
        > "$WORK/motifs/reps_vs_repeat_consensus.tsv"
fi
