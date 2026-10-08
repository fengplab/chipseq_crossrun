#!/usr/bin/env bash
# Step 05: scan every sample with every motif family, and compare with known motifs.
#   work/scan/<s>.peaks.summary.tsv / .background.summary.tsv   prevalence in peaks vs genomic background
#   work/scan/<s>.peaks.spacing.tsv                              inter-site spacing histogram
#   work/scan/<s>.peaks.sites.bed.gz                             genomic motif sites
#   work/scan/<s>.site_labels.tsv                                motif sites falling in repeat/censat labels
#   work/scan/<s>.family_sea.tsv                                 SEA of every family vs shuffled peak sequences
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
    # rescan when the families or the sequences changed, or the results predate the current scan format
    if [[ ! -s "$o.background.summary.tsv" || "$REPS" -nt "$o.background.summary.tsv" || "$q.fa" -nt "$o.background.summary.tsv"
          || "$(cat "$o.scan_version" 2>/dev/null)" != "$SCAN_VERSION" ]]; then
        local pgc; pgc=$(pgc_flag fimo)
        for set in peaks background; do
            local fa="$q.fa" bed="$q.scan.bed"
            [[ $set == background ]] && { fa="$q.bg.fa"; bed="$q.bg.bed"; }
            # shellcheck disable=SC2086
            $MEME_PREFIX fimo --text --skip-matched-sequence $pgc --verbosity 1 --thresh "$FIMO_THRESH" "$REPS" "$fa" \
                | python3 "$ROOT/bin/summarize_fimo.py" --sample "$s" --set "$set" --scan-bed "$bed" \
                    --families "$FEAT" --spacing-max "$SPACING_MAX" --out "$o.$set"
        done
        printf "%s\n" "$SCAN_VERSION" > "$o.scan_version"
    fi
    # Family representatives vs shuffled copies of the same peaks (SEA's default control, the same kind of
    # control XSTREME used to discover the motifs). --thresh 1e6 reports every family, significant or not.
    if [[ ! -s "$o.family_sea.tsv" || "$REPS" -nt "$o.family_sea.tsv" || "$q.fa" -nt "$o.family_sea.tsv" ]]; then
        # shellcheck disable=SC2046
        if ! $MEME_PREFIX sea --text $(pgc_flag sea) --p "$q.fa" --m "$REPS" --thresh 1e6 --verbosity 1 > "$o.family_sea.tsv.tmp" 2> "$o.family_sea.log"; then
            log "ERROR: SEA (motif families) failed for $s; see $o.family_sea.log"; tail -5 "$o.family_sea.log" >&2; return 1
        fi
        grep -q "^RANK" "$o.family_sea.tsv.tmp" || { log "ERROR: SEA (motif families) output for $s has no header"; return 1; }
        mv "$o.family_sea.tsv.tmp" "$o.family_sea.tsv"
    fi
    # motif sites (midpoint) in repeat / censat labels
    printf "sample_id\tfamily_id\tfeature_set\tlabel\tn_sites\n" > "$o.site_labels.tsv"
    for fs in repeatmasker censat; do
        local cb="$WORK/annot/$fs.class.bed"
        [[ -s "$cb" ]] || continue
        zcat "$o.peaks.sites.bed.gz" \
            | awk 'BEGIN{FS=OFS="\t"} { m=int(($2+$3)/2); print $1, m, m+1, $4 }' | sort -k1,1 -k2,2n \
            | isect_sorted -b "$cb" -wa -wb -sorted -g "$(genome_file)" \
            | awk -v s="$s" -v fs="$fs" 'BEGIN{FS=OFS="\t"}
                { k=$1 SUBSEP $2 SUBSEP $4 SUBSEP $8; if (!(k in seen)) { seen[k]=1; n[$4 SUBSEP $8]++ }
                  a=$1 SUBSEP $2 SUBSEP $4; if (!(a in any)) { any[a]=1; na[$4]++ } }
                END { for (k in n) { split(k,p,SUBSEP); print s, p[1], fs, p[2], n[k] }
                      for (f in na) print s, f, fs, "ANY", na[f] }' >> "$o.site_labels.tsv"
    done
    # SEA with every known motif (same databases for every sample).
    # --text writes sea.tsv to stdout and skips the HTML report, so a broken/missing MEME HTML template
    # ("FATAL: Template does not contain data section") cannot prevent the results from being written.
    if [[ -n "$KNOWN_DBS" && ( ! -s "$WORK/known/sea/$s.sea.tsv" || "$q.fa" -nt "$WORK/known/sea/$s.sea.tsv" ) ]]; then
        local m_args="" ctrl=""
        for d in $KNOWN_DBS; do m_args="$m_args --m $d"; done
        [[ "$SEA_CONTROL" == "genomic" ]] && ctrl="--n $q.bg.fa"
        # shellcheck disable=SC2046
        if ! $MEME_PREFIX sea --text $(pgc_flag sea) --p "$q.fa" $ctrl $m_args --verbosity 1 > "$WORK/known/sea/$s.sea.tsv.tmp" 2> "$WORK/known/sea/$s.sea.log"; then
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
# Known-motif results are cached; invalidate them when the databases change (path or content)
DB_SIG=""
if [[ -n "$KNOWN_DBS" ]]; then
    DB_SIG=$(for d in $KNOWN_DBS; do printf "%s %s\n" "$(md5sum < "$d" | cut -d' ' -f1)" "$d"; done)
    if [[ -f "$WORK/known/sea/dbs.sig" && "$(cat "$WORK/known/sea/dbs.sig")" != "$DB_SIG" ]]; then
        log "known-motif databases changed: recomputing SEA known-motif results for every sample"
        rm -f "$WORK"/known/sea/*.sea.tsv
    fi
    printf "%s\n" "$DB_SIG" > "$WORK/known/sea/dbs.sig"
fi
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
    stamp="$WORK/known/tomtom_dbs.sig"
    tt1="$WORK/known/tomtom_reps_vs_known.tsv"; tt2="$WORK/known/tomtom_discovered_vs_known.tsv"
    if [[ -s "$tt1" && -s "$tt2" && ! -f "$stamp" && -f "$WORK/known/tomtom_dbs.txt" \
          && "$(cat "$WORK/known/tomtom_dbs.txt")" == "$KNOWN_DBS" ]]; then
        printf "%s\n" "$DB_SIG" > "$stamp"   # results from an earlier version of the bundle with the same files: adopt them
    fi
    if [[ -s "$tt1" && -s "$tt2" && -f "$stamp" && "$(cat "$stamp")" == "$DB_SIG" \
          && ! "$WORK/motifs/family_reps.meme" -nt "$tt1" && ! "$WORK/motifs/all_discovered.meme" -nt "$tt2" ]]; then
        log "Tomtom vs known databases: reusing existing results (delete $WORK/known/tomtom_*.tsv to force a rerun)"
    else
        log "Tomtom: family representatives and all discovered motifs vs known databases"
        $MEME_PREFIX tomtom -no-ssc -text -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 1 \
            "$REPS" $KNOWN_DBS > "$tt1.tmp" && mv "$tt1.tmp" "$tt1"
        $MEME_PREFIX tomtom -no-ssc -text -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 1 \
            "$WORK/motifs/all_discovered.meme" $KNOWN_DBS > "$tt2.tmp" && mv "$tt2.tmp" "$tt2"
        printf "%s\n" "$DB_SIG" > "$stamp"
    fi
else
    log "no known_motif_dbs configured: known-motif comparisons skipped"
fi

if [[ -n "$REPEAT_CONSENSUS_FASTA" ]]; then
    log "scanning repeat consensus sequences with family representatives"
    # shellcheck disable=SC2046
    $MEME_PREFIX fimo --text --skip-matched-sequence $(pgc_flag fimo) --verbosity 1 --thresh "$FIMO_THRESH" "$REPS" "$REPEAT_CONSENSUS_FASTA" \
        > "$WORK/motifs/reps_vs_repeat_consensus.tsv"
fi
