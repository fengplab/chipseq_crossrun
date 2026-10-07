#!/usr/bin/env bash
# Step 16: per-window features for the regression models (step 17).
#   work/models/win/<s>.windows.tsv.gz   one row per window: every peak passing models.peak_qvalue, plus
#                                        MODEL_BG_DRAWS GC- and width-matched background windows for the peaks
#                                        that enter the peak-level model (at most models.max_peaks_per_sample)
#                                        columns: width, GC, CenSat / RepeatMasker class (largest overlap),
#                                        merged region, IP and control read counts
#   work/models/win/<s>.sites.tsv.gz     window x motif family: sites, sites in tandem arrays
#   work/models/libsize.tsv              total reads per IP and control library
#   work/reads/reads_in_labels.tsv       IP / control reads per repeat class (step 07 is run here if it was skipped)
# Reused from earlier steps: scanned windows (03), family representatives (04), motif sites in peaks (05),
# class BEDs (01), merged regions (08).
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[[ -n "${MODELS_ENABLED:-}" ]] || { log "ERROR: model settings missing from work/config.env; rerun step 0 (./run_all.sh config.yaml 0 0)"; exit 1; }
[[ "$MODELS_ENABLED" == "1" ]] || { log "models disabled (models.enabled: false); skipping"; exit 0; }
M="$WORK/models"; mkdir -p "$M/win" "$M/tmp"
REPS="$WORK/motifs/family_reps.meme"
[[ -s "$REPS" ]] || { log "run step 04 first"; exit 1; }
[[ -s "$WORK/regions/merged_regions.tsv" ]] || { log "run step 08 first"; exit 1; }
for t in bedtools samtools featureCounts; do
    command -v "$t" >/dev/null || { log "ERROR: $t not found on PATH (needed for the model features)"; exit 1; }
done
genome_file >/dev/null
awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1,$2,$3,$4}' "$WORK/regions/merged_regions.tsv" | sort -k1,1 -k2,2n > "$M/regions.bed"
export ROOT M REPS
# settings that decide which windows exist; a change invalidates the per-sample results
SIG="q=$MODEL_PEAK_Q;cap=$MODEL_MAX_PEAKS;draws=$MODEL_BG_DRAWS;seed=$MODEL_SEED;gc=$GC_CANDIDATES;fimo=$FIMO_THRESH;sp=$SPACING_MAX;reps=$(md5sum < "$REPS" | cut -c1-12)"
export SIG

# reads per window for one BAM (paired-end libraries count fragments). SAF names follow the BAM header.
count_bam() {
    local bam="$1" bed="$2" out="$3" t; t=$(mktemp -d "$M/tmp/fc.XXXXXX")
    samtools quickcheck "$bam" || { log "ERROR: cannot read $bam"; return 1; }
    samtools view -H "$bam" | awk -F'\t' '$1=="@SQ"{sub(/^SN:/,"",$2); print $2}' > "$t/chroms"
    awk 'BEGIN{FS=OFS="\t"} FILENAME==ARGV[1] { a[$1]=$2; next }
         FILENAME==ARGV[2] { c = ($1 in a) ? a[$1] : $1; bam[c]=$1; next }
         FNR==1 { print "GeneID","Chr","Start","End","Strand" }
         ($1 in bam) { print $4, bam[$1], $2+1, $3, "+" }' "$WORK/annot/chrom_alias.tsv" "$t/chroms" "$bed" > "$t/w.saf"
    [[ $(wc -l < "$t/w.saf") -gt 1 ]] || { log "ERROR: no windows match the chromosome names of $bam"; return 1; }
    local extra=""
    [[ $(samtools view -f 1 "$bam" | head -n 1000 | wc -l) -gt 0 ]] && extra="-p --countReadPairs"
    # -O: a read overlapping two windows (a peak and a background window, say) counts for both
    featureCounts -a "$t/w.saf" -F SAF -O -T 1 --tmpDir "$t" $extra -o "$t/counts" "$bam" > "$t/log" 2>&1 \
        || { log "ERROR: featureCounts failed for $bam"; tail -5 "$t/log" >&2; return 1; }
    mv "$t/counts" "$out"; mv "$t/counts.summary" "$out.summary"
    rm -rf "$t"
}

features_one() {
    set -euo pipefail
    local s="$1" o="$M/win/$1" q="$WORK/seqs/$1"
    if [[ -s "$o.windows.tsv.gz" && -s "$o.sites.tsv.gz" && -s "$o.sig" && "$(cat "$o.sig")" == "$SIG" ]]; then return 0; fi
    local t; t=$(mktemp -d "$M/tmp/$s.XXXXXX")
    local qarg=""; [[ -n "$MODEL_PEAK_Q" ]] && qarg="--max-q $MODEL_PEAK_Q"
    # shellcheck disable=SC2086
    python3 "$ROOT/bin/model_windows.py" select --sample "$s" --scan-bed "$q.scan.bed" \
        --peaks-file "$(sample_get "$s" peaks_file)" --out "$t/sel" $qarg --max-peaks "$MODEL_MAX_PEAKS" --seed "$MODEL_SEED"
    [[ -s "$t/sel.model.bed" ]] || { log "WARNING: $s has no peaks passing the q-value cutoff; skipped"; rm -rf "$t"; return 0; }
    # GC- and width-matched background windows (always GC-matched here, whatever repeats.null_model says)
    : > "$t/bg.bed"
    local k
    for k in $(seq 1 "$MODEL_BG_DRAWS"); do
        NULL_MODEL=gc_matched shuffle_once "$t/sel.model.bed" $(( MODEL_SEED * 100 + k )) \
            | awk -v k="$k" 'BEGIN{FS=OFS="\t"} { print $1, $2, $3, $4 "#bg" k }' >> "$t/bg.bed"
    done
    { cat "$t/sel.peaks.bed"; awk 'BEGIN{FS=OFS="\t"} {print $0, "NA"}' "$t/bg.bed"; } | sort -k1,1 -k2,2n > "$t/windows.bed"
    cut -f1-4 "$t/windows.bed" > "$t/w4.bed"
    to_fasta_names < "$t/w4.bed" | bedtools nuc -fi "$GENOME_FASTA" -bed stdin > "$t/nuc.tsv"
    local fs args=""
    for fs in censat repeatmasker; do
        if [[ -s "$WORK/annot/$fs.class.bed" ]]; then
            bedtools intersect -a "$t/w4.bed" -b "$WORK/annot/$fs.class.bed" -wo > "$t/$fs.ov"
            args="$args --$fs $t/$fs.ov"
        fi
    done
    cut -f1-4 "$t/sel.peaks.bed" | bedtools intersect -a stdin -b "$M/regions.bed" -wa -wb > "$t/regions.ov"
    # motif sites in the background windows (sites in peaks come from step 05)
    sort -k1,1 -k2,2n "$t/bg.bed" > "$t/bg.sorted.bed"
    getfasta_canonical "$t/bg.sorted.bed" "$t/bg.fa"
    $MEME_PREFIX fimo --text --skip-matched-sequence --verbosity 1 --thresh "$FIMO_THRESH" "$REPS" "$t/bg.fa" > "$t/bg.fimo.tsv" \
        || { log "ERROR: FIMO failed on the background windows of $s"; return 1; }
    local bam ctl
    bam=$(sample_get "$s" bam_file); ctl=$(sample_get "$s" control_bam)
    if [[ -n "$bam" ]]; then count_bam "$bam" "$t/w4.bed" "$t/ip.counts"; args="$args --ip-counts $t/ip.counts"
    else log "WARNING: $s has no BAM; IP/input read ratios are not available for it"; fi
    if [[ -n "$ctl" ]]; then count_bam "$ctl" "$t/w4.bed" "$t/ctl.counts"; args="$args --ctl-counts $t/ctl.counts"
    else log "WARNING: $s has no control BAM; IP/input read ratios are not available for it"; fi
    # shellcheck disable=SC2086
    python3 "$ROOT/bin/model_windows.py" assemble --sample "$s" --windows "$t/windows.bed" --model-bed "$t/sel.model.bed" \
        --nuc "$t/nuc.tsv" --regions "$t/regions.ov" $args --peak-sites "$WORK/scan/$s.peaks.sites.bed.gz" \
        --bg-fimo "$t/bg.fimo.tsv" --spacing-max "$SPACING_MAX" --out "$t/out"
    mv "$t/out.windows.tsv.gz" "$o.windows.tsv.gz"; mv "$t/out.sites.tsv.gz" "$o.sites.tsv.gz"
    mv "$t/out.libsize.tsv" "$o.libsize.tsv"
    printf "%s\n" "$SIG" > "$o.sig"
    rm -rf "$t"
    log "model features: $s"
}
export -f features_one count_bam
sample_col peaks_file | cut -f1 | xargs -P "$THREADS" -I{} bash -c 'features_one "$@"' _ {}
{ printf "sample_id\tip_total_reads\tctl_total_reads\n"; cat "$M"/win/*.libsize.tsv; } > "$M/libsize.tsv"
rm -rf "$M/tmp"

# IP / control reads per repeat class: reuse step 07's result, or run it now
if [[ ! -s "$WORK/reads/reads_in_labels.tsv" || $(wc -l < "$WORK/reads/reads_in_labels.tsv") -lt 2 ]]; then
    log "reads per repeat class not available yet; running step 07"
    sed 's/^export READS_IN_REPEATS=.*/export READS_IN_REPEATS=1/' "$1" > "$M/config.reads.env"
    bash "$ROOT/steps/07_reads_in_repeats.sh" "$M/config.reads.env"
fi
log "model features written for $(ls "$M"/win/*.windows.tsv.gz 2>/dev/null | wc -l) samples"
