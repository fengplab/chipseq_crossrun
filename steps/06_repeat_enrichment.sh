#!/usr/bin/env bash
# Step 06: observed vs expected peak counts per repeat / CenSat label.
#   work/repeats/<s>.<featureset>.counts.tsv   sample_id feature_set iter label n   (iter 0 = observed)
#   work/repeats/<s>.<featureset>.peak_labels.tsv  peak_id chrom label   (observed peaks)
# The same overlap rule (config repeats.overlap_rule) is applied to observed and shuffled peaks.
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
mkdir -p "$WORK/repeats"

enrich_sample() {
    set -euo pipefail
    local s="$1" b="$WORK/peaks/$1.bed5" tmp
    tmp=$(mktemp -d)
    cut -f1-4 "$b" > "$tmp/obs.bed"
    for fs in repeatmasker censat; do
        local cb="$WORK/annot/$fs.class.bed" out="$WORK/repeats/$s.$fs.counts.tsv"
        [[ -s "$cb" ]] || continue
        [[ -s "$out" ]] && continue
        label_peaks "$tmp/obs.bed" "$cb" > "$WORK/repeats/$s.$fs.peak_labels.tsv"
        count_labels "$tmp/obs.bed" "$cb" | awk -v s="$s" -v fs="$fs" 'BEGIN{OFS="\t"} {print s, fs, 0, $1, $2}' > "$out.part"
    done
    for i in $(seq 1 "$N_SHUFFLES"); do
        shuffle_once "$tmp/obs.bed" "$i" > "$tmp/shuf.bed"
        for fs in repeatmasker censat; do
            local cb="$WORK/annot/$fs.class.bed" out="$WORK/repeats/$s.$fs.counts.tsv"
            [[ -s "$out.part" ]] || continue
            count_labels "$tmp/shuf.bed" "$cb" | awk -v s="$s" -v fs="$fs" -v i="$i" 'BEGIN{OFS="\t"} {print s, fs, i, $1, $2}' >> "$out.part"
        done
    done
    for fs in repeatmasker censat; do
        local out="$WORK/repeats/$s.$fs.counts.tsv"
        [[ -s "$out.part" ]] && { printf "sample_id\tfeature_set\titer\tlabel\tn\n"; cat "$out.part"; } > "$out" && rm -f "$out.part"
    done
    rm -rf "$tmp"
    log "repeat enrichment: $s ($N_SHUFFLES shuffles, $NULL_MODEL null, $OVERLAP_RULE rule)"
}
export -f enrich_sample
genome_file >/dev/null
cut -f1 <(sample_col peaks_file) | xargs -P "$THREADS" -I{} bash -c 'enrich_sample "$@"' _ {}
printf "null_model\t%s\nn_shuffles\t%s\noverlap_rule\t%s\noverlap_fraction\t%s\nsame_chrom\t%s\ninclude_bed\t%s\nexclude_bed\t%s\n" \
    "$NULL_MODEL" "$N_SHUFFLES" "$OVERLAP_RULE" "$OVERLAP_FRACTION" "$SAME_CHROM" "$INCLUDE_BED" "$EXCLUDE_BED" \
    > "$WORK/repeats/null_settings.tsv"
