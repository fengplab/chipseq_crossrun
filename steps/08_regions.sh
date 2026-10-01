#!/usr/bin/env bash
# Step 08: cross-sample region recurrence, overlap and hotspots.
#   work/regions/merged_regions.tsv   chrom start end region_id n_samples samples(comma)
#   work/regions/region_labels.tsv    region_id feature_set label      (any overlap)
#   work/regions/region_nearest_gene.tsv region_id gene distance
#   work/regions/region_input_overlap.tsv region_id n_input_peak_sets
#   work/regions/jaccard.tsv          sample_a sample_b jaccard
#   work/regions/bins.tsv             chrom start end sample_id n_peaks   (non-zero bins only)
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
mkdir -p "$WORK/regions"
D="$WORK/regions"
G=$(genome_file)

mapfile -t S < <(cut -f1 <(sample_col peaks_file))
for s in "${S[@]}"; do
    awk -v s="$s" 'BEGIN{OFS="\t"} {print $1,$2,$3,s}' "$WORK/peaks/$s.bed5"
done | sort -k1,1 -k2,2n > "$D/all_peaks.bed"

bedtools merge -i "$D/all_peaks.bed" -c 4,4 -o count_distinct,distinct \
    | awk 'BEGIN{OFS="\t"; print "chrom","start","end","region_id","n_samples","samples"}
           {print $1,$2,$3,"R" NR,$4,$5}' > "$D/merged_regions.tsv"
awk 'BEGIN{OFS="\t"} NR>1 {print $1,$2,$3,$4}' "$D/merged_regions.tsv" > "$D/regions.bed"

printf "region_id\tfeature_set\tlabel\n" > "$D/region_labels.tsv"
for fs in repeatmasker censat; do
    cb="$WORK/annot/$fs.class.bed"; [[ -s "$cb" ]] || continue
    bedtools intersect -a "$D/regions.bed" -b "$cb" -wa -wb -sorted -g "$G" \
        | awk -v fs="$fs" 'BEGIN{OFS="\t"} { k=$4 SUBSEP $8; if (!(k in s)) { s[k]=1; print $4, fs, $8 } }' \
        >> "$D/region_labels.tsv"
done

if [[ -s "$WORK/annot/genes.bed" ]]; then
    printf "region_id\tgene\tdistance\n" > "$D/region_nearest_gene.tsv"
    bedtools closest -a "$D/regions.bed" -b "$WORK/annot/genes.bed" -d -t first \
        | awk 'BEGIN{OFS="\t"} {print $4,$8,$9}' >> "$D/region_nearest_gene.tsv"
fi

# Peaks called on input/IgG controls: recurrent regions that also appear there are likely artifacts
printf "region_id\tn_input_peak_sets\n" > "$D/region_input_overlap.tsv"
if [[ -n "$INPUT_PEAKS_DIR" ]]; then
    shopt -s nullglob
    for f in "$INPUT_PEAKS_DIR"/*Peak "$INPUT_PEAKS_DIR"/*.bed; do
        rename_chroms "input_peaks:$(basename "$f")" < "$f" \
            | awk 'BEGIN{OFS="\t"} $2 ~ /^[0-9]+$/ {print $1,$2,$3}' | sort -k1,1 -k2,2n \
            | bedtools intersect -a "$D/regions.bed" -b stdin -u | cut -f4
    done | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $2,$1}' >> "$D/region_input_overlap.tsv"
fi

# Pairwise Jaccard on per-sample merged peaks
for s in "${S[@]}"; do cut -f1-3 "$WORK/peaks/$s.bed5" | bedtools merge > "$D/$s.merged.bed"; done
printf "sample_a\tsample_b\tjaccard\n" > "$D/jaccard.tsv"
for ((i=0; i<${#S[@]}; i++)); do
    for ((j=i; j<${#S[@]}; j++)); do
        jac=$(bedtools jaccard -a "$D/${S[i]}.merged.bed" -b "$D/${S[j]}.merged.bed" | awk 'NR==2{print $3}')
        printf "%s\t%s\t%s\n" "${S[i]}" "${S[j]}" "${jac:-0}" >> "$D/jaccard.tsv"
    done
done

# Hotspot bins
bedtools makewindows -g "$G" -w "$BIN_SIZE" > "$D/bins.bed"
printf "chrom\tstart\tend\tsample_id\tn_peaks\n" > "$D/bins.tsv"
for s in "${S[@]}"; do
    bedtools intersect -a "$D/bins.bed" -b "$D/$s.merged.bed" -c -sorted -g "$G" \
        | awk -v s="$s" 'BEGIN{OFS="\t"} $4>0 {print $1,$2,$3,s,$4}' >> "$D/bins.tsv"
done
log "regions: $(($(wc -l < "$D/merged_regions.tsv") - 1)) merged regions across ${#S[@]} samples"
