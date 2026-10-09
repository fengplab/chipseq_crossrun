#!/usr/bin/env bash
# Step 01: build the chromosome alias table, then normalise annotations and peak files.
# Every output uses canonical UCSC-style chromosome names (chr1..chr22, chrX, chrY, chrM).
#   work/annot/chrom_alias.tsv, fasta_names.tsv, genome.sizes, chrom_naming_report.tsv
#   work/annot/{repeatmasker,censat}.class.bed   chrom start end label (sorted)
#   work/annot/genes.bed, include.bed, exclude.bed (if configured)
#   work/annot/rename_log.tsv                    per input: records, renamed, unrecognised names
#   work/annot/chrom_overlap.tsv                 chromosomes shared between peaks and each annotation
#   work/peaks/<sample>.bed5                     chrom start end sample|peak summit
#   work/peaks/peak_stats.tsv                    sample, n_peaks, median_width, bp_covered
#   work/peaks/<sample>.info.tsv                 per peak: width, signal value (MACS fold enrichment), -log10 p, -log10 q
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
BIN="$(cd "$(dirname "$0")/../bin" && pwd)"
mkdir -p "$WORK/annot" "$WORK/peaks"
printf "input\trecords\trenamed\tunrecognised\texamples\n" > "$WORK/annot/rename_log.tsv"
G=$(genome_file)
A="$WORK/annot/chrom_alias.tsv"

to_bed() {  # bigBed -> BED on the fly, anything else passed through
    local f="$1"
    if [[ "$f" == *.bb || "$f" == *.bigBed || "$f" == *.bigbed ]]; then
        command -v bigBedToBed >/dev/null || { echo "bigBedToBed needed for $f" >&2; exit 1; }
        bigBedToBed "$f" /dev/stdout
    else
        zcat -f "$f"
    fi
}

for kind in repeatmasker censat; do
    src=$([[ $kind == repeatmasker ]] && echo "$REPEATMASKER" || echo "$CENSAT")
    out="$WORK/annot/$kind.class.bed"
    if [[ -z "$src" ]]; then log "no $kind annotation configured; skipping"; rm -f "$out"; continue; fi
    if [[ ! -s "$out" ]]; then
        log "preparing $kind labels"
        to_bed "$src" | python3 "$BIN/features_to_classbed.py" --kind "$kind" --genome "$G" --alias "$A" - \
            | sort -k1,1 -k2,2n > "$out"
    fi
done

if [[ -n "$GENES_BED" ]]; then
    zcat -f "$GENES_BED" | rename_chroms genes \
        | awk 'BEGIN{FS=OFS="\t"} NR==FNR{ok[$1]=1;next} ($1 in ok) && $2 ~ /^[0-9]+$/ {print $1,$2,$3,($4==""?"gene":$4)}' "$G" - \
        | sort -k1,1 -k2,2n > "$WORK/annot/genes.bed"
fi
for which in include exclude; do
    src=$([[ $which == include ]] && echo "$INCLUDE_BED" || echo "$EXCLUDE_BED")
    [[ -n "$src" ]] || continue
    zcat -f "$src" | rename_chroms "${which}_bed" \
        | awk 'BEGIN{FS=OFS="\t"} NR==FNR{ok[$1]=1;next} ($1 in ok) && $2 ~ /^[0-9]+$/ {print $1,$2,$3}' "$G" - \
        | sort -k1,1 -k2,2n | bedtools merge > "$WORK/annot/$which.bed"
    [[ -s "$WORK/annot/$which.bed" ]] || { log "ERROR: ${which}_bed has no regions on genome chromosomes after renaming"; exit 1; }
done

printf "sample_id\tn_peaks\tmedian_width\tbp_covered\n" > "$WORK/peaks/peak_stats.tsv"
while IFS=$'\t' read -r s f; do
    peaks_to_bed5 "$s" "$f" > "$WORK/peaks/$s.bed5"
    { printf "chrom\tstart\tend\twidth\tsignal_value\tneglog10_p\tneglog10_q\n"; peaks_info "$s" "$f"; } > "$WORK/peaks/$s.info.tsv"
    b="$WORK/peaks/$s.bed5"
    [[ -s "$b" ]] || log "WARNING: $s has no peaks on genome chromosomes after renaming (see rename_log.tsv)"
    med=$(awk '{print $3-$2}' "$b" | sort -n | awk '{a[NR]=$1} END{ if (NR==0) print 0; else if (NR%2) print a[(NR+1)/2]; else print (a[NR/2]+a[NR/2+1])/2 }')
    printf "%s\t%s\t%s\t%s\n" "$s" "$(wc -l < "$b" | tr -d ' ')" "$med" "$(awk '{t+=$3-$2} END{print t+0}' "$b")" \
        >> "$WORK/peaks/peak_stats.tsv"
done < <(sample_col peaks_file)

# Consistency check: peaks and each annotation must share chromosomes, otherwise every overlap is 0
cut -f1 "$WORK"/peaks/*.bed5 2>/dev/null | sort -u > "$WORK/annot/peak_chroms.txt"
printf "annotation\tchroms_in_annotation\tchroms_with_peaks\tshared\n" > "$WORK/annot/chrom_overlap.tsv"
for f in "$WORK"/annot/{repeatmasker,censat}.class.bed "$WORK/annot/genes.bed"; do
    [[ -s "$f" ]] || continue
    cut -f1 "$f" | sort -u > "$f.chroms"
    sh=$(comm -12 "$f.chroms" "$WORK/annot/peak_chroms.txt" | wc -l | tr -d ' ')
    printf "%s\t%s\t%s\t%s\n" "$(basename "$f")" "$(wc -l < "$f.chroms" | tr -d ' ')" \
        "$(wc -l < "$WORK/annot/peak_chroms.txt" | tr -d ' ')" "$sh" >> "$WORK/annot/chrom_overlap.tsv"
    [[ "$sh" -gt 0 ]] || log "ERROR: $(basename "$f") shares no chromosome names with the peaks; add a chrom_alias file"
    rm -f "$f.chroms"
done
awk -F'\t' 'NR>1 && $4>0 { log_ = log_ sprintf("  %s: %d records, %d renamed, %d unrecognised %s\n", $1, $2, $3, $4, ($5 != "" ? "(e.g. " $5 ")" : "")) }
            END { if (log_ != "") printf "[step01] inputs with unrecognised chromosome names:\n%s", log_ > "/dev/stderr" }' "$WORK/annot/rename_log.tsv"
log "prepared $(($(wc -l < "$WORK/peaks/peak_stats.tsv") - 1)) peak sets; chromosome names: $(head -3 "$G" | cut -f1 | paste -sd, -), ..."
