#!/usr/bin/env bash
# Step 07 (optional, repeats.reads_in_repeats: true): read-level signal in each repeat / CenSat label
# for every IP BAM and its control BAM, so broad satellite enrichment that never forms discrete peaks
# is still measured.
#   work/reads/<featureset>.<pe|se>.counts(.summary)   featureCounts output (label = meta-feature)
#   work/reads/reads_in_labels.tsv                     bam, feature_set, label, reads, total_reads
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
genome_file >/dev/null
[[ "$READS_IN_REPEATS" == "1" ]] || { log "reads_in_repeats disabled; skipping"; exit 0; }
command -v featureCounts >/dev/null || { log "featureCounts (subread) not found"; exit 1; }
mkdir -p "$WORK/reads"
R="$WORK/reads"

{ sample_col bam_file; sample_col control_bam; } | cut -f2 | sort -u > "$R/bams.txt"
: > "$R/pe.txt"; : > "$R/se.txt"
while read -r bam; do
    if [[ $(samtools view -f 1 "$bam" 2>/dev/null | head -n 1000 | wc -l) -gt 0 ]]; then echo "$bam" >> "$R/pe.txt"; else echo "$bam" >> "$R/se.txt"; fi
done < "$R/bams.txt"

printf "bam\tfeature_set\tlabel\treads\ttotal_reads\n" > "$R/reads_in_labels.tsv"
for fs in repeatmasker censat; do
    cb="$WORK/annot/$fs.class.bed"; [[ -s "$cb" ]] || continue
    for mode in pe se; do
        [[ -s "$R/$mode.txt" ]] || continue
        # SAF chromosome names must match the BAM header (which may use e.g. NC_060925.1)
        samtools view -H "$(head -1 "$R/$mode.txt")" | awk -F'\t' '$1=="@SQ"{sub(/^SN:/,"",$2); print $2}' > "$R/$mode.bam_chroms"
        awk 'BEGIN{FS=OFS="\t"} FILENAME==ARGV[1] { a[$1]=$2; next }
             FILENAME==ARGV[2] { c = ($1 in a) ? a[$1] : $1; bam[c]=$1; next }
             FNR==1 { print "GeneID","Chr","Start","End","Strand" }
             ($1 in bam) { print $4, bam[$1], $2+1, $3, "+" }' \
            "$WORK/annot/chrom_alias.tsv" "$R/$mode.bam_chroms" "$cb" > "$WORK/annot/$fs.$mode.saf"
        [[ $(wc -l < "$WORK/annot/$fs.$mode.saf") -gt 1 ]] || { log "ERROR: no $fs features match BAM chromosome names"; exit 1; }
        extra=""; [[ $mode == pe ]] && extra="-p --countReadPairs"
        featureCounts -a "$WORK/annot/$fs.$mode.saf" -F SAF -O -T "$THREADS" $extra \
            -o "$R/$fs.$mode.counts" $(cat "$R/$mode.txt") 2> "$R/$fs.$mode.log"
        python3 - "$R/$fs.$mode.counts" "$fs" >> "$R/reads_in_labels.tsv" <<'PY'
import sys
cnt, fs = sys.argv[1], sys.argv[2]
tot = {}
with open(cnt + ".summary") as fh:
    bams = fh.readline().rstrip("\n").split("\t")[1:]
    for line in fh:
        v = line.rstrip("\n").split("\t")[1:]
        for b, x in zip(bams, v):
            tot[b] = tot.get(b, 0) + int(x)
with open(cnt) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] == "Geneid":
            cols = f[6:]
            continue
        for b, x in zip(cols, f[6:]):
            print(f"{b}\t{fs}\t{f[0]}\t{x}\t{tot[b]}")
PY
    done
done
log "reads in repeat labels written"
