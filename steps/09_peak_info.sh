#!/usr/bin/env bash
# Per-peak details for the peak overview (R/05_peaks.R): work/peaks/<sample>.info.tsv
# Step 01 writes these; this fills in any that are missing or older than their peak file, so outputs made
# before the peak overview existed do not need step 01 rerun.
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
mkdir -p "$WORK/peaks"
n=0
while IFS=$'\t' read -r s f; do
    o="$WORK/peaks/$s.info.tsv"
    if [[ ! -s "$o" || "$f" -nt "$o" ]]; then
        { printf "chrom\tstart\tend\twidth\tsignal_value\tneglog10_p\tneglog10_q\n"; peaks_info "$s" "$f"; } > "$o.tmp" && mv "$o.tmp" "$o"
        n=$((n + 1))
    fi
done < <(sample_col peaks_file)
(( n == 0 )) || log "peak details written for $n samples"

# Reads per contig and MAPQ for every IP and control BAM (Figure P5, table00c, qc07, pies). Cached per BAM.
genome_file >/dev/null
if command -v samtools >/dev/null; then
    { sample_col bam_file | awk 'BEGIN{OFS="\t"} {print $1, "IP", $2}'
      sample_col control_bam | awk 'BEGIN{OFS="\t"} {print $1, "control", $2}'; } > "$WORK/contigs.list.tmp" 2>/dev/null || true
    if [[ -s "$WORK/contigs.list.tmp" ]]; then
        mkdir -p "$WORK/contigs"
        cut -f3 "$WORK/contigs.list.tmp" | sort -u | xargs -P "$THREADS" -I{} bash -c 'contig_counts "$1" > /dev/null' _ {}
        { printf "sample_id\tlibrary\tbam\tcounts_file\tpaired\n"
          while IFS=$'\t' read -r s lib b; do
              [[ -s "$b" ]] || continue
              pe=$([[ $(samtools view -f 1 "$b" 2>/dev/null | head -n 1000 | wc -l) -gt 0 ]] && echo 1 || echo 0)
              printf "%s\t%s\t%s\t%s\t%s\n" "$s" "$lib" "$b" "$(contig_counts "$b")" "$pe"
          done < "$WORK/contigs.list.tmp"; } > "$WORK/contigs/bams.tsv"
        log "reads per contig counted for $(cut -f3 "$WORK/contigs.list.tmp" | sort -u | wc -l) BAM files"
    fi
    rm -f "$WORK/contigs.list.tmp"
else
    log "WARNING: samtools not found; reads per contig (Figure P5, pies) are skipped"
fi
