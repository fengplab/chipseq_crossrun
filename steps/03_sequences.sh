#!/usr/bin/env bash
# Step 03: sequences for motif scanning.
#   work/seqs/<sample>.scan.bed    windows scanned (whole peak, or summit +/- flank)
#   work/seqs/<sample>.fa          their sequences (header: <sample>|<peak>::chrom:start-end)
#   work/seqs/<sample>.bg.bed/.fa  one null relocation of the same windows (genomic background)
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
mkdir -p "$WORK/seqs"

one_sample() {
    set -euo pipefail
    local s="$1" b="$WORK/peaks/$1.bed5" o="$WORK/seqs/$1"
    # reuse finished sequences only if they were made with the same settings and sequence-name format
    local stamp="v=$SCAN_VERSION;mode=$SEQ_MODE;flank=$SUMMIT_FLANK;null=$NULL_MODEL"
    [[ -s "$o.fa" && -s "$o.bg.fa" && "$(cat "$o.stamp" 2>/dev/null)" == "$stamp" && ! "$b" -nt "$o.fa" ]] && return 0
    if [[ "$SEQ_MODE" == "summit" ]]; then
        awk -v f="$SUMMIT_FLANK" 'BEGIN{FS=OFS="\t"} { st=$5-f; if (st<0) st=0; print $1, st, $5+f, $4 }' "$b"
    else
        cut -f1-4 "$b"
    fi | sort -k1,1 -k2,2n > "$o.scan.bed"
    getfasta_canonical "$o.scan.bed" "$o.fa"
    shuffle_once "$o.scan.bed" 999 > "$o.bg.bed"
    getfasta_canonical "$o.bg.bed" "$o.bg.fa"
    printf "%s\n" "$stamp" > "$o.stamp"
    log "sequences: $s"
}
export -f one_sample
genome_file >/dev/null
cut -f1 <(sample_col peaks_file) | xargs -P "$THREADS" -I{} bash -c 'one_sample "$@"' _ {}
