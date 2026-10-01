#!/usr/bin/env bash
# Optional: call peaks on each control (input/IgG) BAM with MACS3, without a control, so recurrent
# ChIP regions that are also "peaks" in inputs can be flagged as likely artifacts.
# Output: <outdir>/work/input_peaks/<control>_peaks.narrowPeak  -> set regions.input_peaks_dir to this dir.
# Usage: optional_call_input_peaks.sh <work>/config.env <macs_gsize, e.g. hs or 3.1e9>
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
GSIZE="${2:?genome size for MACS3 -g}"
mkdir -p "$WORK/input_peaks"
awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["control_bam"]!="" {print $h["control_id"] "\t" $h["control_bam"]}' "$SAMPLES" \
    | sort -u | while IFS=$'\t' read -r cid bam; do
    fmt=BAM; [[ $(samtools view -f 1 "$bam" | head -n 1000 | wc -l) -gt 0 ]] && fmt=BAMPE
    macs3 callpeak -t "$bam" -f "$fmt" -g "$GSIZE" -n "$cid" --outdir "$WORK/input_peaks" --keep-dup all -q 0.05 \
        2> "$WORK/input_peaks/$cid.log"
    log "input peaks: $cid"
done
