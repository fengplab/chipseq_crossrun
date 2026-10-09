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
