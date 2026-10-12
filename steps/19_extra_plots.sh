#!/usr/bin/env bash
# Step 19: extra figures, written to <outdir>/extra/ and NOT part of the report (extra_plots in config.yaml).
#   extra/coverage/<sample>/<contig>_<norm>.png   coverage along every host chromosome, chrM and viral genome, for each
#                                                 normalisation in extra_plots.coverage_norms (binned BAM counts are cached
#                                                 in work/coverage/, so a rerun only redraws)
#   extra/pies/                                   reads and peaks per contig, per sample and pooled (R/plot_pies.R)
#   extra/insert_sizes/                           per-sample insert-size distributions (made in step 15)
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[[ "${EXTRA_ENABLED:-1}" == "1" ]] || { log "extra plots disabled (extra_plots.enabled: false); skipping"; exit 0; }
X="$OUTDIR/extra"; mkdir -p "$X"
genome_file >/dev/null
CFGY="$CONFIG_YAML"

if [[ "${EXTRA_COVERAGE:-1}" == "1" ]]; then
    command -v samtools >/dev/null || { log "ERROR: samtools not found (needed for the coverage plots)"; exit 1; }
    mkdir -p "$X/coverage" "$WORK/coverage"
    cut -f1 "$WORK/annot/genome_all.sizes" > "$X/coverage/contigs.txt"
    awk -F'\t' 'NR > 1 && $3 == "viral" {print $1}' "$WORK/annot/contigs.tsv" > "$X/coverage/viral.txt"
    export ROOT X
    cov_sample() {
        local s="$1" bam c n_ok=0 n_skip=0 args=() norms
        bam=$(sample_get "$s" bam_file); [[ -s "$bam" ]] || { log "WARNING: $s has no BAM; no coverage plots"; return 0; }
        mkdir -p "$X/coverage/$s"
        [[ -n "$EXTRA_BANDS" ]] && args+=(--bands "$EXTRA_BANDS")
        [[ -n "$EXTRA_CAP" ]] && args+=(--cap "$EXTRA_CAP")
        while read -r c; do
            # up to date when every requested plot exists and is newer than the BAM
            local done=1 nrm
            for nrm in ${EXTRA_NORMS//,/ }; do
                [[ "$nrm" == input && -z "$(sample_get "$s" control_bam)" ]] && continue
                [[ -s "$X/coverage/$s/${c}_$nrm.png" && "$X/coverage/$s/${c}_$nrm.png" -nt "$bam" ]] || done=0
            done
            if (( done )); then n_skip=$((n_skip + 1)); continue; fi
            # viral genomes: same MAPQ cutoff as the read counts (contigs.viral_mapq); host contigs unfiltered
            local mq=0; grep -qxF "$c" "$X/coverage/viral.txt" && mq="${VIRAL_MAPQ:-10}"
            if Rscript "$ROOT/R/plot_chrom_coverage.R" --chrom "$c" --crossrun "$OUTDIR" --samples "$s" --norm "$EXTRA_NORMS" --mapq "$mq" \
                   --cache "$WORK/coverage" --threads 1 --title "$s: $c" --out "$X/coverage/$s/{chrom}_{norm}.png" "${args[@]}" \
                   > "$X/coverage/$s/.log_$c" 2>&1; then
                n_ok=$((n_ok + 1)); rm -f "$X/coverage/$s/.log_$c"
            else
                log "WARNING: coverage plot failed for $s $c (see extra/coverage/$s/.log_$c)"
            fi
        done < "$X/coverage/contigs.txt"
        log "coverage plots: $s ($n_ok contigs drawn, $n_skip up to date)"
    }
    export -f cov_sample
    sample_col bam_file | cut -f1 | xargs -P "$THREADS" -I{} bash -c 'cov_sample "$@"' _ {}
fi

if [[ "${EXTRA_PIES:-1}" == "1" ]]; then
    Rscript "$ROOT/R/plot_pies.R" "$CFGY"
fi
[[ "${EXTRA_INSERT_SIZES:-1}" == "1" ]] && log "insert-size figures (formerly Figure QC3) are written by step 15 to extra/insert_sizes/"
log "extra figures in $X"
