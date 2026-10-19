#!/usr/bin/env bash
# Step 19: extra figures, written to <outdir>/extra/ and NOT part of the report (extra_plots in config.yaml).
#   extra/coverage/<sample>/<contig>_<norm>.png   coverage along every host chromosome, chrM and viral genome, for each
#                                                 normalisation in extra_plots.coverage_norms (binned BAM counts are cached
#                                                 in work/coverage/, so a rerun only redraws). Viral genomes are drawn from
#                                                 the pipeline's duplicate-marked viral BAM when it exists; chrM and viral
#                                                 genomes with no reads in the BAM used (the pipeline removes them before
#                                                 the filtered BAM) are skipped with a note.
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
    cut -f1 "$WORK/annot/genome.sizes" > "$WORK/annot/genome.sizes.list"     # host nuclear chromosomes: always drawn
    export ROOT X
    cov_sample() {
        local s="$1" bam c n_ok=0 n_skip=0 args=() norms
        bam=$(sample_get "$s" bam_file); [[ -s "$bam" ]] || { log "WARNING: $s has no BAM; no coverage plots"; return 0; }
        local vbam skipped=() src
        vbam=$(sample_get "$s" viral_bam)
        mkdir -p "$X/coverage/$s"
        # reads per contig in the filtered BAM, and in the viral BAM
        samtools idxstats "$bam" | awk 'BEGIN{OFS="\t"} {print $1, $3}' | from_any_names > "$X/coverage/$s/.idx_main"
        : > "$X/coverage/$s/.idx_viral"
        [[ -n "$vbam" && -s "$vbam" ]] && samtools idxstats "$vbam" | awk 'BEGIN{OFS="\t"} {print $1, $3}' | from_any_names > "$X/coverage/$s/.idx_viral"
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
            # viral genomes: same MAPQ cutoff as the read counts (contigs.viral_mapq), from the viral BAM when there is one;
            # host contigs unfiltered. Contigs without reads in the BAM used are skipped (e.g. chrM removed by the pipeline).
            local mq=0 vflag=() nr
            src="$bam"
            if grep -qxF "$c" "$X/coverage/viral.txt"; then
                mq="${VIRAL_MAPQ:-10}"
                if [[ -s "$X/coverage/$s/.idx_viral" ]] && awk -F'\t' -v c="$c" '$1 == c && $2 > 0 {f=1} END {exit !f}' "$X/coverage/$s/.idx_viral"; then
                    vflag=(--viral-bam); src="$vbam"
                fi
            fi
            if [[ -z "${vflag[*]:-}" ]] && ! grep -qxF "$c" "$WORK/annot/genome.sizes.list"; then
                nr=$(awk -F'\t' -v c="$c" '$1 == c {n += $2} END {print n + 0}' "$X/coverage/$s/.idx_main")
                if (( nr == 0 )); then skipped+=("$c"); continue; fi
            fi
            if Rscript "$ROOT/R/plot_chrom_coverage.R" --chrom "$c" --crossrun "$OUTDIR" --samples "$s" --norm "$EXTRA_NORMS" --mapq "$mq" "${vflag[@]}" \
                   --cache "$WORK/coverage" --threads 1 --title "$s: $c" --out "$X/coverage/$s/{chrom}_{norm}.png" "${args[@]}" \
                   > "$X/coverage/$s/.log_$c" 2>&1; then
                n_ok=$((n_ok + 1)); rm -f "$X/coverage/$s/.log_$c"
            else
                log "WARNING: coverage plot failed for $s $c (see extra/coverage/$s/.log_$c)"
            fi
        done < "$X/coverage/contigs.txt"
        rm -f "$X/coverage/$s/.idx_main" "$X/coverage/$s/.idx_viral"
        log "coverage plots: $s ($n_ok contigs drawn, $n_skip up to date${skipped[*]:+; no reads, not drawn: ${skipped[*]}})"
    }
    export -f cov_sample
    sample_col bam_file | cut -f1 | xargs -P "$THREADS" -I{} bash -c 'cov_sample "$@"' _ {}
fi

if [[ "${EXTRA_PIES:-1}" == "1" ]]; then
    Rscript "$ROOT/R/plot_pies.R" "$CFGY"
fi
[[ "${EXTRA_INSERT_SIZES:-1}" == "1" ]] && log "insert-size figures (formerly Figure QC3) are written by step 15 to extra/insert_sizes/"
log "extra figures in $X"
