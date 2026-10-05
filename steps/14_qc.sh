#!/usr/bin/env bash
# Step 14: quality-control computations (report section "Quality control").
#  1. Mappability: single-copy k-mer mappability with k = qc.mappability_k (GenMap, both strands, 0
#     mismatches; or a supplied bedGraph), converted to per-base mappability and summarised for every
#     CenSat class and every peak (bin/qc_mappability.py). Reuses work/peaks and work/annot from step 01.
#  2. Insert sizes: for every IP and control BAM (paired-end only), the fragment-length distribution of
#     fragments whose midpoint lies in each CenSat class, plus a random sample of fragments outside CenSat.
#   work/qc/censat_mappability_hist.tsv, peak_mappability.tsv.gz, mappability_k<k>.bedgraph.gz
#   work/qc/insert/<library>.tsv   class  insert_size  n      (<library>.status = paired | single_end)
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[[ "$QC_ENABLED" == "1" ]] || { log "qc.enabled is false; skipping"; exit 0; }
Q="$WORK/qc"; mkdir -p "$Q/insert"
G=$(genome_file)
CEN="$WORK/annot/censat.class.bed"

# GenMap's -O is an output FOLDER when several FASTA files were indexed, but a file-name PREFIX when only
# one was (the usual case): it then writes <prefix>.bedgraph (or .bg in some versions) next to it.
find_genmap_output() {
    local p="$Q/genmap_k$QC_K" f
    for f in "$p.bedgraph" "$p.bg" "$p".*.bedgraph "$p".*.bg "$p"/*.bedgraph "$p"/*.bg; do
        [[ -s "$f" ]] && { echo "$f"; return 0; }
    done
    return 0
}

# ---------------------------------------------------------------- 1. mappability
if [[ ! -s "$Q/censat_mappability_hist.tsv" || ! -s "$Q/peak_mappability.tsv.gz" || "$Q/peak_mappability.tsv.gz" -ot "$WORK/peaks/peak_stats.tsv" ]]; then
    BG="$QC_MAPPABILITY_BEDGRAPH"
    if [[ -z "$BG" ]]; then
        BG=$(find_genmap_output)
        if [[ -z "$BG" ]]; then
            command -v genmap >/dev/null || { log "ERROR: genmap not found (conda install -c bioconda genmap), or set qc.mappability_bedgraph"; exit 1; }
            IDX="${QC_GENMAP_INDEX:-$Q/genmap_index}"
            if [[ ! -d "$IDX" ]]; then
                log "GenMap: indexing the genome (one-off; reuse it later via qc.genmap_index)"
                genmap index -F "$GENOME_FASTA" -I "$IDX" > "$Q/genmap_index.log" 2>&1
            fi
            log "GenMap: single-copy $QC_K-mers (0 mismatches, both strands)"
            genmap map -K "$QC_K" -E 0 -I "$IDX" -O "$Q/genmap_k$QC_K" -bg -T "$THREADS" > "$Q/genmap_map.log" 2>&1
            BG=$(find_genmap_output)
            if [[ -z "$BG" ]]; then
                log "ERROR: GenMap finished but no bedGraph was found. Contents of $Q:"
                ls -la "$Q" >&2
                log "Last lines of $Q/genmap_map.log:"
                tr '\r' '\n' < "$Q/genmap_map.log" | sed 's/\x1b\[K//g' | grep -v '^Progress' | tail -5 >&2
                exit 1
            fi
        else
            log "GenMap: reusing $BG"
        fi
    fi
    log "per-base mappability from $BG"
    python3 "$ROOT/bin/qc_mappability.py" --bedgraph "$BG" --alias "$WORK/annot/chrom_alias.tsv" --sizes "$G" \
        --k "$QC_K" --censat "$CEN" --peaks-dir "$WORK/peaks" --outdir "$Q" $([[ "$QC_WRITE_TRACK" == "1" ]] && echo --track)
else
    log "mappability: reusing $Q/peak_mappability.tsv.gz"
fi

# ---------------------------------------------------------------- 2. insert sizes per CenSat class
[[ -s "$CEN" ]] || { log "no CenSat annotation: insert-size QC skipped"; exit 0; }
insert_one() {
    set -euo pipefail
    local id="$1" bam="$2" out="$Q/insert/$1.tsv" frac
    [[ -s "$out" ]] && return 0
    if [[ $(samtools view -f 1 "$bam" | head -n 1000 | wc -l) -eq 0 ]]; then
        echo single_end > "$Q/insert/$id.status"; printf "class\tinsert_size\tn\n" > "$out"; return 0
    fi
    echo paired > "$Q/insert/$id.status"
    # CenSat regions in this BAM's own chromosome names, for indexed retrieval
    samtools view -H "$bam" | awk -F'\t' '$1=="@SQ"{sub(/^SN:/,"",$2); print $2}' > "$Q/insert/$id.chroms"
    awk 'BEGIN{FS=OFS="\t"} FILENAME==ARGV[1] { a[$1]=$2; next } FILENAME==ARGV[2] { c=($1 in a)?a[$1]:$1; b[c]=$1; next }
         ($1 in b) { print b[$1], $2, $3 }' "$WORK/annot/chrom_alias.tsv" "$Q/insert/$id.chroms" "$CEN" > "$Q/insert/$id.censat.bed"
    frag() {   # read-1 records of proper pairs -> chrom  mid  mid+1  insert_size  (canonical names)
        awk -v mx="$QC_INSERT_MAX" 'BEGIN{OFS="\t"} $7=="=" && $9 != 0 { t = ($9 < 0) ? -$9 : $9; if (t > mx) next
            l = (($4 < $8) ? $4 : $8) - 1; m = l + int(t/2); print $3, m, m + 1, t }' | from_any_names \
            | awk 'BEGIN{FS=OFS="\t"} NR==FNR{ok[$1]=1; next} ($1 in ok)' "$G" - | sort -k1,1 -k2,2n
    }
    {
        printf "class\tinsert_size\tn\n"
        samtools view -M -L "$Q/insert/$id.censat.bed" -f 66 -F 3852 "$bam" | frag \
            | bedtools intersect -a stdin -b "$CEN" -wa -wb -sorted -g "$G" \
            | awk 'BEGIN{OFS="\t"} { n[$8 "\t" $4]++ } END { for (x in n) print x, n[x] }'
        frac=$(awk -v f="$QC_INSERT_GENOME_FRACTION" 'BEGIN{ if (f >= 1) print ""; else printf "13%s", substr(sprintf("%.9f", f), 2) }')
        samtools view ${frac:+-s $frac} -f 66 -F 3852 "$bam" | frag \
            | bedtools intersect -v -a stdin -b "$CEN" -sorted -g "$G" \
            | awk 'BEGIN{OFS="\t"} { n[$4]++ } END { for (x in n) print "outside_censat", x, n[x] }'
    } > "$out.tmp"
    mv "$out.tmp" "$out"; rm -f "$Q/insert/$id.chroms" "$Q/insert/$id.censat.bed"
    log "insert sizes: $id"
}
export -f insert_one
export Q G CEN
{ sample_col bam_file; sample_col control_bam; } | cut -f2 | sort -u \
    | awk '{ id=$0; sub(/^.*\//, "", id); sub(/\.mLb\.clN\.sorted\.bam$/, "", id); sub(/\.bam$/, "", id); print id "\t" $0 }' \
    > "$Q/insert/libraries.tsv"
awk -F'\t' '{print $1; print $2}' "$Q/insert/libraries.tsv" | xargs -n 2 -P "$THREADS" bash -c 'insert_one "$@"' _
