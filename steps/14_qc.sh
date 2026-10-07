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
    if [[ ! -r "$bam" ]]; then
        log "ERROR: cannot read BAM for $id: $bam"
        log "  (if the pipeline results moved: update results_dir in the manifest and rerun step 0)"
        return 1
    fi
    # reuse paired-end results; a cached single-end verdict is re-checked because it may come from a failed read
    [[ -s "$out" && "$(cat "$Q/insert/$id.status" 2>/dev/null)" == paired ]] && return 0
    if [[ $(samtools view -f 1 "$bam" | head -n 1000 | wc -l) -eq 0 ]]; then
        samtools quickcheck "$bam" || { log "ERROR: $bam is not a valid BAM"; return 1; }
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
            | isect_sorted -b "$CEN" -wa -wb -sorted -g "$G" \
            | awk 'BEGIN{OFS="\t"} { n[$8 "\t" $4]++ } END { for (x in n) print x, n[x] }'
        frac=$(awk -v f="$QC_INSERT_GENOME_FRACTION" 'BEGIN{ if (f >= 1) print ""; else printf "13%s", substr(sprintf("%.9f", f), 2) }')
        samtools view ${frac:+-s $frac} -f 66 -F 3852 "$bam" | frag \
            | isect_sorted -v -b "$CEN" -sorted -g "$G" \
            | awk 'BEGIN{OFS="\t"} { n[$4]++ } END { for (x in n) print "outside_censat", x, n[x] }'
    } > "$out.tmp.$$"
    mv "$out.tmp.$$" "$out"; rm -f "$Q/insert/$id.chroms" "$Q/insert/$id.censat.bed"
    log "insert sizes: $id"
}
export -f insert_one
export Q G CEN
# One entry per distinct BAM path. Library IDs must be unique because each names its output file: the same
# control is often present as a separate copy in several runs' results folders under an identical file name,
# so a name shared by different paths is prefixed with the run ID (e.g. run03__THP1input_REP1).
python3 - "$SAMPLES" > "$Q/insert/libraries.tsv" <<'PY'
import csv, os, re, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
paths = {}
for r in rows:
    for col in ("bam_file", "control_bam"):
        p = (r.get(col) or "").strip()
        if p and p not in paths:
            paths[p] = r.get("run_id", "") or "run"
base = lambda p: re.sub(r"(\.mLb\.clN\.sorted)?\.bam$", "", os.path.basename(p))
count = {}
for p in paths:
    count[base(p)] = count.get(base(p), 0) + 1
ids = {}
for p, run in paths.items():
    i = base(p) if count[base(p)] == 1 else f"{re.sub(r'[^A-Za-z0-9._-]+', '_', run)}__{base(p)}"
    j, n = i, 2
    while j in ids.values():          # same run and same name but different paths: number them
        j, n = f"{i}_{n}", n + 1
    ids[p] = j
for p, i in ids.items():
    print(f"{i}\t{p}")
PY
awk -F'\t' '{print $1; print $2}' "$Q/insert/libraries.tsv" | xargs -n 2 -P "$THREADS" bash -c 'insert_one "$@"' _

# ---------------------------------------------------------------- 3. peak summits within CenSat arrays
# Each peak summit inside a CenSat feature (array) with that array's coordinates, for distance-to-edge QC.
#   work/qc/array_position.tsv   peak_id  summit  array_start  array_end  class
if [[ ! -s "$Q/array_position.tsv" || "$Q/array_position.tsv" -ot "$WORK/peaks/peak_stats.tsv" ]]; then
    {
        printf "peak_id\tsummit\tarray_start\tarray_end\tclass\n"
        cat "$WORK"/peaks/*.bed5 | awk 'BEGIN{OFS="\t"} {print $1, $5, $5 + 1, $4}' | sort -k1,1 -k2,2n \
            | isect_sorted -b "$CEN" -wa -wb -sorted -g "$G" \
            | awk 'BEGIN{OFS="\t"} {print $4, $2, $6, $7, $8}'
    } > "$Q/array_position.tsv.tmp" && mv "$Q/array_position.tsv.tmp" "$Q/array_position.tsv"
    log "array positions: $(($(wc -l < "$Q/array_position.tsv") - 1)) peak summits inside CenSat arrays"
fi

# ---------------------------------------------------------------- 4. strand balance and read stacking per peak
# For every peak of every IP sample: read 5' start positions (read 1 only for paired-end, so each fragment
# counts once) -> number of reads, forward-strand reads, distinct start positions (position + strand) and
# reads at the most-used start position.
#   work/qc/peak_shape/<sample>.tsv   peak_id  n_reads  n_forward  n_distinct_starts  n_top_start
mkdir -p "$Q/peak_shape"
peak_shape_one() {
    set -euo pipefail
    local s="$1" bam="$2" out="$Q/peak_shape/$1.tsv" pk="$WORK/peaks/$1.bed5" t
    [[ -r "$bam" ]] || { log "ERROR: cannot read BAM for $s: $bam"; return 1; }
    if [[ -s "$out" && "$out" -nt "$pk" && "$out" -nt "$bam" ]]; then return 0; fi
    t=$(mktemp -d)
    samtools view -H "$bam" | awk -F'\t' '$1=="@SQ"{sub(/^SN:/,"",$2); print $2}' > "$t/chroms"
    # peaks in this BAM's chromosome names (for indexed retrieval) and in canonical names (for intersect)
    cut -f1-4 "$pk" > "$t/peaks.bed"
    awk 'BEGIN{FS=OFS="\t"} FILENAME==ARGV[1] { a[$1]=$2; next } FILENAME==ARGV[2] { c=($1 in a)?a[$1]:$1; b[c]=$1; next }
         ($1 in b) { print b[$1], $2, $3 }' "$WORK/annot/chrom_alias.tsv" "$t/chroms" "$t/peaks.bed" > "$t/peaks_bam.bed"
    local flt="-F 3844"
    [[ $(samtools view -f 1 "$bam" | head -n 1000 | wc -l) -gt 0 ]] && flt="-f 64 -F 3852"
    {
        printf "peak_id\tn_reads\tn_forward\tn_distinct_starts\tn_top_start\n"
        samtools view -M -L "$t/peaks_bam.bed" $flt "$bam" \
            | awk 'BEGIN{OFS="\t"} {
                  rl = 0; c = $6
                  while (match(c, /[0-9]+[MDN=X]/)) { rl += substr(c, RSTART, RLENGTH - 1); c = substr(c, RSTART + RLENGTH) }
                  if (int($2 / 16) % 2 == 1) { p = $4 - 1 + rl - 1; st = "-" } else { p = $4 - 1; st = "+" }
                  print $3, p, p + 1, st }' \
            | from_any_names | awk 'BEGIN{FS=OFS="\t"} NR==FNR{ok[$1]=1; next} ($1 in ok)' "$G" - | sort -k1,1 -k2,2n \
            | isect_sorted -b "$t/peaks.bed" -wa -wb -sorted -g "$G" \
            | awk 'BEGIN{OFS="\t"} { id = $8; n[id]++; if ($4 == "+") f[id]++
                                      k = id SUBSEP $2 SUBSEP $4; if (!(k in c)) d[id]++; c[k]++
                                      if (c[k] > top[id]) top[id] = c[k] }
                   END { for (id in n) print id, n[id], f[id] + 0, d[id], top[id] }'
    } > "$out.tmp.$$"
    mv "$out.tmp.$$" "$out"; rm -rf "$t"
    log "peak strand/stacking: $s"
}
export -f peak_shape_one
sample_col bam_file | awk -F'\t' '{print $1; print $2}' | xargs -n 2 -P "$THREADS" bash -c 'peak_shape_one "$@"' _
