#!/usr/bin/env bash
# Step 16: bin table for the binding models.
#   work/bins/w<size>/manifest.tsv, classes.tsv, arrays.tsv, <chrom>/{class,array_id,edge_dist,repeat_frac,mappability,gc,n_frac}.*
#   work/bins/fimo/<family>.tsv.gz        genome-wide FIMO sites of each family representative (cached per motif + genome)
#   work/bins/w<size>/sites.tsv.gz        family, chrom, bin, n (sites whose midpoint lies in the bin)
#   work/bins/w<size>/<library>/          per-library signal: all, uniq, h1, h2, fl per chromosome + library.json
#   work/bins/w<size>/libraries.tsv       sample_id, library (IP/control), name, dir
#   work/bins/drift_check.tsv             signal retained at the chosen size vs the tuning floor (bins.drift_check)
# Bin size: bins.size, or chosen_size from bins.tuning_record (tune/run_tuning.sh).
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[[ "${MODELS_ENABLED:-1}" == "1" ]] || { log "binding models disabled (binding_models.enabled: false); skipping"; exit 0; }
if [[ -z "${BINS_SIZE:-}" ]]; then
    log "ERROR: no bin size. Run the bin-size tuning once (tune/run_tuning.sh <config>) and set bins.tuning_record,"
    log "       or set bins.size in the config, then rerun step 0."
    exit 1
fi
W="$BINS_SIZE"; B="$WORK/bins"; D="$B/w$W"
mkdir -p "$D" "$B/fimo"
genome_file >/dev/null
SIZES="$WORK/annot/genome.sizes"
MAPPA="${BINS_MAPPABILITY:-}"
if [[ -z "$MAPPA" ]]; then
    for f in "$WORK"/qc/mappability_k*.bedgraph.gz; do if [[ -s "$f" ]]; then MAPPA="$f"; break; fi; done
fi
if [[ -z "$MAPPA" || ! -s "$MAPPA" ]]; then
    log "ERROR: no per-base mappability track. Run step 14 with qc.write_mappability_track: true, or set bins.mappability_bedgraph"
    exit 1
fi

# ---------------------------------------------------------------- 1. annotation (cached per bin size and inputs)
sig="w=$W;map=$(stat -c %Y "$MAPPA");cen=$(md5sum < "$WORK/annot/censat.class.bed" 2>/dev/null | cut -c1-12);rm=$(md5sum < "$WORK/annot/repeatmasker.class.bed" 2>/dev/null | cut -c1-12);fa=$GENOME_FASTA;sz=$(md5sum < "$SIZES" | cut -c1-12)"
if [[ ! -s "$D/manifest.tsv" || ! -f "$D/.annot.sig" || "$(cat "$D/.annot.sig")" != "$sig" ]]; then
    log "bin annotation at $W bp"
    python3 "$ROOT/bin/bin_table.py" annotate --sizes "$SIZES" --bin "$W" --out "$B" \
        --censat "$WORK/annot/censat.class.bed" --rmsk "$WORK/annot/repeatmasker.class.bed" --mappability "$MAPPA" \
        --fasta "$GENOME_FASTA" --fasta-names "$WORK/annot/fasta_names.tsv"
    printf "%s\n" "$sig" > "$D/.annot.sig"
else
    log "bin annotation at $W bp: reusing $D"
fi

# ---------------------------------------------------------------- 2. motif sites, genome-wide (one FIMO run per family)
REPS="$WORK/motifs/family_reps.meme"
[[ -s "$REPS" ]] || { log "ERROR: $REPS missing; run steps 02-04 first"; exit 1; }
export ROOT B FIMO_THRESH MEME_PREFIX GENOME_FASTA
python3 - <<'PY'
import os, sys, hashlib
sys.path.insert(0, os.path.join(os.environ["ROOT"], "bin"))
from memelib import read_meme, write_meme
B = os.environ["B"]
bg, ms = read_meme(os.path.join(os.path.dirname(B), "motifs", "family_reps.meme"))
os.makedirs(os.path.join(B, "fimo", "motifs"), exist_ok=True)
fa = os.environ["GENOME_FASTA"]; th = os.environ.get("FIMO_THRESH", "1e-4")
with open(os.path.join(B, "fimo", "families.tsv"), "w") as out:
    out.write("family\tkey\tmeme\n")
    for m in ms:
        key = hashlib.md5((repr(m.rows) + fa + str(os.path.getsize(fa)) + th).encode()).hexdigest()[:12]
        f = os.path.join(B, "fimo", "motifs", f"{m.id}.meme")
        write_meme(f, [m], bg)
        out.write(f"{m.id}\t{key}\t{f}\n")
PY
export WORK
fimo_one() {
    set -euo pipefail
    local fam="$1" key="$2" meme="$3" out="$B/fimo/$1.$2.tsv.gz"
    [[ -s "$out" ]] && return 0
    # FIMO text mode streams every site (no stored-score cap); --skip-matched-sequence keeps the output small
    $MEME_PREFIX fimo --text --skip-matched-sequence $(pgc_flag fimo) --verbosity 1 --thresh "$FIMO_THRESH" "$meme" "$GENOME_FASTA" \
        | gzip -c > "$out.tmp" && mv "$out.tmp" "$out"
}
export -f fimo_one
n_fam=$(($(wc -l < "$B/fimo/families.tsv") - 1))
log "genome-wide FIMO for $n_fam families (cached in $B/fimo)"
tail -n +2 "$B/fimo/families.tsv" | xargs -P "$THREADS" -n 3 bash -c 'fimo_one "$@"' _
mapfile -t FIMO_FILES < <(tail -n +2 "$B/fimo/families.tsv" | awk -v d="$B/fimo" '{print d "/" $1 "." $2 ".tsv.gz"}')
if [[ ! -s "$D/sites.tsv.gz" ]] || for f in "${FIMO_FILES[@]}"; do [[ "$f" -nt "$D/sites.tsv.gz" ]] && echo new; done | grep -q new; then
    python3 "$ROOT/bin/bin_table.py" sites --sizes "$SIZES" --bin "$W" --out "$B" --alias "$WORK/annot/chrom_alias.tsv" \
        --fimo "${FIMO_FILES[@]}"
fi

# ---------------------------------------------------------------- 3. signal per library (IP and control BAMs, once each)
{ sample_col bam_file | awk 'BEGIN{OFS="\t"} {print $1, "IP", $2}'
  sample_col control_bam | awk 'BEGIN{OFS="\t"} {print $1, "control", $2}'; } > "$D/.libs.tmp"
lib_name() { basename "$1" | sed -E 's/\.bam$//; s/\.mLb\.clN(\.sorted)?$//'; }
export -f lib_name
signal_one() {
    set -euo pipefail
    local bam="$1" name; name=$(lib_name "$bam")
    local o="$D/$name"
    if [[ -s "$o/library.json" && ! "$bam" -nt "$o/library.json" ]]; then return 0; fi
    python3 "$ROOT/bin/bin_table.py" signal --sizes "$SIZES" --bin "$W" --out "$B" --alias "$WORK/annot/chrom_alias.tsv" \
        --bam "$bam" --name "$name" --unique-mapq "$BINS_UNIQUE_MAPQ" --insert-max "$BINS_INSERT_MAX" --threads 2
}
export -f signal_one
export D W SIZES BINS_UNIQUE_MAPQ BINS_INSERT_MAX
cut -f3 "$D/.libs.tmp" | sort -u | xargs -P "$THREADS" -I{} bash -c 'signal_one "$@"' _ {}
{ printf "sample_id\tlibrary\tname\tdir\n"
  while IFS=$'\t' read -r s lib bam; do n=$(lib_name "$bam"); printf "%s\t%s\t%s\t%s\n" "$s" "$lib" "$n" "$D/$n"; done < "$D/.libs.tmp"; } > "$D/libraries.tsv"
rm -f "$D/.libs.tmp"
log "bin signal: $(($(cut -f3 "$D/libraries.tsv" | sort -u | wc -l) - 1)) libraries at $W bp"

# ---------------------------------------------------------------- 4. drift check against the tuning record
if [[ "${BINS_DRIFT_CHECK:-1}" == "1" && -n "${BINS_RECORD:-}" ]]; then
    floor=$(awk -F'\t' '$1 == "floor_size" {print $2}' "$BINS_RECORD")
    dchr=$(awk -F'\t' '$1 == "drift_chromosome" {print $2}' "$BINS_RECORD")
    if [[ -n "$floor" && -n "$dchr" ]] && grep -qx "$dchr" <(cut -f1 "$SIZES"); then
        g=$(python3 -c "import math,sys; print(math.gcd(int(sys.argv[1]), int(sys.argv[2])))" "$floor" "$W")   # both sizes are sums of these
        DF="$B/drift"
        log "drift check: $dchr at the tuning floor ($floor bp) vs $W bp (computed at $g bp)"
        python3 "$ROOT/bin/bin_table.py" annotate --sizes "$SIZES" --bin "$g" --chroms "$dchr" --out "$DF" \
            --censat "$WORK/annot/censat.class.bed" --rmsk "$WORK/annot/repeatmasker.class.bed"
        drift_one() {
            set -euo pipefail
            local bam="$1" name; name=$(lib_name "$bam")
            [[ -s "$DF/w$g/$name/library.json" && ! "$bam" -nt "$DF/w$g/$name/library.json" ]] && return 0
            python3 "$ROOT/bin/bin_table.py" signal --sizes "$SIZES" --bin "$g" --chroms "$dchr" --out "$DF" \
                --alias "$WORK/annot/chrom_alias.tsv" --bam "$bam" --name "$name" --unique-mapq "$BINS_UNIQUE_MAPQ" --insert-max "$BINS_INSERT_MAX"
        }
        export -f drift_one; export DF floor dchr g
        { sample_col bam_file | cut -f2; sample_col control_bam | cut -f2; } | sort -u | xargs -P "$THREADS" -I{} bash -c 'drift_one "$@"' _ {}
        Rscript "$ROOT/R/16_drift_check.R" "$CONFIG_YAML" "$DF/w$g" "$D/libraries.tsv" "$BINS_RECORD" "$W"
    else
        log "drift check skipped: the tuning record has no floor_size/drift_chromosome, or that chromosome is not in genome.sizes"
    fi
fi
log "bin table done: $D"
