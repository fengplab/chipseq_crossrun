#!/usr/bin/env bash
# Bin-size tuning: run once, before the first run that uses the binding models; every later run reuses the result
# through bins.tuning_record (step 16 then checks for drift on tuning.drift_chromosome).
# Usage: tune/run_tuning.sh config.yaml
#   1. setup and input preparation (run_all.sh steps 0-1) into <tuning.outdir>/setup
#   2. bin tables at tuning.base bp (default 10) on tuning.chromosomes only: annotation, and per-library signal
#      for every IP and control BAM in the manifest (in parallel, THREADS at a time)
#   3. optional: FIMO sites of tuning.reference_motifs (a MEME file) on the tuning chromosomes
#   4. R/tune_bins.R: signal retained and split-half correlation at every size in tuning.sizes -> chosen size
# Output (tuning.outdir; relative paths are under the config's outdir):
#   bin_size.tsv           the tuning record: chosen_size, floor_size, tolerance, ... (set bins.tuning_record to this)
#   tuning_curves.tsv      per sample, group (satellite = CenSat bins, other) and size
#   tuning_summary.tsv     medians per size and group, and whether each size passes
#   tuning_motif_effects.tsv  (with reference motifs) per family, sample and size
#   tuning_curves.png
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CFG="$(cd "$(dirname "${1:?config.yaml}")" && pwd)/$(basename "$1")"
read -r OUT TO < <(python3 - "$CFG" <<'PY'
import os, sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
out = os.path.abspath(c.get("outdir", "crossrun_out"))
to = (c.get("tuning") or {}).get("outdir") or "bin_tuning"
print(out, to if os.path.isabs(to) else os.path.join(out, to))
PY
)
mkdir -p "$TO"
TCFG="$TO/tuning_config.yaml"
python3 - "$CFG" "$TCFG" "$TO/setup" <<'PY'
import sys, yaml
c = yaml.safe_load(open(sys.argv[1]))
c["outdir"] = sys.argv[3]
yaml.safe_dump(c, open(sys.argv[2], "w"), sort_keys=False)
PY
echo "[tuning] setup in $TO/setup" >&2
bash "$ROOT/run_all.sh" "$TCFG" 0 1
ENV="$TO/setup/work/config.env"
source "$ENV"; source "$ROOT/bin/lib.sh"
genome_file >/dev/null
SIZES="$WORK/annot/genome.sizes"

eval "$(python3 - "$CFG" <<'PY'
import sys, yaml, shlex
t = (yaml.safe_load(open(sys.argv[1])).get("tuning") or {})
b = (yaml.safe_load(open(sys.argv[1])).get("bins") or {})
sizes = [int(x) for x in t.get("sizes", [60, 100, 150, 200, 250, 300, 400, 500, 750, 1000, 2000])]
base = int(t.get("base", 10))
bad = [s for s in sizes if s % base]
if bad:
    sys.exit(f"[tuning] sizes must be multiples of tuning.base ({base}): {bad}")
print("TCHROMS=" + shlex.quote(" ".join(t.get("chromosomes", ["chr1", "chr9", "chr16", "chr21"]))))
print(f"TBASE={base}")
print("TREF=" + shlex.quote(t.get("reference_motifs") or ""))
print("TMAP=" + shlex.quote(b.get("mappability_bedgraph") or ""))
PY
)"
for c in $TCHROMS; do
    grep -qx "$c" <(cut -f1 "$SIZES") || { log "ERROR: tuning chromosome $c is not a host nuclear chromosome in $SIZES"; exit 1; }
done
# mappability is optional here (used only to stratify the reference motif effects): bins.mappability_bedgraph,
# else the per-base track that step 14 of the main run wrote
if [[ -z "$TMAP" ]]; then
    for f in "$OUT"/work/qc/mappability_k*.bedgraph.gz; do if [[ -s "$f" ]]; then TMAP="$f"; break; fi; done
fi
B="$TO/bins"; D="$B/w$TBASE"
read -r -a CH <<< "$TCHROMS"

log "tuning annotation at $TBASE bp on ${CH[*]}${TMAP:+ (mappability: $TMAP)}"
python3 "$ROOT/bin/bin_table.py" annotate --sizes "$SIZES" --bin "$TBASE" --chroms "${CH[@]}" --out "$B" \
    --censat "$WORK/annot/censat.class.bed" --rmsk "$WORK/annot/repeatmasker.class.bed" \
    ${TMAP:+--mappability "$TMAP"} --fasta "$GENOME_FASTA" --fasta-names "$WORK/annot/fasta_names.tsv"

{ sample_col bam_file | awk 'BEGIN{OFS="\t"} {print $1, "IP", $2}'
  sample_col control_bam | awk 'BEGIN{OFS="\t"} {print $1, "control", $2}'; } > "$D/.libs.tmp"
lib_name() { basename "$1" | sed -E 's/\.bam$//; s/\.mLb\.clN(\.sorted)?$//'; }
tune_signal() {
    set -euo pipefail
    local bam="$1" name; name=$(lib_name "$bam")
    if [[ -s "$D/$name/library.json" && ! "$bam" -nt "$D/$name/library.json" ]]; then return 0; fi
    python3 "$ROOT/bin/bin_table.py" signal --sizes "$SIZES" --bin "$TBASE" --chroms $TCHROMS --out "$B" \
        --alias "$WORK/annot/chrom_alias.tsv" --bam "$bam" --name "$name" \
        --unique-mapq "$BINS_UNIQUE_MAPQ" --insert-max "$BINS_INSERT_MAX"
}
export -f lib_name tune_signal
export ROOT B D TBASE TCHROMS SIZES WORK BINS_UNIQUE_MAPQ BINS_INSERT_MAX
n_lib=$(cut -f3 "$D/.libs.tmp" | sort -u | wc -l)
log "tuning signal: $n_lib libraries at $TBASE bp, $THREADS at a time"
cut -f3 "$D/.libs.tmp" | sort -u | xargs -P "$THREADS" -I{} bash -c 'tune_signal "$@"' _ {}
{ printf "sample_id\tlibrary\tname\tdir\n"
  while IFS=$'\t' read -r s lib bam; do n=$(lib_name "$bam"); printf "%s\t%s\t%s\t%s\n" "$s" "$lib" "$n" "$D/$n"; done < "$D/.libs.tmp"; } > "$D/libraries.tsv"
rm -f "$D/.libs.tmp"

if [[ -n "$TREF" ]]; then
    [[ -s "$TREF" ]] || { log "ERROR: tuning.reference_motifs not found: $TREF"; exit 1; }
    FA="$TO/tuning_chroms.fa"
    if [[ ! -s "$FA" || "$GENOME_FASTA" -nt "$FA" ]]; then
        # FASTA names of the tuning chromosomes (fasta_names.tsv: canonical -> FASTA name)
        mapfile -t FN < <(awk -F'\t' -v c=" $TCHROMS " 'index(c, " " $1 " ") {print $2}' "$WORK/annot/fasta_names.tsv")
        [[ ${#FN[@]} -eq ${#CH[@]} ]] || { log "ERROR: not every tuning chromosome is in $WORK/annot/fasta_names.tsv"; exit 1; }
        samtools faidx "$GENOME_FASTA" "${FN[@]}" > "$FA.tmp" && mv "$FA.tmp" "$FA"
    fi
    log "FIMO of the reference motifs on the tuning chromosomes"
    $MEME_PREFIX fimo --text --skip-matched-sequence $(pgc_flag fimo) --verbosity 1 --thresh "$FIMO_THRESH" "$TREF" "$FA" \
        | gzip -c > "$TO/ref_fimo.tsv.gz.tmp"
    mv "$TO/ref_fimo.tsv.gz.tmp" "$TO/ref_fimo.tsv.gz"
    python3 "$ROOT/bin/bin_table.py" sites --sizes "$SIZES" --bin "$TBASE" --chroms "${CH[@]}" --out "$B" \
        --alias "$WORK/annot/chrom_alias.tsv" --fimo "$TO/ref_fimo.tsv.gz" --name ref_sites.tsv.gz
fi

Rscript "$ROOT/R/tune_bins.R" "$TCFG" "$TO"
log "tuning done. Set this in the config of every run that uses the binding models:"
log "  bins: {tuning_record: $TO/bin_size.tsv}"
