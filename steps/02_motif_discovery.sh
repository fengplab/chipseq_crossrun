#!/usr/bin/env bash
# Step 02: motif discovery.
#   motifs.discovery = pipeline     -> parse the pipeline's XSTREME results (shuffled-sequence controls)
#   motifs.discovery = control_bam  -> rerun XSTREME per sample with control sequences from the sample's
#                                      control (input/IgG) alignment, then parse those results
# Control sequences: control-library fragments are sampled from the control BAM (properly paired read-1
# records; read alignments for single-end data). Windows are centred on the fragment midpoints, with widths
# drawn from the sample's own peak windows (whole peaks or summit +/- flank, as in step 03), so primary
# and control sequences have matching length distributions. Windows overlapping the sample's peaks are
# discarded. Sequences come from the genome FASTA at those positions.
#   work/control_seqs/<sample>.bed / .fa      control windows and sequences
#   work/xstreme/<sample>_xstreme/            XSTREME output (+ control_info.tsv)
source "$1"; source "$(dirname "$0")/../bin/lib.sh"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export ROOT

if [[ "$DISCOVERY" == "pipeline" ]]; then
    log "motif discovery: reusing the pipeline's XSTREME results"
    python3 "$ROOT/bin/02_parse_xstreme.py" "$SAMPLES" "$OUTDIR"
    exit 0
fi

bash "$ROOT/steps/03_sequences.sh" "$1"          # peak windows (cached; step 03 reuses them)
mkdir -p "$WORK/control_seqs" "$WORK/xstreme"

# --- XSTREME always writes an HTML report: make sure this MEME installation can, before 24 long runs --------
if [[ "${CROSSRUN_SKIP_MEME_CHECK:-0}" != "1" && ! -f "$WORK/xstreme/.meme_html_ok" ]]; then
    t=$(mktemp -d)
    printf ">a\nACGTACGTTTGACGTCAACGTACGTAAAA\n>b\nTTGACGTCAAGGCCTTAACCGGTTGACGTCA\n" > "$t/p.fa"
    printf "MEME version 5\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF m1\nletter-probability matrix: alength= 4 w= 4 nsites= 20 E= 0\n0.97 0.01 0.01 0.01\n0.01 0.97 0.01 0.01\n0.01 0.01 0.97 0.01\n0.01 0.01 0.01 0.97\n\n" > "$t/m.meme"
    if $MEME_PREFIX sea --p "$t/p.fa" --m "$t/m.meme" --oc "$t/out" --verbosity 1 > "$t/log" 2>&1 && [[ -s "$t/out/sea.html" ]]; then
        touch "$WORK/xstreme/.meme_html_ok"; rm -rf "$t"
    else
        log "ERROR: this MEME installation cannot write HTML reports, which XSTREME requires:"
        tail -3 "$t/log" >&2
        log "Fix the MEME install (e.g. conda install -c bioconda -c conda-forge --force-reinstall meme=5.5.9),"
        log "or set meme_prefix in config.yaml to run MEME from a working container (e.g. the one the pipeline used),"
        log "or set motifs.discovery: pipeline to reuse the pipeline's XSTREME results."
        rm -rf "$t"; exit 1
    fi
fi

control_sequences() {   # <sample> -> work/control_seqs/<sample>.fa ; prints the control type
    set -euo pipefail
    local s="$1" bam scan="$WORK/seqs/$1.scan.bed" o="$WORK/control_seqs/$1" n want total frac pe flt
    bam=$(sample_get "$s" control_bam)
    if [[ -z "$bam" || ! -s "$bam" ]]; then echo shuffled; return 0; fi
    if [[ -s "$o.fa" && "$o.fa" -nt "$scan" ]]; then echo control_bam; return 0; fi
    n=$(wc -l < "$scan"); [[ "$n" -gt 0 ]] || { echo shuffled; return 0; }
    want=$(awk -v n="$n" -v r="$CONTROL_RATIO" 'BEGIN{printf "%d", n*r + 0.5}')
    pe=0; [[ $(samtools view -f 1 "$bam" | head -n 1000 | wc -l) -gt 0 ]] && pe=1
    total=$(samtools idxstats "$bam" 2>/dev/null | awk '{t+=$3} END{print t+0}')
    [[ "$total" -gt 0 ]] || total=$(samtools view -c -F 3844 "$bam")
    # oversample ~6x (fragments are lost to peak overlap and chromosome filters), then draw exactly 'want'
    frac=$(awk -v w="$want" -v t="$total" -v pe="$pe" 'BEGIN{f=6*w*(pe?2:1)/t; if (f>=1) print ""; else printf "11%s", substr(sprintf("%.9f", f), 2)}')
    flt="-F 3844"; [[ $pe == 1 ]] && flt="-f 66 -F 3852"
    awk '{print $3-$2}' "$scan" | shuf --random-source=<(yes 7) > "$o.widths"
    samtools view ${frac:+-s $frac} $flt "$bam" \
        | awk -v pe="$pe" -v mx="$QC_INSERT_MAX" 'BEGIN{OFS="\t"}
            pe==1 { if ($7 != "=" || $9 == 0) next; t = ($9 < 0) ? -$9 : $9; if (t > mx) next
                    l = (($4 < $8) ? $4 : $8) - 1; m = l + int(t / 2); print $3, m, m + 1; next }
            pe==0 { m = $4 - 1 + int(length($10) / 2); print $3, m, m + 1 }' \
        | from_any_names \
        | awk 'BEGIN{FS=OFS="\t"} NR==FNR{ok[$1]=1; next} ($1 in ok)' "$(genome_file)" - \
        | sort -k1,1 -k2,2n | bedtools intersect -v -a stdin -b "$scan" \
        | shuf --random-source=<(yes 11) | awk -v n="$want" 'NR <= n' \
        | awk 'BEGIN{FS=OFS="\t"} NR==FNR { w[++nw]=$1; next } FILENAME==ARGV[2] { L[$1]=$2; next }
               { i++; wd = w[(i-1) % nw + 1]; st = $2 - int(wd/2); if (st < 0) st = 0; en = st + wd
                 if (en > L[$1]) next; print $1, st, en, "ctrl|" i }' "$o.widths" "$(genome_file)" - \
        | sort -k1,1 -k2,2n > "$o.bed"
    getfasta_canonical "$o.bed" "$o.fa"
    rm -f "$o.widths"
    echo control_bam
}

discover_sample() {
    set -euo pipefail
    local s="$1" q="$WORK/seqs/$1.fa" out="$WORK/xstreme/$1_xstreme" ctype maxw px m_args="" n_ctrl=0
    if [[ -s "$out/xstreme.tsv" && -s "$out/control_info.tsv" && ! "$q" -nt "$out/xstreme.tsv" ]]; then return 0; fi
    ctype=$(control_sequences "$s")
    px=$(sample_get "$s" xstreme_dir)
    maxw=$( [[ -n "$px" && -s "$px/pipeline_maxw.txt" ]] && awk '{v=$NF} END{print v}' "$px/pipeline_maxw.txt" || echo "$XSTREME_MAXW")
    for d in $KNOWN_DBS; do m_args="$m_args --m $d"; done
    local ctrl_arg=""
    if [[ $ctype == control_bam ]]; then
        ctrl_arg="--n $WORK/control_seqs/$s.fa"; n_ctrl=$(grep -c '^>' "$WORK/control_seqs/$s.fa")
    else
        log "WARNING: $s has no control BAM; XSTREME uses shuffled peak sequences as controls"
    fi
    rm -rf "$out"
    # shellcheck disable=SC2086
    if ! $MEME_PREFIX xstreme --p "$q" $ctrl_arg --maxw "$maxw" $m_args --fimo-skip --sea-noseqs \
            --verbosity 1 --oc "$out" $XSTREME_EXTRA_ARGS > "$WORK/xstreme/$s.log" 2>&1 \
       || [[ ! -s "$out/xstreme.tsv" ]]; then
        log "ERROR: XSTREME failed for $s; see $WORK/xstreme/$s.log"; tail -5 "$WORK/xstreme/$s.log" >&2; return 1
    fi
    [[ -s "$out/combined.meme" || -s "$out/xstreme.txt" ]] || { log "ERROR: XSTREME for $s wrote no motif file"; return 1; }
    printf "maxw_used\t%s\n" "$maxw" > "$out/pipeline_maxw.txt"
    printf "control_type\tn_control_seqs\tn_primary_seqs\n%s\t%s\t%s\n" "$ctype" "$n_ctrl" "$(grep -c '^>' "$q")" > "$out/control_info.tsv"
    log "XSTREME: $s ($ctype controls: $n_ctrl)"
}
export -f control_sequences discover_sample
genome_file >/dev/null
cut -f1 <(sample_col peaks_file) | xargs -P "$THREADS" -I{} bash -c 'discover_sample "$@"' _ {}
python3 "$ROOT/bin/02_parse_xstreme.py" "$SAMPLES" "$OUTDIR" "$WORK/xstreme"
