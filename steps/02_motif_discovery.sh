#!/usr/bin/env bash
# Step 02: motif discovery.
#   motifs.discovery_engine = text     -> STREME --text, MEME -text and SEA --text run directly (no HTML
#                                         output, so a MEME install with broken HTML templates still works)
#   motifs.discovery_engine = xstreme  -> XSTREME (always writes an HTML report)
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

# STREME's own width limits, from `streme --help` (MEME Suite 5.5.9): --minw must be >= 3 and --maxw must be
# <= 30. MEME has no width cap. One helper decides for both engines, so a width range that STREME cannot
# search is handled the same way everywhere instead of failing inside the tool.
STREME_MAX_WIDTH=30

# streme_maxw <minw> <maxw> -> the --maxw to hand STREME, or "skip" when STREME cannot search this range at all
streme_maxw() {
    local mn="$1" mx="$2"
    if (( mn > STREME_MAX_WIDTH )); then echo skip
    elif (( mx > STREME_MAX_WIDTH )); then echo "$STREME_MAX_WIDTH"
    else echo "$mx"; fi
}
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export ROOT

if [[ "$DISCOVERY" == "pipeline" ]]; then
    log "motif discovery: reusing the pipeline's XSTREME results"
    python3 "$ROOT/bin/02_parse_xstreme.py" "$SAMPLES" "$OUTDIR"
    exit 0
fi

bash "$ROOT/steps/03_sequences.sh" "$1"          # peak windows (cached; step 03 reuses them)
mkdir -p "$WORK/control_seqs" "$WORK/xstreme"

# --- installation check on a tiny input, before many long runs ----------------------------------------
#   text engine: STREME --text, MEME -text and SEA --text never touch the HTML templates
#   xstreme engine: XSTREME always writes an HTML report, so HTML output must work
stamp="$WORK/xstreme/.meme_ok_$DISCOVERY_ENGINE"
if [[ "${CROSSRUN_SKIP_MEME_CHECK:-0}" != "1" && ! -f "$stamp" ]]; then
    t=$(mktemp -d)
    python3 - "$t" <<'PY'
import random, sys
random.seed(1); t = sys.argv[1]
with open(f"{t}/p.fa", "w") as p, open(f"{t}/n.fa", "w") as n:
    for i in range(40):
        s = "".join(random.choice("ACGT") for _ in range(60))
        p.write(f">p{i}\n{s[:25]}TTGACGTCAA{s[35:]}\n"); n.write(f">n{i}\n{''.join(random.choice('ACGT') for _ in range(60))}\n")
open(f"{t}/m.meme", "w").write("MEME version 5\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\n"
    "A 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF m1\nletter-probability matrix: alength= 4 w= 4 nsites= 20 E= 0\n"
    "0.97 0.01 0.01 0.01\n0.01 0.97 0.01 0.01\n0.01 0.01 0.97 0.01\n0.01 0.01 0.01 0.97\n\n")
PY
    ok=1
    if [[ "$DISCOVERY_ENGINE" == "xstreme" ]]; then
        $MEME_PREFIX sea --p "$t/p.fa" --m "$t/m.meme" --oc "$t/out" --verbosity 1 > "$t/log" 2>&1 && [[ -s "$t/out/sea.html" ]] || ok=0
    else
        { $MEME_PREFIX fasta-get-markov -m 2 -nostatus -nosummary "$t/n.fa" "$t/bg" \
          && $MEME_PREFIX streme --text --dna --p "$t/p.fa" --n "$t/n.fa" --nmotifs 1 --verbosity 1 > "$t/s.txt" \
          && $MEME_PREFIX meme "$t/p.fa" -dna -text -revcomp -nmotifs 1 -minw 6 -maxw 10 -nostatus > "$t/m.txt" \
          && $MEME_PREFIX sea --text --p "$t/p.fa" --m "$t/m.meme" --verbosity 1 > "$t/sea.tsv"; } > "$t/log" 2>&1 || ok=0
        grep -q "^MOTIF" "$t/s.txt" 2>/dev/null || ok=0
    fi
    if [[ $ok == 1 ]]; then
        touch "$stamp"; rm -rf "$t"
    else
        log "ERROR: MEME Suite check failed for the '$DISCOVERY_ENGINE' discovery engine:"
        tail -5 "$t/log" >&2
        if [[ "$DISCOVERY_ENGINE" == "xstreme" ]]; then
            log "XSTREME needs working HTML output. Use motifs.discovery_engine: text (no HTML needed),"
        else
            log "Text-mode STREME/MEME/SEA failed, so the problem is not the HTML templates."
        fi
        log "fix the MEME install (conda install -c bioconda -c conda-forge --force-reinstall meme=5.5.9),"
        log "set meme_prefix to a working MEME container, or set motifs.discovery: pipeline."
        rm -rf "$t"; exit 1
    fi
fi

control_sequences() {   # <sample> -> work/control_seqs/<sample>.fa ; prints the control type
    set -euo pipefail
    local s="$1" bam scan="$WORK/seqs/$1.scan.bed" o="$WORK/control_seqs/$1" n want total frac pe flt
    bam=$(sample_get "$s" control_bam)
    if [[ -z "$bam" ]]; then echo shuffled; return 0; fi
    if [[ ! -r "$bam" ]]; then
        log "ERROR: control BAM for $s cannot be read: $bam (if results moved, update the manifest and rerun step 0)"
        return 1
    fi
    if [[ -s "$o.fa" && "$o.fa" -nt "$scan" ]]; then echo control_bam; return 0; fi
    n=$(wc -l < "$scan"); [[ "$n" -gt 0 ]] || { echo shuffled; return 0; }
    want=$(awk -v n="$n" -v r="$CONTROL_RATIO" 'BEGIN{printf "%d", n*r + 0.5}')
    pe=0; [[ $(samtools view -f 1 "$bam" | head -n 1000 | wc -l) -gt 0 ]] && pe=1
    # reads on the host nuclear chromosomes only (viral genomes and chrM are not sampled), so the draw fraction is right
    total=$(host_total "$bam" reads 2>/dev/null || echo 0)
    [[ "$total" -gt 0 ]] || total=$(samtools idxstats "$bam" 2>/dev/null | awk '{t+=$3} END{print t+0}')
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

# Text-mode chain mirroring XSTREME's discovery: order-2 background from the controls (as XSTREME does),
# STREME and MEME (differential enrichment) against the same controls, then SEA to score every motif.
# Each command checks its own status (errexit is not active inside a function called from 'if !').
text_discovery() {
    local s="$1" q="$2" out="$3" ctype="$4" maxw="$5" smaxw="$6" bgsrc="$2" n_arg="" neg_arg=""
    if [[ $ctype == control_bam ]]; then
        bgsrc="$WORK/control_seqs/$s.fa"; n_arg="--n $bgsrc"; neg_arg="-neg $bgsrc"
    fi
    $MEME_PREFIX fasta-get-markov -m 2 -pseudo 1 -nostatus -nosummary "$bgsrc" "$out/background.txt" || return 1
    # STREME only when it can search the range (see streme_maxw). assemble_discovery.py treats a missing
    # streme.txt as "no STREME motifs", so MEME alone carries discovery and nothing downstream changes shape.
    if [[ "$smaxw" == skip ]]; then
        rm -f "$out/streme.txt"
        echo "[discovery] STREME skipped: motifs.min_width $MOTIF_MIN_WIDTH is above STREME's maximum width" \
             "$STREME_MAX_WIDTH; discovering with MEME alone over ${MOTIF_MIN_WIDTH}-${maxw} bp" >&2
    else
        # shellcheck disable=SC2086
        $MEME_PREFIX streme --text --dna --p "$q" $n_arg --minw "$MOTIF_MIN_WIDTH" --maxw "$smaxw" --bfile "$out/background.txt" \
            --verbosity 1 $STREME_EXTRA_ARGS > "$out/streme.txt" || return 1
    fi
    # shellcheck disable=SC2086
    $MEME_PREFIX meme "$q" -dna -text -revcomp -mod zoops -objfun de $neg_arg -minw "$MOTIF_MIN_WIDTH" -maxw "$maxw" \
        -nmotifs "$MEME_NMOTIFS" -searchsize "$MEME_SEARCHSIZE" -bfile "$out/background.txt" -nostatus $MEME_EXTRA_ARGS > "$out/meme.txt" || return 1
    python3 "$ROOT/bin/assemble_discovery.py" --streme "$out/streme.txt" --meme "$out/meme.txt" --out "$out" || return 1
    if grep -q "^MOTIF" "$out/combined.meme"; then
        # shellcheck disable=SC2086
        $MEME_PREFIX sea --text --p "$q" $n_arg --m "$out/combined.meme" --bfile "$out/background.txt" \
            --thresh 1e6 --verbosity 1 > "$out/sea_all.tsv" || return 1
        grep -q "^RANK" "$out/sea_all.tsv" || { echo "SEA output has no header" >&2; return 1; }
    else
        printf "RANK\tDB\tID\tALT_ID\tCONSENSUS\tTP\tTP%%\tFP\tFP%%\tENR_RATIO\tSCORE_THR\tPVALUE\tLOG_PVALUE\tEVALUE\tLOG_EVALUE\tQVALUE\tLOG_QVALUE\n" > "$out/sea_all.tsv"
    fi
    python3 "$ROOT/bin/assemble_discovery.py" --streme "$out/streme.txt" --meme "$out/meme.txt" --sea "$out/sea_all.tsv" --out "$out" || return 1
}

discover_sample() {
    set -euo pipefail
    local s="$1" q="$WORK/seqs/$1.fa" out="$WORK/xstreme/$1_xstreme" ctype maxw m_args="" n_ctrl=0
    # reuse finished discovery unless the sequences, the engine or the minimum motif width changed
    # (results from before motifs.min_width existed were made with the old fixed minimum of 6)
    if [[ -s "$out/xstreme.tsv" && -s "$out/control_info.tsv" && ! "$q" -nt "$out/xstreme.tsv" ]] \
       && [[ "$(awk -F'\t' 'NR==2{print $4}' "$out/control_info.tsv")" == "$DISCOVERY_ENGINE" ]] \
       && [[ "$(awk -F'\t' 'NR==2{print ($5 == "" ? 6 : $5)}' "$out/control_info.tsv")" == "$MOTIF_MIN_WIDTH" ]]; then return 0; fi
    ctype=$(control_sequences "$s") || return 1
    # The motif width range comes from the config and nothing else.
    #
    # This function only runs when motifs.discovery is control_bam -- the 'pipeline' branch returns at the top
    # of this file -- so crossrun is running its own STREME and MEME here and has no reason to inherit a width
    # from the pipeline's XSTREME. It used to read <xstreme_dir>/pipeline_maxw.txt and let that REPLACE
    # motifs.xstreme_maxw. When that file said 20 and motifs.min_width was 20, the search collapsed to a single
    # width, every motif came out exactly 20 bp, and neither branch below logged anything about it.
    local maxw_src="motifs.xstreme_maxw"
    maxw=$XSTREME_MAXW
    if (( maxw < MOTIF_MIN_WIDTH )); then
        log "WARNING: $s: motifs.xstreme_maxw ($maxw) is below motifs.min_width $MOTIF_MIN_WIDTH;" \
            "using $MOTIF_MIN_WIDTH for both"
        maxw=$MOTIF_MIN_WIDTH
    fi
    if (( maxw == MOTIF_MIN_WIDTH )); then
        log "WARNING: $s: the motif width search is pinned to $MOTIF_MIN_WIDTH bp exactly (motifs.min_width ==" \
            "motifs.xstreme_maxw, from $maxw_src). Every motif discovered will be $MOTIF_MIN_WIDTH bp long."
    fi
    local smaxw; smaxw=$(streme_maxw "$MOTIF_MIN_WIDTH" "$maxw")
    if [[ "$smaxw" == skip ]]; then
        log "$s: motif width ${MOTIF_MIN_WIDTH}-${maxw} bp (max from $maxw_src); STREME skipped (its maximum is $STREME_MAX_WIDTH bp), MEME alone"
    elif (( smaxw < maxw )); then
        log "$s: motif width ${MOTIF_MIN_WIDTH}-${maxw} bp (max from $maxw_src); STREME searches ${MOTIF_MIN_WIDTH}-${smaxw} bp (its maximum), MEME the full range"
    fi
    for d in $KNOWN_DBS; do m_args="$m_args --m $d"; done
    local ctrl_arg=""
    if [[ $ctype == control_bam ]]; then
        ctrl_arg="--n $WORK/control_seqs/$s.fa"; n_ctrl=$(grep -c '^>' "$WORK/control_seqs/$s.fa")
    else
        log "WARNING: $s has no control BAM; XSTREME uses shuffled peak sequences as controls"
    fi
    rm -rf "$out"; mkdir -p "$out"
    local L="$WORK/xstreme/$s.log"
    if [[ "$DISCOVERY_ENGINE" == "xstreme" ]]; then
        if [[ "$smaxw" == skip ]]; then
            log "ERROR: $s: motifs.min_width $MOTIF_MIN_WIDTH is above STREME's maximum width $STREME_MAX_WIDTH, and"
            log "       motifs.discovery_engine 'xstreme' runs STREME internally, so discovery cannot run."
            log "       Set motifs.discovery_engine: text, which uses MEME alone above $STREME_MAX_WIDTH bp."
            return 1
        fi
        # shellcheck disable=SC2086
        if ! $MEME_PREFIX xstreme --p "$q" $ctrl_arg --minw "$MOTIF_MIN_WIDTH" --maxw "$maxw" $m_args --fimo-skip --sea-noseqs \
                --verbosity 1 --oc "$out" $XSTREME_EXTRA_ARGS > "$L" 2>&1 || [[ ! -s "$out/xstreme.tsv" ]]; then
            log "ERROR: XSTREME failed for $s; see $L"; tail -5 "$L" >&2; return 1
        fi
    else
        if ! text_discovery "$s" "$q" "$out" "$ctype" "$maxw" "$smaxw" > "$L" 2>&1; then
            log "ERROR: motif discovery (text mode) failed for $s; see $L"; tail -5 "$L" >&2; return 1
        fi
    fi
    [[ -s "$out/combined.meme" || -s "$out/xstreme.txt" ]] || { log "ERROR: motif discovery for $s wrote no motif file"; return 1; }
    # Kept under this name because bin/02_parse_xstreme.py reads it for xstreme_run_info.maxw_used. It records
    # the width THIS step used; it is never read back as an input (see the width block above).
    printf "maxw_used\t%s\n" "$maxw" > "$out/pipeline_maxw.txt"
    printf "control_type\tn_control_seqs\tn_primary_seqs\tengine\tmin_width\n%s\t%s\t%s\t%s\t%s\n" "$ctype" "$n_ctrl" "$(grep -c '^>' "$q")" "$DISCOVERY_ENGINE" "$MOTIF_MIN_WIDTH" > "$out/control_info.tsv"
    log "motif discovery ($DISCOVERY_ENGINE): $s ($ctype controls: $n_ctrl; $(grep -c '^MOTIF' "$out/combined.meme" 2>/dev/null || echo 0) motifs)"
}
# discover_sample runs in an xargs subshell, so everything it calls must cross into it: streme_maxw is a
# function and STREME_MAX_WIDTH a plain variable, and neither is inherited without being exported.
export -f control_sequences discover_sample text_discovery streme_maxw
export STREME_MAX_WIDTH
genome_file >/dev/null
cut -f1 <(sample_col peaks_file) | xargs -P "$THREADS" -I{} bash -c 'discover_sample "$@"' _ {}
python3 "$ROOT/bin/02_parse_xstreme.py" "$SAMPLES" "$OUTDIR" "$WORK/xstreme"
