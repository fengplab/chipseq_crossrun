#!/usr/bin/env bash
# Shared helpers for the bash steps. Source after <outdir>/work/config.env.
set -euo pipefail
export LC_ALL=C

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }


# Print "<sample_id>\t<value of column>" for every sample whose column is non-empty.
# Usage: sample_col peaks_file
sample_col() {
    awk -v col="$1" 'BEGIN{FS=OFS="\t"}
        NR==1 { for (i=1;i<=NF;i++) h[$i]=i; if (!(col in h)) exit 1; next }
        $h[col] != "" { print $h["sample_id"], $h[col] }' "$SAMPLES"
}

# Value of one column for one sample.
sample_get() {
    awk -v s="$1" -v col="$2" 'BEGIN{FS="\t"}
        NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
        $h["sample_id"]==s { print $h[col] }' "$SAMPLES"
}

# MACS3 narrowPeak/broadPeak -> BED5: chrom start end <sample>|<peak_name> summit_abs
# broadPeak has no summit column, so the midpoint is used.
peaks_to_bed5() {
    local sample="$1" peaks="$2" g
    g=$(genome_file)
    # chromosome names are made canonical first (e.g. NC_060925.1 -> chr1)
    awk -v s="$sample" 'BEGIN{FS=OFS="\t"}
        NR==FNR { ok[$1]=$2; next }
        !($1 in ok) { next }
        $3 > ok[$1] { $3 = ok[$1] }
        $2 >= $3 { next }
        $0 !~ /^(track|browser|#)/ && NF>=3 && $2 ~ /^[0-9]+$/ {
            n++; name = (NF>=4 && $4!="") ? $4 : "peak_" n
            summit = (NF>=10 && $10 ~ /^[0-9]+$/ && $10 >= 0) ? $2+$10 : int(($2+$3)/2)
            print $1, $2, $3, s "|" name, summit
        }' "$g" <(rename_chroms "peaks:$sample" < "$peaks") | sort -k1,1 -k2,2n
}

# Per-peak details for the peak overview (same peaks, names and clipping as peaks_to_bed5):
#   chrom start end width signal_value neglog10_p neglog10_q
# signal_value, -log10 p and -log10 q are columns 7-9 of narrowPeak/broadPeak files (for MACS, column 7 is the
# fold enrichment at the summit, or over the broad region); NA when the file has fewer columns or MACS wrote -1.
peaks_info() {
    local sample="$1" peaks="$2" g
    g=$(genome_file)
    awk 'BEGIN{FS=OFS="\t"}
        function v(x) { return (x ~ /^-?[0-9.]+([eE][-+]?[0-9]+)?$/ && x + 0 >= 0) ? x : "NA" }
        NR==FNR { ok[$1]=$2; next }
        !($1 in ok) { next }
        $3 > ok[$1] { $3 = ok[$1] }
        $2 >= $3 { next }
        $0 !~ /^(track|browser|#)/ && NF>=3 && $2 ~ /^[0-9]+$/ {
            print $1, $2, $3, $3 - $2, (NF>=7 ? v($7) : "NA"), (NF>=8 ? v($8) : "NA"), (NF>=9 ? v($9) : "NA")
        }' "$g" <(rename_chroms "peaks:$sample" < "$peaks") | sort -k1,1 -k2,2n
}

# ---------------------------------------------------------------------------------------------
# Chromosome naming. Every file the bundle writes uses canonical UCSC-style names (chr1..chrX, chrY,
# chrM). Inputs in any other convention (e.g. RefSeq NC_060925.1 for T2T-CHM13) are renamed on the
# way in with work/annot/chrom_alias.tsv, and renamed back only where a tool must read the original
# FASTA or BAM.
# ---------------------------------------------------------------------------------------------
LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export LIB_DIR

# Canonical chrom sizes (built once, together with the alias table)
genome_file() {
    local out="$WORK/annot/genome.sizes"
    if [[ ! -s "$out" || ! -s "$WORK/annot/chrom_alias.tsv" ]]; then
        mkdir -p "$WORK/annot"
        # shellcheck disable=SC2086
        python3 "$LIB_DIR/chrom_alias.py" --fasta "$GENOME_FASTA" --sizes "$CHROM_SIZES" \
            --outdir "$WORK/annot" ${CHROM_ALIAS_FILES:+--alias $CHROM_ALIAS_FILES} >&2
        sort -k1,1 -o "$out" "$out"
    fi
    echo "$out"
}

# Rename column 1 of a BED-like stream to canonical names.
# Usage: rename_chroms <label> < in > out
# Header/comment lines pass through. Names that are neither in the alias table nor chr-style are kept,
# counted and logged to work/annot/rename_log.tsv under <label>.
rename_chroms() {
    local label="$1"
    genome_file >/dev/null
    awk -v lab="$label" -v logf="$WORK/annot/rename_log.tsv" 'BEGIN{FS=OFS="\t"}
        NR==FNR { m[$1]=$2; next }
        /^(#|track|browser)/ { print; next }
        { n++
          if ($1 in m) { if (m[$1] != $1) r++; $1 = m[$1] }
          else if ($1 !~ /^chr/) { u[$1]++ }
          print }
        END { ex=""; k=0; for (x in u) { nu += u[x]; if (k++ < 3) ex = ex (ex==""?"":",") x }
              printf "%s\t%d\t%d\t%d\t%s\n", lab, n, r, nu, ex >> logf }' \
        "$WORK/annot/chrom_alias.tsv" -
}

# Canonical -> genome FASTA names (for bedtools getfasta / nuc), column 1 only.
to_fasta_names() {
    awk 'BEGIN{FS=OFS="\t"} NR==FNR { m[$1]=$2; next } ($1 in m) { $1 = m[$1] } { print }' \
        "$WORK/annot/fasta_names.tsv" -
}

# Any names (e.g. FASTA names coming back from bedtools nuc) -> canonical, column 1 only.
from_any_names() {
    awk 'BEGIN{FS=OFS="\t"} NR==FNR { m[$1]=$2; next } ($1 in m) { $1 = m[$1] } { print }' \
        "$WORK/annot/chrom_alias.tsv" -
}

# getfasta on a canonical BED; FASTA headers are exactly the BED name (sample|peak).
# bedtools -name appends "::chrom:start-end"; that suffix is removed because MEME Suite 5.5.4+ (FIMO, SEA, STREME)
# parses such headers as genomic coordinates and then reports the chromosome instead of the sequence name.
getfasta_canonical() {
    local bed="$1" out="$2"
    to_fasta_names < "$bed" | bedtools getfasta -fi "$GENOME_FASTA" -bed stdin -name \
        | awk '/^>/ { sub(/::.*$/, "") } { print }' > "$out"
}

# Version of the scanned-sequence format; stored next to outputs of steps 03, 05 and 16 so that a bundle update that
# changes how sequences are named or scanned reruns them instead of reusing stale results.
SCAN_VERSION=2

# "--no-pgc" for a MEME Suite tool (fimo, sea) whose version parses genomic coordinates in FASTA headers by default
# (FIMO >= 5.5.4, SEA >= 5.5.5); empty for older versions. Cached per tool in the environment.
pgc_flag() {
    local tool="$1" var="PGC_FLAG_${1^^}" v min
    if [[ -z "${!var+x}" ]]; then
        v=$($MEME_PREFIX "$tool" --version 2>/dev/null | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' || true)
        [[ "$tool" == sea ]] && min=5.5.5 || min=5.5.4
        if [[ -n "$v" && "$(printf '%s\n%s\n' "$min" "$v" | sort -V | head -1)" == "$min" ]]; then
            printf -v "$var" '%s' "--no-pgc"
        else
            printf -v "$var" '%s' ""
        fi
        export "${var?}"
    fi
    echo "${!var}"
}

# Canonical-named copies of the null-model region files (made in step 01)
include_bed() { [[ -n "$INCLUDE_BED" ]] && echo "$WORK/annot/include.bed" || true; }
exclude_bed() { [[ -n "$EXCLUDE_BED" ]] && echo "$WORK/annot/exclude.bed" || true; }

_shuffle_args() {
    local g; g=$(genome_file)
    local a="-g $g -noOverlapping -maxTries 10000"
    [[ -n "$INCLUDE_BED" ]] && a="$a -incl $(include_bed)"
    [[ -n "$EXCLUDE_BED" ]] && a="$a -excl $(exclude_bed)"
    [[ "$SAME_CHROM" == "1" ]] && a="$a -chrom"
    echo "$a"
}

# bedtools intersect -sorted on a STREAM. With -sorted, bedtools may stop reading -a as soon as it passes
# the last chromosome present in -b; the upstream writer is then killed by SIGPIPE, which pipefail treats as
# a failure even though the result is complete. Reading the whole stream into a file first avoids that.
# Usage: ... | isect_sorted -b features.bed [other intersect options]
isect_sorted() {
    mkdir -p "$WORK/tmp"
    local t rc=0; t=$(mktemp "$WORK/tmp/isect.XXXXXX")
    cat > "$t"
    bedtools intersect -a "$t" "$@" || rc=$?
    rm -f "$t"; return $rc
}

# Clip intervals to chromosome ends (bedtools shuffle -incl only constrains the start)
_clip() {
    awk 'BEGIN{FS=OFS="\t"} NR==FNR { sz[$1]=$2+0; next } { if ($3+0 > sz[$1]) $3 = sz[$1]; if ($2+0 < $3+0) print }' "$(genome_file)" -
}

# One null relocation of a BED4+ file (widths and names preserved) -> stdout BED4, sorted.
# uniform:    bedtools shuffle once
# gc_matched: GC_CANDIDATES relocations per peak; keep the one whose GC is closest to the real peak
shuffle_once() {
    local bed="$1" seed="$2" tmp
    local args; args=$(_shuffle_args)
    if [[ "$NULL_MODEL" == "uniform" ]]; then
        cut -f1-4 "$bed" | bedtools shuffle -i stdin $args -seed "$seed" | _clip | sort -k1,1 -k2,2n
        return
    fi
    tmp=$(mktemp -d)
    cut -f1-4 "$bed" > "$tmp/in.bed"
    to_fasta_names < "$tmp/in.bed" | bedtools nuc -fi "$GENOME_FASTA" -bed stdin | awk 'BEGIN{FS=OFS="\t"} NR>1{print $4,$6}' > "$tmp/obs.gc"
    for k in $(seq 1 "$GC_CANDIDATES"); do
        bedtools shuffle -i "$tmp/in.bed" $args -seed $(( seed * 1000 + k ))
    done > "$tmp/cand.bed"
    # candidates that run past a chromosome end are skipped by bedtools with a warning per window; hide only that warning
    to_fasta_names < "$tmp/cand.bed" | bedtools nuc -fi "$GENOME_FASTA" -bed stdin 2> >(grep -v "beyond the length of" >&2 || true) \
        | awk 'BEGIN{FS=OFS="\t"}
            NR==FNR { gc[$1]=$2; next }
            FNR==1 { next }
            { d = $6 - gc[$4]; if (d < 0) d = -d
              if (!($4 in best) || d < best[$4]) { best[$4]=d; line[$4]=$1 OFS $2 OFS $3 OFS $4 } }
            END { for (k in line) print line[k] }' "$tmp/obs.gc" - \
        | from_any_names | _clip | sort -k1,1 -k2,2n
    rm -rf "$tmp"
}

# Count peaks per feature label under the configured overlap rule.
# in:  BED4 peaks (name unique), class BED4 (chrom start end label), sorted
# out: label \t n_peaks   (plus a row "ANY")
count_labels() {
    local bed="$1" classbed="$2"
    if [[ "$OVERLAP_RULE" == "midpoint" ]]; then
        awk 'BEGIN{FS=OFS="\t"} { m=int(($2+$3)/2); print $1, m, m+1, $4 }' "$bed" | sort -k1,1 -k2,2n \
            | isect_sorted -b "$classbed" -wa -wb -sorted -g "$(genome_file)" \
            | awk 'BEGIN{FS=OFS="\t"}
                { k=$4 SUBSEP $8; if (!(k in s)) { s[k]=1; n[$8]++ } if (!($4 in a)) { a[$4]=1; any++ } }
                END { for (l in n) print l, n[l]; print "ANY", any+0 }'
    else
        cut -f1-4 "$bed" | isect_sorted -b "$classbed" -wo -sorted -g "$(genome_file)" \
            | awk -v f="$OVERLAP_FRACTION" 'BEGIN{FS=OFS="\t"}
                { len[$4]=$3-$2; ov[$4 SUBSEP $8]+=$9; tot[$4]+=$9 }
                END {
                    for (k in ov) { split(k, p, SUBSEP); if (ov[k] >= f*len[p[1]]) n[p[2]]++ }
                    for (q in tot) if (tot[q] >= f*len[q]) any++
                    for (l in n) print l, n[l]; print "ANY", any+0 }'
    fi
}

# Per-peak labels (observed peaks only): peak_id \t chrom \t label
label_peaks() {
    local bed="$1" classbed="$2"
    if [[ "$OVERLAP_RULE" == "midpoint" ]]; then
        awk 'BEGIN{FS=OFS="\t"} { m=int(($2+$3)/2); print $1, m, m+1, $4 }' "$bed" | sort -k1,1 -k2,2n \
            | isect_sorted -b "$classbed" -wa -wb -sorted -g "$(genome_file)" \
            | awk 'BEGIN{FS=OFS="\t"} { k=$4 SUBSEP $8; if (!(k in s)) { s[k]=1; print $4, $1, $8 } }'
    else
        cut -f1-4 "$bed" | isect_sorted -b "$classbed" -wo -sorted -g "$(genome_file)" \
            | awk -v f="$OVERLAP_FRACTION" 'BEGIN{FS=OFS="\t"}
                { len[$4]=$3-$2; chr[$4]=$1; ov[$4 SUBSEP $8]+=$9 }
                END { for (k in ov) { split(k, p, SUBSEP); if (ov[k] >= f*len[p[1]]) print p[1], chr[p[1]], p[2] } }'
    fi
}

export SCAN_VERSION
export -f pgc_flag isect_sorted _clip log sample_col sample_get peaks_to_bed5 peaks_info genome_file rename_chroms to_fasta_names from_any_names getfasta_canonical include_bed exclude_bed _shuffle_args shuffle_once count_labels label_peaks
