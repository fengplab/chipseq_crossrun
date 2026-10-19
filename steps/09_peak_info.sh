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

# Reads per contig and MAPQ for every IP and control BAM (Figure P5, table00c, qc07, pies). Cached per BAM.
genome_file >/dev/null
if command -v samtools >/dev/null; then
    { sample_col bam_file | awk 'BEGIN{OFS="\t"} {print $1, "IP", $2}'
      sample_col control_bam | awk 'BEGIN{OFS="\t"} {print $1, "control", $2}'; } > "$WORK/contigs.list.tmp" 2>/dev/null || true
    if [[ -s "$WORK/contigs.list.tmp" ]]; then
        mkdir -p "$WORK/contigs"
        cut -f3 "$WORK/contigs.list.tmp" | sort -u | xargs -P "$THREADS" -I{} bash -c 'contig_counts "$1" > /dev/null' _ {}
        { printf "sample_id\tlibrary\tbam\tcounts_file\tpaired\n"
          while IFS=$'\t' read -r s lib b; do
              [[ -s "$b" ]] || continue
              pe=$([[ $(samtools view -f 1 "$b" 2>/dev/null | head -n 1000 | wc -l) -gt 0 ]] && echo 1 || echo 0)
              printf "%s\t%s\t%s\t%s\t%s\n" "$s" "$lib" "$b" "$(contig_counts "$b")" "$pe"
          done < "$WORK/contigs.list.tmp"; } > "$WORK/contigs/bams.tsv"
        log "reads per contig counted for $(cut -f3 "$WORK/contigs.list.tmp" | sort -u | wc -l) BAM files"
    fi
    rm -f "$WORK/contigs.list.tmp"
    # The pipeline's own counts on chrM and the viral genomes (contig_reads.tsv, before host filtering) and its
    # duplicate-marked viral alignments. The viral BAMs give the MAPQ table (qc07) and, with contigs.viral_mapq > 0,
    # the viral counts at that cutoff. work/contigs/tables.tsv: sample_id library contig_table viral_bam viral_counts_file
    : > "$WORK/contigs.tab.tmp"; : > "$WORK/contigs.vbam.tmp"
    if head -n 1 "$SAMPLES" | tr '\t' '\n' | grep -qx contig_table; then     # samples.tsv from step 0 of this version
        { sample_col contig_table | awk 'BEGIN{OFS="\t"} {print $1, "IP", $2}'
          sample_col control_contig_table | awk 'BEGIN{OFS="\t"} {print $1, "control", $2}'; } > "$WORK/contigs.tab.tmp"
        { sample_col viral_bam | awk 'BEGIN{OFS="\t"} {print $1, "IP", $2}'
          sample_col control_viral_bam | awk 'BEGIN{OFS="\t"} {print $1, "control", $2}'; } > "$WORK/contigs.vbam.tmp"
    else
        log "NOTE: work/samples.tsv predates the contig_reads.tsv support; rerun step 0 to use the pipeline's chrM/viral counts"
    fi
    if [[ -s "$WORK/contigs.tab.tmp" ]]; then
        mkdir -p "$WORK/contigs"
        cut -f3 "$WORK/contigs.vbam.tmp" | sort -u | while read -r b; do
            [[ -s "$b" ]] && echo "$b"
        done | xargs -r -P "$THREADS" -I{} bash -c 'contig_counts "$1" > /dev/null' _ {}
        { printf "sample_id\tlibrary\tcontig_table\tviral_bam\tviral_counts_file\n"
          while IFS=$'\t' read -r s lib t; do
              [[ -s "$t" ]] || continue
              vb=$(awk -F'\t' -v s="$s" -v l="$lib" '$1 == s && $2 == l {print $3}' "$WORK/contigs.vbam.tmp")
              vc=""; [[ -n "$vb" && -s "$vb" ]] && vc=$(contig_counts "$vb")
              printf "%s\t%s\t%s\t%s\t%s\n" "$s" "$lib" "$t" "$vb" "$vc"
          done < "$WORK/contigs.tab.tmp"; } > "$WORK/contigs/tables.tsv"
        log "pipeline contig read tables: $(($(wc -l < "$WORK/contigs/tables.tsv") - 1)) libraries ($(cut -f3 "$WORK/contigs.vbam.tmp" | sort -u | grep -c . || true) viral BAMs)"
    else
        rm -f "$WORK/contigs/tables.tsv"
        log "no pipeline contig_reads.tsv tables found: chrM and viral reads are counted from the filtered BAMs"
    fi
    rm -f "$WORK/contigs.tab.tmp" "$WORK/contigs.vbam.tmp"
else
    log "WARNING: samtools not found; reads per contig (Figure P5, pies) are skipped"
fi
