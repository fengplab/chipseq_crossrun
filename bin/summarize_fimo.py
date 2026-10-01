#!/usr/bin/env python3
"""
Summarise a FIMO --text stream (motif-major order, as FIMO writes it) for one sample and one set.

Outputs (prefix = --out):
  <out>.summary.tsv   sample_id set family_id n_seqs n_seqs_hit n_sites max_sites_per_seq
  <out>.spacing.tsv   sample_id family_id spacing n      (peaks set only; same-strand adjacent sites)
  <out>.sites.bed.gz  chrom start end family_id peak_id strand   (peaks set only; genomic coordinates)

Sites of one motif with identical start/stop on opposite strands (palindromes) are counted once.

Usage: fimo --text ... | summarize_fimo.py --sample S --set peaks|background --scan-bed S.scan.bed \
          --families family_rep_features.tsv --spacing-max 400 --out prefix
"""
import argparse
import gzip
import sys
from collections import Counter, defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True)
ap.add_argument("--set", required=True, choices=["peaks", "background"])
ap.add_argument("--scan-bed", required=True)
ap.add_argument("--families", required=True)
ap.add_argument("--spacing-max", type=int, default=400)
ap.add_argument("--out", required=True)
a = ap.parse_args()

coords = {}
with open(a.scan_bed) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        coords[f[3]] = (f[0], int(f[1]))
n_seqs = len(coords)
families = [l.split("\t")[0] for l in open(a.families).read().splitlines()[1:] if l.strip()]

peaks = a.set == "peaks"
summary = {}
spacing = defaultdict(Counter)
bed = gzip.open(a.out + ".sites.bed.gz", "wt") if peaks else None


def flush(fam, per_seq):
    n_hit = len(per_seq)
    n_sites = 0
    mx = 0
    for seq, sites in per_seq.items():
        uniq = {}
        for st, sp, strand in sites:
            uniq.setdefault((st, sp), strand)
        k = len(uniq)
        n_sites += k
        mx = max(mx, k)
        if not peaks:
            continue
        chrom, off = coords.get(seq, (None, 0))
        by_strand = defaultdict(list)
        for (st, sp), strand in uniq.items():
            by_strand[strand].append(st)
            if chrom is not None:
                bed.write(f"{chrom}\t{off + st - 1}\t{off + sp}\t{fam}\t{seq}\t{strand}\n")
        for pos in by_strand.values():
            pos.sort()
            for x, y in zip(pos, pos[1:]):
                d = y - x
                if 0 < d <= a.spacing_max:
                    spacing[fam][d] += 1
    prev = summary.get(fam, (0, 0, 0))
    summary[fam] = (prev[0] + n_hit, prev[1] + n_sites, max(prev[2], mx))


idx = None
cur, per_seq = None, defaultdict(list)
for line in sys.stdin:
    if not line.strip() or line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if idx is None:
        idx = {h: i for i, h in enumerate(f)}
        need = ("motif_id", "sequence_name", "start", "stop", "strand")
        if not all(k in idx for k in need):
            sys.exit(f"unexpected FIMO header: {f}")
        continue
    fam = f[idx["motif_id"]]
    if fam != cur:
        if cur is not None:
            flush(cur, per_seq)
        cur, per_seq = fam, defaultdict(list)
    seq = f[idx["sequence_name"]].split("::")[0]
    per_seq[seq].append((int(f[idx["start"]]), int(f[idx["stop"]]), f[idx["strand"]]))
if cur is not None:
    flush(cur, per_seq)
if bed:
    bed.close()

with open(a.out + ".summary.tsv", "w") as fh:
    fh.write("sample_id\tset\tfamily_id\tn_seqs\tn_seqs_hit\tn_sites\tmax_sites_per_seq\n")
    for fam in families:
        h, s, m = summary.get(fam, (0, 0, 0))
        fh.write(f"{a.sample}\t{a.set}\t{fam}\t{n_seqs}\t{h}\t{s}\t{m}\n")
if peaks:
    with open(a.out + ".spacing.tsv", "w") as fh:
        fh.write("sample_id\tfamily_id\tspacing\tn\n")
        for fam, c in spacing.items():
            for d in sorted(c):
                fh.write(f"{a.sample}\t{fam}\t{d}\t{c[d]}\n")
