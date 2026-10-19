#!/usr/bin/env python3
"""
Per-base mappability from single-copy k-mer mappability.

Input is a k-mer mappability bedGraph (GenMap or Umap style): the value at position i describes the
k-mer starting at i; a value of 1 means the k-mer occurs exactly once in the genome (both strands).
Per-base mappability is then the fraction of the k k-mers overlapping a base that are single-copy:

    M(i) = #{ j in [i-k+1, i] : k-mer starting at j is single-copy } / k

(the Umap "multi-read" definition; with k = insert size, M(i) is the probability that a fragment of
that length covering base i can be placed uniquely). Computed one chromosome at a time.

Writes (to --outdir):
  censat_mappability_hist.tsv   class, mappable_kmers (0..k), n_bases   for each CenSat class,
                                'outside_censat' and 'genome'
  peak_mappability.tsv.gz       sample_id, peak_id, chrom, start, end, mean_mappability,
                                frac_unique_bases, min_mappability     for every peak of every sample
  mappability_k<k>.bedgraph.gz  per-base mappability track (optional, --track)
"""
import argparse
import glob
import gzip
import os
import sys
from collections import defaultdict

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bedgraph", required=True)
    ap.add_argument("--alias", required=True)
    ap.add_argument("--sizes", required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--censat", default="")
    ap.add_argument("--peaks-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--track", action="store_true")
    a = ap.parse_args()
    k = a.k
    dt = np.uint8 if k <= 255 else np.uint16

    alias = dict(l.rstrip("\n").split("\t")[:2] for l in open(a.alias) if l.strip())
    sizes = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.sizes) if l.strip()}

    # single-copy k-mer start intervals per canonical chromosome
    uniq = defaultdict(list)
    opener = gzip.open if a.bedgraph.endswith(".gz") else open
    n_rec = 0
    with opener(a.bedgraph, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("track", "#", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            c = f[0].split()[0]
            c = alias.get(c, c)
            if c not in sizes:
                continue
            try:
                v = float(f[3])
            except ValueError:
                continue
            if v >= 0.999:   # value 1 = k-mer occurs once
                uniq[c].append((int(f[1]), int(f[2])))
            n_rec += 1
    if not uniq:
        sys.exit("[qc_mappability] no single-copy k-mers found; check the bedGraph and chromosome names")
    print(f"[qc_mappability] read {n_rec} bedGraph records; single-copy k-mers on {len(uniq)} chromosomes", file=sys.stderr)

    censat = defaultdict(list)
    if a.censat and os.path.exists(a.censat):
        for l in open(a.censat):
            c, s, e, lab = l.rstrip("\n").split("\t")[:4]
            censat[c].append((int(s), int(e), lab))
    peaks = defaultdict(list)
    for f in sorted(glob.glob(os.path.join(a.peaks_dir, "*.bed5"))):
        for l in open(f):
            c, s, e, name = l.rstrip("\n").split("\t")[:4]
            peaks[c].append((int(s), int(e), name))

    hist = defaultdict(lambda: np.zeros(k + 1, dtype=np.int64))
    os.makedirs(a.outdir, exist_ok=True)
    pk_out = gzip.open(os.path.join(a.outdir, "peak_mappability.tsv.gz.tmp"), "wt")
    pk_out.write("sample_id\tpeak_id\tchrom\tstart\tend\tmean_mappability\tfrac_unique_bases\tmin_mappability\n")
    tr_out = gzip.open(os.path.join(a.outdir, f"mappability_k{k}.bedgraph.gz.tmp"), "wt") if a.track else None

    for c in sorted(sizes):
        L = sizes[c]
        u = np.zeros(L, dtype=np.uint8)
        for s, e in uniq.get(c, ()):
            u[max(0, s):min(L, e)] = 1
        cs = np.zeros(L + 1, dtype=np.int32)
        np.cumsum(u, out=cs[1:])
        del u
        mc = cs[1:].copy()
        if L > k - 1:
            mc[k - 1:] -= cs[:L - k + 1]
        del cs
        mc = mc.astype(dt)          # number of single-copy k-mers covering each base (0..k)

        hist["genome"] += np.bincount(mc, minlength=k + 1)
        mask = np.zeros(L, dtype=bool)
        for s, e, lab in censat.get(c, ()):
            s, e = max(0, s), min(L, e)
            if e > s:
                hist[lab] += np.bincount(mc[s:e], minlength=k + 1)
                mask[s:e] = True
        hist["outside_censat"] += np.bincount(mc[~mask], minlength=k + 1)
        del mask

        rows = []
        for s, e, name in peaks.get(c, ()):
            s, e = max(0, s), min(L, e)
            if e <= s:
                continue
            seg = mc[s:e]
            rows.append(f"{name.split('|', 1)[0]}\t{name}\t{c}\t{s}\t{e}\t{seg.mean() / k:.4f}\t"
                        f"{np.count_nonzero(seg == k) / seg.size:.4f}\t{seg.min() / k:.4f}\n")
        pk_out.write("".join(rows))

        if tr_out is not None:
            ch = np.flatnonzero(mc[1:] != mc[:-1]) + 1
            st = np.concatenate(([0], ch)); en = np.concatenate((ch, [L]))
            vals = mc[st].astype(np.float64) / k
            for i in range(0, st.size, 500000):
                tr_out.write("".join(f"{c}\t{x}\t{y}\t{v:.4g}\n" for x, y, v in
                                     zip(st[i:i + 500000].tolist(), en[i:i + 500000].tolist(), vals[i:i + 500000].tolist())))
        print(f"[qc_mappability] {c}: done", file=sys.stderr)
        del mc

    pk_out.close()
    os.replace(os.path.join(a.outdir, "peak_mappability.tsv.gz.tmp"), os.path.join(a.outdir, "peak_mappability.tsv.gz"))
    if tr_out is not None:
        tr_out.close()
        os.replace(os.path.join(a.outdir, f"mappability_k{k}.bedgraph.gz.tmp"), os.path.join(a.outdir, f"mappability_k{k}.bedgraph.gz"))
    with open(os.path.join(a.outdir, "censat_mappability_hist.tsv"), "w") as fh:
        fh.write("class\tmappable_kmers\tn_bases\n")
        for lab in sorted(hist):
            for i, n in enumerate(hist[lab].tolist()):
                if n:
                    fh.write(f"{lab}\t{i}\t{n}\n")


if __name__ == "__main__":
    main()
