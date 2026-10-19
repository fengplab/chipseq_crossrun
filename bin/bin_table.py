#!/usr/bin/env python3
"""Genome-wide bin table for the binding models (step 16) and the bin-size tuning (tune/).

Subcommands (all coordinates 0-based half-open, BED convention; bins tile each host nuclear chromosome from 0):

  annotate  per bin: repeat class, CenSat array ID, distance to the array edge, repeat fraction, mappability, GC, N fraction
  sites     per bin and motif family: number of FIMO sites whose midpoint lies in the bin
  signal    per bin and library: mean per-base fragment coverage (all reads; unique reads; two read-name halves) and
            the coverage-weighted mean fragment length

Storage: <out>/<chrom>/<column>.<dtype> raw little-endian arrays (one value per bin), read in R with readBin
(R/bin_io.R). Every directory carries a manifest.tsv (chrom, n_bins, bin_size) and the bin size is part of every path,
so tables made at different bin sizes never mix.

Only numpy and pandas are needed. Reads are streamed from samtools; nothing is held per read.
"""
import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd


def log(*a):
    print("[bin_table]", *a, file=sys.stderr, flush=True)


# ------------------------------------------------------------------------------------------------ helpers
def read_sizes(path, chroms=None):
    s = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"], dtype={"chrom": str})
    if chroms:
        missing = [c for c in chroms if c not in set(s.chrom)]
        if missing:
            sys.exit(f"[bin_table] chromosomes not in {path}: {', '.join(missing)}")
        s = s[s.chrom.isin(chroms)]
    return dict(zip(s.chrom, s.length.astype(np.int64)))


def read_alias(path):
    a = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], names=["from", "to"], dtype=str)
    return dict(zip(a["from"], a["to"]))


def n_bins(length, w):
    return int((length + w - 1) // w)


def bin_lengths(length, w):
    nb = n_bins(length, w)
    bl = np.full(nb, w, dtype=np.float64)
    bl[-1] = length - (nb - 1) * w
    return bl


def accumulate(starts, ends, weights, nb, w, length):
    """Sum over intervals of weight x (bp of the interval inside each bin). Intervals are clipped to [0, length).
    Exact: partial first/last bins get their overlap, bins fully inside get w (difference array)."""
    s = np.clip(np.asarray(starts, dtype=np.int64), 0, length)
    e = np.clip(np.asarray(ends, dtype=np.int64), 0, length)
    wt = np.asarray(weights, dtype=np.float64)
    ok = e > s
    s, e, wt = s[ok], e[ok], wt[ok]
    out = np.zeros(nb, dtype=np.float64)
    if not len(s):
        return out
    b0 = s // w
    b1 = (e - 1) // w
    same = b0 == b1
    out += np.bincount(b0[same], weights=wt[same] * (e[same] - s[same]), minlength=nb)
    d = ~same
    if d.any():
        out += np.bincount(b0[d], weights=wt[d] * ((b0[d] + 1) * w - s[d]), minlength=nb)
        out += np.bincount(b1[d], weights=wt[d] * (e[d] - b1[d] * w), minlength=nb)
        full = d & (b1 > b0 + 1)
        if full.any():
            diff = np.bincount(b0[full] + 1, weights=wt[full] * w, minlength=nb + 1)[: nb + 1]
            diff -= np.bincount(b1[full], weights=wt[full] * w, minlength=nb + 1)[: nb + 1]
            out += np.cumsum(diff)[:nb]
    return out


def write_col(d, name, arr, dtype):
    os.makedirs(d, exist_ok=True)
    arr = np.asarray(arr).astype(dtype)
    tmp = os.path.join(d, f".{name}.{dtype}.{os.getpid()}")
    arr.astype(f"<{np.dtype(dtype).str[1:]}").tofile(tmp)
    os.replace(tmp, os.path.join(d, f"{name}.{dtype}"))


def write_manifest(out, sizes, w, extra=None):
    rows = [{"chrom": c, "length": L, "n_bins": n_bins(L, w), "bin_size": w} for c, L in sizes.items()]
    m = pd.DataFrame(rows)
    if extra:
        for k, v in extra.items():
            m[k] = v
    m.to_csv(os.path.join(out, "manifest.tsv"), sep="\t", index=False)


def read_bed(path, cols=4):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=["chrom", "start", "end", "label"][:cols])
    return pd.read_csv(path, sep="\t", header=None, usecols=list(range(cols)), names=["chrom", "start", "end", "label"][:cols],
                       dtype={"chrom": str, "label": str}, comment="#")


def merge_intervals(s, e):
    """Union of intervals (sorted output)."""
    if not len(s):
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    o = np.lexsort((e, s))
    s, e = s[o], e[o]
    run_end = np.maximum.accumulate(e)
    new = np.empty(len(s), dtype=bool)
    new[0] = True
    new[1:] = s[1:] > run_end[:-1]
    gid = np.cumsum(new) - 1
    ms = s[new]
    me = np.zeros(len(ms), dtype=np.int64)
    np.maximum.at(me, gid, e)
    return ms, me


def point_lookup(pos, starts, ends):
    """Index of the interval containing each position (intervals sorted by start, non-overlapping), else -1."""
    if not len(starts):
        return np.full(len(pos), -1, dtype=np.int64)
    i = np.searchsorted(starts, pos, side="right") - 1
    ok = (i >= 0) & (pos < ends[np.clip(i, 0, None)])
    return np.where(ok, i, -1)


# ------------------------------------------------------------------------------------------------ annotate
def cmd_annotate(a):
    sizes = read_sizes(a.sizes, a.chroms)
    w = a.bin
    out = os.path.join(a.out, f"w{w}")
    os.makedirs(out, exist_ok=True)
    cen = read_bed(a.censat)
    rm = read_bed(a.rmsk)
    rm["label"] = rm["label"].str.split("/").str[0]                     # RepeatMasker class (SINE/Alu -> SINE)
    classes = sorted(cen["label"].dropna().unique().tolist()) + \
        [c for c in sorted(rm["label"].dropna().unique().tolist()) if c not in set(cen["label"])]
    code = {c: i + 1 for i, c in enumerate(classes)}                    # 0 = no repeat
    pd.DataFrame({"code": [0] + [code[c] for c in classes], "class": ["none"] + classes,
                  "source": ["none"] + ["censat" if c in set(cen["label"]) else "repeatmasker" for c in classes]}) \
        .to_csv(os.path.join(out, "classes.tsv"), sep="\t", index=False)
    fai = None
    if a.fasta:
        fai = pd.read_csv(a.fasta + ".fai", sep="\t", header=None, usecols=[0, 1, 2, 3, 4],
                          names=["name", "length", "offset", "linebases", "linewidth"], dtype={"name": str}).set_index("name")
        fnames = dict(pd.read_csv(a.fasta_names, sep="\t", header=None, names=["canon", "fasta"], dtype=str).values) \
            if a.fasta_names else {}
    mp = None
    if a.mappability:
        log("reading mappability", a.mappability)
        mp = pd.read_csv(a.mappability, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "v"],
                         dtype=str)
        mp = mp[mp.chrom.isin(sizes)]                                     # also drops 'track'/'browser' lines
        for k in ("start", "end"):
            mp[k] = mp[k].astype(np.int64)
        mp["v"] = mp["v"].astype(np.float64)
    arr_rows = []
    next_array = 0
    for c, L in sizes.items():
        nb = n_bins(L, w)
        bl = bin_lengths(L, w)
        mid = np.arange(nb, dtype=np.int64) * w + (bl // 2).astype(np.int64)
        d = os.path.join(out, c)
        # repeat class by bin midpoint: CenSat first, then RepeatMasker
        cls = np.zeros(nb, dtype=np.int16)
        ce = cen[cen.chrom == c].sort_values("start")
        arr_id = np.full(nb, -1, dtype=np.int32)
        edge = np.full(nb, -1, dtype=np.int32)
        if len(ce):
            # arrays: CenSat records of the same class that touch or overlap are one array
            cs, cend, cl = ce.start.values.astype(np.int64), ce.end.values.astype(np.int64), ce.label.values
            grp = np.zeros(len(cs), dtype=np.int64)
            g = 0
            for i in range(1, len(cs)):
                if not (cl[i] == cl[i - 1] and cs[i] <= cend[i - 1]):
                    g += 1
                grp[i] = g
            A = pd.DataFrame({"s": cs, "e": cend, "cl": cl, "g": grp}).groupby("g").agg(s=("s", "min"), e=("e", "max"), cl=("cl", "first"))
            As, Ae = A.s.values, A.e.values
            k = point_lookup(mid, As, Ae)
            hit = k >= 0
            cls[hit] = [code[x] for x in A.cl.values[k[hit]]]
            arr_id[hit] = (k[hit] + next_array).astype(np.int32)
            edge[hit] = np.minimum(mid[hit] - As[k[hit]], Ae[k[hit]] - mid[hit]).astype(np.int32)
            for j in range(len(As)):
                arr_rows.append((next_array + j, c, int(As[j]), int(Ae[j]), A.cl.values[j]))
            next_array += len(As)
        r = rm[rm.chrom == c]
        if len(r):
            rs, re_ = merge_intervals(r.start.values.astype(np.int64), r.end.values.astype(np.int64))
            # class of the RepeatMasker feature holding the midpoint (largest-start feature at or before it)
            rr = r.sort_values("start")
            i = np.searchsorted(rr.start.values, mid, side="right") - 1
            ok = (i >= 0) & (mid < rr.end.values[np.clip(i, 0, None)]) & (cls == 0)
            cls[ok] = [code[x] for x in rr.label.values[i[ok]]]
        else:
            rs, re_ = np.array([], dtype=np.int64), np.array([], dtype=np.int64)
        # repeat fraction: union of CenSat and RepeatMasker intervals over the bin
        us, ue = merge_intervals(np.concatenate([rs, ce.start.values.astype(np.int64)]) if len(ce) else rs,
                                 np.concatenate([re_, ce.end.values.astype(np.int64)]) if len(ce) else re_)
        repfrac = accumulate(us, ue, np.ones(len(us)), nb, w, L) / bl
        write_col(d, "class", cls, "int16")
        write_col(d, "array_id", arr_id, "int32")
        write_col(d, "edge_dist", edge, "int32")
        write_col(d, "repeat_frac", repfrac, "float32")
        # mappability: mean over the bases the track covers (NaN where it covers none)
        if mp is not None:
            m = mp[mp.chrom == c]
            num = accumulate(m.start.values, m.end.values, m.v.values, nb, w, L)
            den = accumulate(m.start.values, m.end.values, np.ones(len(m)), nb, w, L)
            with np.errstate(invalid="ignore", divide="ignore"):
                write_col(d, "mappability", np.where(den > 0, num / den, np.nan), "float32")
        # GC and N fraction from the FASTA
        if fai is not None:
            fn = fnames.get(c, c)
            if fn not in fai.index:
                sys.exit(f"[bin_table] {c} ({fn}) not in {a.fasta}.fai")
            f = fai.loc[fn]
            with open(a.fasta, "rb") as fh:
                fh.seek(int(f.offset))
                nlines = (int(f.length) + int(f.linebases) - 1) // int(f.linebases)
                raw = fh.read(int(f.length) + nlines * (int(f.linewidth) - int(f.linebases)))
            seq = np.frombuffer(raw, dtype=np.uint8)
            seq = seq[(seq != 10) & (seq != 13)][:L]
            if len(seq) != L:
                sys.exit(f"[bin_table] {c}: read {len(seq)} bases from the FASTA, expected {L}")
            up = seq & 0xDF                                               # upper case
            gc = (up == ord("G")) | (up == ord("C"))
            nn = up == ord("N")
            idx = np.arange(0, L, w)
            gcs = np.add.reduceat(gc.astype(np.int64), idx)
            ns = np.add.reduceat(nn.astype(np.int64), idx)
            with np.errstate(invalid="ignore", divide="ignore"):
                write_col(d, "gc", np.where(bl - ns > 0, gcs / (bl - ns), np.nan), "float32")
            write_col(d, "n_frac", ns / bl, "float32")
        log(f"annotated {c}: {nb:,} bins of {w} bp")
    pd.DataFrame(arr_rows, columns=["array_id", "chrom", "start", "end", "class"]).to_csv(
        os.path.join(out, "arrays.tsv"), sep="\t", index=False)
    write_manifest(out, sizes, w)


# ------------------------------------------------------------------------------------------------ sites
def cmd_sites(a):
    """FIMO --text output (any number of families) -> per family and chromosome: sites per bin (sparse TSV).
    Sites reported on both strands at the same position (palindromes) count once."""
    sizes = read_sizes(a.sizes, a.chroms)
    alias = read_alias(a.alias) if a.alias else {}
    w = a.bin
    out = os.path.join(a.out, f"w{w}")
    os.makedirs(out, exist_ok=True)
    parts = []
    for f in a.fimo:
        op = gzip.open if f.endswith(".gz") else open
        with op(f, "rt") as fh:
            first = fh.readline()
        df = pd.read_csv(f, sep="\t", comment="#", usecols=[0, 2, 3, 4], header=0 if first.startswith("motif_id") else None,
                         names=["family", "seq", "start", "stop"], dtype={"family": str, "seq": str})
        if not len(df):
            continue
        df["chrom"] = df.seq.map(lambda x: alias.get(x, x))
        df = df[df.chrom.isin(sizes)].drop_duplicates(["family", "chrom", "start", "stop"])
        # FIMO is 1-based closed [start, stop]; midpoint in 0-based coordinates = (start - 1 + stop) // 2
        df["bin"] = ((df.start.astype(np.int64) - 1 + df.stop.astype(np.int64)) // 2) // w
        parts.append(df.groupby(["family", "chrom", "bin"]).size().rename("n").reset_index())
    res = pd.concat(parts) if parts else pd.DataFrame(columns=["family", "chrom", "bin", "n"])
    res = res.groupby(["family", "chrom", "bin"], as_index=False)["n"].sum()
    res.to_csv(os.path.join(out, a.name), sep="\t", index=False, compression="gzip")
    log(f"sites: {len(res):,} non-empty bin x family cells, {res.n.sum() if len(res) else 0:,} sites -> {a.name}")


# ------------------------------------------------------------------------------------------------ signal
def cmd_signal(a):
    """One pass over a BAM: fragments (proper pairs, read 1 only; or single reads) on host nuclear chromosomes.
    Columns written per chromosome (mean per-base fragment coverage, i.e. sum of overlap bp / bin length):
      all, uniq (MAPQ >= --unique-mapq), h1, h2 (halves by a hash of the read name), fl (sum of fragment length x
      overlap / bin length; divide by 'all' for the coverage-weighted mean fragment length)."""
    sizes = read_sizes(a.sizes, a.chroms)
    alias = read_alias(a.alias) if a.alias else {}
    w = a.bin
    out = os.path.join(a.out, f"w{w}", a.name)
    os.makedirs(out, exist_ok=True)
    acc = {c: {k: np.zeros(n_bins(L, w)) for k in ("all", "uniq", "h1", "h2", "fl")} for c, L in sizes.items()}
    tot = {k: 0 for k in ("fragments", "unique", "h1", "h2", "single_end")}
    # chromosomes as named in the BAM
    hdr = subprocess.run(["samtools", "view", "-H", a.bam], capture_output=True, text=True, check=True).stdout
    bam_names = [l.split("\tSN:")[1].split("\t")[0] for l in hdr.splitlines() if l.startswith("@SQ") and "\tSN:" in l]
    regions = [n for n in bam_names if alias.get(n, n) in sizes]
    if not regions:
        sys.exit(f"[bin_table] {a.bam}: none of the requested chromosomes is in the BAM header")
    paired = subprocess.run(f"samtools view -f 1 '{a.bam}' {regions[0]} | head -n 1000 | wc -l", shell=True,
                            capture_output=True, text=True, check=True).stdout.strip() != "0"
    # -F 0xF0C: unmapped, mate unmapped, secondary, QC fail, duplicate, supplementary; paired: proper pair, read 1 only
    flt = ["-f", "66", "-F", "3852"] if paired else ["-F", "3844"]
    cmd = ["samtools", "view", "-@", str(a.threads)] + flt + [a.bam] + regions
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    cols = [0, 2, 3, 4, 7, 8, 9] if paired else [0, 2, 3, 4, 9]
    names = ["qname", "rname", "pos", "mapq", "pnext", "tlen", "seq"] if paired else ["qname", "rname", "pos", "mapq", "seq"]
    reader = pd.read_csv(p.stdout, sep="\t", header=None, usecols=cols, names=names, dtype={"qname": str, "rname": str},
                         chunksize=2_000_000, quoting=3, engine="c")
    for ch in reader:
        ch["chrom"] = ch.rname.map(lambda x: alias.get(x, x))
        ch = ch[ch.chrom.isin(sizes)]
        if paired:
            fl = ch.tlen.abs().astype(np.int64)
            ch = ch[(fl > 0) & (fl <= a.insert_max)]
            fl = ch.tlen.abs().astype(np.int64)
            st = np.minimum(ch.pos.values, ch.pnext.values).astype(np.int64) - 1   # SAM 1-based -> BED 0-based
        else:
            fl = ch.seq.str.len().astype(np.int64)
            st = ch.pos.values.astype(np.int64) - 1
            tot["single_end"] += len(ch)
        en = st + fl.values
        h = pd.util.hash_array(ch.qname.values.astype(object), hash_key="crossrunhalves16") % 2 == 0
        uq = ch.mapq.values >= a.unique_mapq
        tot["fragments"] += len(ch); tot["unique"] += int(uq.sum()); tot["h1"] += int(h.sum()); tot["h2"] += int((~h).sum())
        for c, gi in ch.groupby("chrom").indices.items():
            L = sizes[c]; nb = n_bins(L, w); A = acc[c]
            s, e, f = st[gi], en[gi], fl.values[gi].astype(np.float64)
            one = np.ones(len(gi))
            A["all"] += accumulate(s, e, one, nb, w, L)
            A["fl"] += accumulate(s, e, f, nb, w, L)
            if uq[gi].any():
                A["uniq"] += accumulate(s[uq[gi]], e[uq[gi]], one[uq[gi]], nb, w, L)
            hh = h[gi]
            A["h1"] += accumulate(s[hh], e[hh], one[hh], nb, w, L)
            A["h2"] += accumulate(s[~hh], e[~hh], one[~hh], nb, w, L)
    if p.wait() != 0:
        sys.exit(f"[bin_table] samtools failed on {a.bam}")
    for c, L in sizes.items():
        bl = bin_lengths(L, w)
        for k, v in acc[c].items():
            write_col(os.path.join(out, c), k, v / bl, "float32")
    with open(os.path.join(out, "library.json"), "w") as fh:
        json.dump({"bam": os.path.abspath(a.bam), "paired": paired, "bin_size": w, "chroms": list(sizes), **tot}, fh, indent=1)
    write_manifest(out, sizes, w)
    log(f"signal {a.name}: {tot['fragments']:,} fragments ({tot['unique']:,} unique) on {len(sizes)} chromosomes, {w} bp bins")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    for name in ("annotate", "sites", "signal"):
        s = sp.add_parser(name)
        s.add_argument("--sizes", required=True, help="host nuclear chromosome sizes (canonical names)")
        s.add_argument("--bin", type=int, required=True)
        s.add_argument("--chroms", nargs="*", default=None, help="restrict to these chromosomes")
        s.add_argument("--out", required=True)
        s.add_argument("--alias", help="two columns: name in the file -> canonical name")
    s = sp.choices["annotate"]
    s.add_argument("--censat"); s.add_argument("--rmsk"); s.add_argument("--mappability")
    s.add_argument("--fasta"); s.add_argument("--fasta-names", help="canonical name -> FASTA name (work/annot/fasta_names.tsv)")
    s = sp.choices["sites"]
    s.add_argument("--fimo", nargs="+", required=True); s.add_argument("--name", default="sites.tsv.gz")
    s = sp.choices["signal"]
    s.add_argument("--bam", required=True); s.add_argument("--name", required=True)
    s.add_argument("--unique-mapq", type=int, default=255); s.add_argument("--insert-max", type=int, default=2000)
    s.add_argument("--threads", type=int, default=1)
    a = ap.parse_args()
    if a.bin < 1:
        sys.exit("[bin_table] --bin must be >= 1")
    {"annotate": cmd_annotate, "sites": cmd_sites, "signal": cmd_signal}[a.cmd](a)


if __name__ == "__main__":
    main()
