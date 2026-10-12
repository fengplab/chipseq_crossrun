#!/usr/bin/env python3
"""
Per-window features for the regression models (step 16).

  select    choose the peaks of one sample that enter the models
            - all peaks passing the q-value cutoff            -> <out>.peaks.bed  (chrom start end window_id neglog10_q)
            - a random subset of at most --max-peaks of them  -> <out>.model.bed  (chrom start end window_id)
              (only these get matched background windows and enter the peak-level mixed model)
  assemble  join everything computed for the windows of one sample
            - <out>.windows.tsv.gz  one row per window (peak or background)
            - <out>.sites.tsv.gz    window x motif family: number of sites and of sites in tandem arrays (non-zero only)
            - <out>.libsize.tsv     total reads of the IP and control libraries

A window is a peak (or summit window, as scanned in step 03) or a background window: the same window moved to
a GC-matched position elsewhere in the genome. Background window ids are "<peak window id>#bg<draw>".

Tandem arrays: sites of one family on one strand whose consecutive spacings are equal (+/- 1 bp) and at most
--spacing-max; a run of >= 3 such sites is an array, and every site in it counts as a tandem site.
"""
import argparse
import gzip
import random
import sys
from collections import defaultdict


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def cmd_select(a):
    q = {}
    n_noq = 0
    with opener(a.peaks_file) as fh:
        for line in fh:
            if line.startswith(("track", "browser", "#")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[3]:
                try:
                    q[f"{a.sample}|{f[3]}"] = float(f[8])
                except ValueError:
                    pass
    cut = None if a.max_q is None else a.max_q
    rows = []
    with open(a.scan_bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            v = q.get(f[3])
            if v is None:
                n_noq += 1
            # column 9 of MACS peak files is -log10(q); a negative value means "not available"
            qv = None if v is None or v < 0 else 10 ** (-v)
            if cut is not None and qv is not None and qv > cut:
                continue
            rows.append((f[0], f[1], f[2], f[3], "NA" if v is None or v < 0 else f"{v:.4g}"))
    if n_noq and cut is not None:
        print(f"[model_windows] {a.sample}: {n_noq} windows have no q-value in the peak file and were kept", file=sys.stderr)
    with open(a.out + ".peaks.bed", "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    sub = rows
    if a.max_peaks and len(rows) > a.max_peaks:
        rnd = random.Random(f"{a.seed}:{a.sample}")
        keep = set(rnd.sample(range(len(rows)), a.max_peaks))
        sub = [r for i, r in enumerate(rows) if i in keep]
    with open(a.out + ".model.bed", "w") as fh:
        for r in sub:
            fh.write("\t".join(r[:4]) + "\n")
    print(f"[model_windows] {a.sample}: {len(rows)} windows pass the q-value cutoff; {len(sub)} enter the peak-level model",
          file=sys.stderr)


def tandem_sites(starts, spacing_max):
    """Number of sites (sorted start positions on one strand) that belong to a run of >= 3 regularly spaced sites."""
    n = len(starts)
    if n < 3:
        return 0
    marked = [False] * n
    run_start = 0          # index of the first site of the current run
    prev_gap = None
    for i in range(1, n):
        gap = starts[i] - starts[i - 1]
        ok = 0 < gap <= spacing_max
        if ok and prev_gap is not None and abs(gap - prev_gap) <= 1:
            pass                                    # run continues
        else:
            if i - 1 - run_start + 1 >= 3:
                for j in range(run_start, i):
                    marked[j] = True
            run_start = i - 1 if ok else i
        prev_gap = gap if ok else None
    if n - run_start >= 3:
        for j in range(run_start, n):
            marked[j] = True
    return sum(marked)


def best_overlap(path):
    """bedtools intersect -wo of windows (BED4) with a class BED4 -> window_id -> (label, bp) of the largest overlap.
    Overlapping features of the same class (common in RepeatMasker) are merged, so bp never exceeds the window width."""
    iv = defaultdict(lambda: defaultdict(list))
    if path:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                s0, e0 = max(int(f[1]), int(f[5])), min(int(f[2]), int(f[6]))
                if e0 > s0:
                    iv[f[3]][f[7]].append((s0, e0))
    out = {}
    for w, d in iv.items():
        best = None
        for lab, ivs in d.items():
            ivs.sort()
            tot, cs, ce = 0, None, None
            for s0, e0 in ivs:
                if cs is None or s0 > ce:
                    if cs is not None:
                        tot += ce - cs
                    cs, ce = s0, e0
                else:
                    ce = max(ce, e0)
            tot += ce - cs
            if best is None or (tot, lab) > (best[1], best[0]):
                best = (lab, tot)
        out[w] = best
    return out


def read_counts(path):
    """featureCounts output -> (window_id -> reads, total reads over all summary categories)."""
    if not path:
        return None, None
    cnt = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid\t"):
                continue
            f = line.rstrip("\n").split("\t")
            cnt[f[0]] = int(f[6])
    try:                                   # host nuclear library size (no viral genomes, no chrM), written by step 16
        with open(path + ".host_total") as fh:
            return cnt, int(fh.read().strip())
    except (OSError, ValueError):
        pass
    tot = 0
    with open(path + ".summary") as fh:
        next(fh)
        for line in fh:
            tot += int(line.rstrip("\n").split("\t")[1])
    return cnt, tot


def cmd_assemble(a):
    win = []              # (chrom, start, end, id, set, draw, peak_id, neglog10_q)
    with open(a.windows) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            wid = f[3]
            if "#bg" in wid:
                pk, dr = wid.rsplit("#bg", 1)
                win.append((f[0], int(f[1]), int(f[2]), wid, "background", dr, pk, "NA"))
            else:
                win.append((f[0], int(f[1]), int(f[2]), wid, "peak", "0", wid, f[4] if len(f) > 4 else "NA"))
    in_model = {l.split("\t")[3].rstrip("\n") for l in open(a.model_bed)}
    gc = {}
    with open(a.nuc) as fh:
        hdr = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                gi = next(i for i, h in enumerate(f) if h.endswith("pct_gc"))
                continue
            gc[f[3]] = f[gi]
    cs, rm = best_overlap(a.censat), best_overlap(a.repeatmasker)
    region = {}
    if a.regions:
        with open(a.regions) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                region.setdefault(f[3], f[-1])
    ip, ip_tot = read_counts(a.ip_counts)
    ct, ct_tot = read_counts(a.ctl_counts)

    with gzip.open(a.out + ".windows.tsv.gz", "wt") as fh:
        fh.write("sample_id\twindow_id\tset\tbg_draw\tpeak_id\tin_model\tchrom\tstart\tend\twidth\tgc\tneglog10_q\tregion_id\t"
                 "censat_class\tcensat_frac\trm_class\trm_frac\tip_reads\tctl_reads\n")
        for c, s, e, wid, st, dr, pk, q in win:
            w = e - s
            cl, cb = cs.get(wid, ("none", 0))
            rl, rb = rm.get(wid, ("none", 0))
            reg = region.get(wid, "NA") if st == "peak" else "NA"
            fh.write("\t".join(map(str, [
                a.sample, wid, st, dr, pk, int(pk in in_model), c, s, e, w, gc.get(wid, "NA"), q, reg,
                cl, f"{cb / w:.4f}", rl, f"{rb / w:.4f}",
                "NA" if ip is None else ip.get(wid, 0), "NA" if ct is None else ct.get(wid, 0)])) + "\n")
    with open(a.out + ".libsize.tsv", "w") as fh:
        fh.write(f"{a.sample}\t{'NA' if ip_tot is None else ip_tot}\t{'NA' if ct_tot is None else ct_tot}\n")

    # ---- motif sites per window. Both inputs are motif-major (all rows of one family together), so one family is held
    # in memory at a time; satellite motifs can have millions of sites.
    wins = {w[3] for w in win}

    def peak_rows():                                    # chrom start end family window strand (genomic coordinates)
        with gzip.open(a.peak_sites, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[4] in wins:
                    yield f[3], f[4], int(f[1]), int(f[2]), f[5]

    def bg_rows():                                      # FIMO --text on the background windows
        idx = None
        with opener(a.bg_fimo) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if idx is None:
                    idx = {h: i for i, h in enumerate(f)}
                    if not all(k in idx for k in ("motif_id", "sequence_name", "start", "stop", "strand")):
                        sys.exit(f"unexpected FIMO header: {f}")
                    continue
                yield (f[idx["motif_id"]], f[idx["sequence_name"]].split("::")[0], int(f[idx["start"]]),
                       int(f[idx["stop"]]), f[idx["strand"]])

    # every background site must belong to a background window of this sample (see summarize_fimo.py)
    bg_names = {w for w in wins if "#bg" in w}
    bad = set()
    for _fam, w, *_ in bg_rows():
        if w not in bg_names:
            bad.add(w)
    if bad:
        sys.exit(f"[model_windows] ERROR ({a.sample}): FIMO reported {len(bad)} sequence names that are not background "
                 f"windows (e.g. {', '.join(sorted(bad)[:3])}). FIMO probably parsed the FASTA headers as genomic "
                 "coordinates; use a FIMO that accepts --no-pgc.")
    with gzip.open(a.out + ".sites.tsv.gz", "wt") as out:
        out.write("sample_id\twindow_id\tfamily_id\tn_sites\tn_tandem_sites\n")

        def flush(fam, per):
            for w, sites in per.items():
                by = defaultdict(list)
                for (st, _sp), strand in sites.items():
                    by[strand].append(st)
                t = sum(tandem_sites(sorted(v), a.spacing_max) for v in by.values())
                out.write(f"{a.sample}\t{w}\t{fam}\t{len(sites)}\t{t}\n")

        for rows in (peak_rows(), bg_rows()):
            cur, per = None, defaultdict(dict)
            for fam, w, st, sp, strand in rows:
                if fam != cur:
                    if cur is not None:
                        flush(cur, per)
                    cur, per = fam, defaultdict(dict)
                per[w].setdefault((st, sp), strand)         # the same position on both strands (palindrome) counts once
            if cur is not None:
                flush(cur, per)


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
sub = ap.add_subparsers(dest="cmd", required=True)
s = sub.add_parser("select")
s.add_argument("--sample", required=True); s.add_argument("--scan-bed", required=True)
s.add_argument("--peaks-file", required=True); s.add_argument("--out", required=True)
s.add_argument("--max-q", type=float, default=None); s.add_argument("--max-peaks", type=int, default=0)
s.add_argument("--seed", default="1")
s.set_defaults(fn=cmd_select)
s = sub.add_parser("assemble")
s.add_argument("--sample", required=True); s.add_argument("--windows", required=True)
s.add_argument("--model-bed", required=True); s.add_argument("--nuc", required=True)
s.add_argument("--censat", default=""); s.add_argument("--repeatmasker", default=""); s.add_argument("--regions", default="")
s.add_argument("--ip-counts", default=""); s.add_argument("--ctl-counts", default="")
s.add_argument("--peak-sites", required=True); s.add_argument("--bg-fimo", required=True)
s.add_argument("--spacing-max", type=int, default=400); s.add_argument("--out", required=True)
s.set_defaults(fn=cmd_assemble)
a = ap.parse_args()
a.fn(a)
