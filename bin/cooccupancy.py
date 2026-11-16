#!/usr/bin/env python3
"""Site co-occupancy between motif families: do two families mark the same positions?

Tomtom compares PWMs. Two motifs can be the same repeat read in a different frame or at a different length and
still be left apart by Tomtom, yet their sites sit on top of each other. This measures that directly.

For families A and B (sites pooled over samples, a site on both strands counted once, strand ignored):
    a site of A "co-occupies" with B when it overlaps a site of B by at least
        k = ceil(--overlap x min(width_A, width_B)) bp                      (default --overlap 0.5)
    frac_a_in_b = share of A's sites that co-occupy with B;  frac_b_in_a likewise
    mutual      = min(frac_a_in_b, frac_b_in_a)   <- the clustering similarity: BOTH families must mostly coincide
    containment = max(frac_a_in_b, frac_b_in_a)   (reported only: high containment with low mutual = a short or
                                                   general motif sitting inside a longer, more specific one)
Families are merged by average-linkage clustering on d = 1 - mutual, cut at 1 - threshold: the same linkage the
Tomtom stage uses, so one loose pair cannot chain two groups together.

Subcommands
  assess  existing crossrun output -> pair table, threshold sweep and the groups each threshold would make.
          Reads <outdir>/work/scan/*.peaks.sites.bed.gz (BED, 0-based half-open) and family_rep_features.tsv.
          Changes nothing in the outdir.
              python3 bin/cooccupancy.py assess --crossrun <outdir> [--out DIR] [--overlap 0.5]
  merge   used by step 04 (motifs.family_cooccupancy_min): FIMO --text output of the provisional family
          representatives (1-based closed start/stop, converted here to 0-based half-open) -> families.tsv with
          co-occupying families merged and renumbered.

Needs numpy only.
"""
import argparse
import csv
import glob
import gzip
import math
import os
import sys

import numpy as np

BIG = 10 ** 10            # spacing between sequences on the shared coordinate axis (longer than any chromosome)


def log(*a):
    print("[cooccupancy]", *a, file=sys.stderr, flush=True)


def die(m):
    sys.exit(f"[cooccupancy] ERROR: {m}")


# ------------------------------------------------------------------------------------------------ loading
class Sites:
    """Per family: unique (start, end) on one shared axis, 0-based half-open, sorted by start."""

    def __init__(self):
        self.raw = {}             # family -> list of (starts array, ends array) chunks
        self.keys = {}            # sequence / chromosome name -> index on the shared axis
        self.n_records = 0

    def key_offset(self, name):
        if name not in self.keys:
            self.keys[name] = len(self.keys)
        return self.keys[name] * BIG

    def add(self, fam, starts, ends):
        self.raw.setdefault(fam, []).append((np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64)))
        self.n_records += len(starts)

    def finish(self):
        out = {}
        for fam, chunks in self.raw.items():
            s = np.concatenate([c[0] for c in chunks]); e = np.concatenate([c[1] for c in chunks])
            se = np.unique(np.stack([s, e], axis=1), axis=0)          # same site in two samples / on two strands: once
            out[fam] = (se[:, 0].copy(), se[:, 1].copy())
        self.raw = None
        return out


def load_sites_bed(files):
    """work/scan/<sample>.peaks.sites.bed.gz: chrom start end family_id peak_id strand (BED, 0-based half-open)."""
    S = Sites()
    for f in files:
        per = {}
        with gzip.open(f, "rt") as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) < 4:
                    continue
                d = per.setdefault((p[3], p[0]), ([], []))
                d[0].append(int(p[1])); d[1].append(int(p[2]))
        for (fam, chrom), (st, en) in per.items():
            off = S.key_offset(chrom)
            S.add(fam, np.array(st, dtype=np.int64) + off, np.array(en, dtype=np.int64) + off)
        log(f"read {os.path.basename(f)}")
    return S


def load_sites_fimo(files):
    """FIMO --text: motif_id motif_alt_id sequence_name start stop ... ; start/stop are 1-based closed.
    Converted to 0-based half-open: [start - 1, stop). Each sequence of each file is its own stretch of the axis."""
    S = Sites()
    for f in files:
        op = gzip.open if f.endswith(".gz") else open
        per = {}
        with op(f, "rt") as fh:
            idx = None
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if idx is None:
                    idx = {h: i for i, h in enumerate(p)}
                    if not all(k in idx for k in ("motif_id", "sequence_name", "start", "stop")):
                        die(f"unexpected FIMO header in {f}: {p}")
                    continue
                d = per.setdefault((p[idx["motif_id"]], p[idx["sequence_name"]]), ([], []))
                d[0].append(int(p[idx["start"]]) - 1); d[1].append(int(p[idx["stop"]]))
        for (fam, seq), (st, en) in per.items():
            off = S.key_offset(f + "\t" + seq)
            S.add(fam, np.array(st, dtype=np.int64) + off, np.array(en, dtype=np.int64) + off)
    return S


# ------------------------------------------------------------------------------------------------ pairwise
def frac_overlapping(a, b, k):
    """Share of A's sites that overlap some B site by >= k bp. a, b = (starts, ends) sorted by start."""
    a_s, a_e = a; b_s, b_e = b
    if not len(a_s) or not len(b_s):
        return 0.0, 0
    b_e_sorted = np.sort(b_e)
    # B overlaps [s, e) by >= k bp  <=>  b_start <= e - k  and  b_end >= s + k.
    # (# with b_start <= e - k) - (# with b_end < s + k) counts them, because no site can have both
    # b_start > e - k and b_end < s + k while k <= half the shorter width.
    n = np.searchsorted(b_s, a_e - k, side="right") - np.searchsorted(b_e_sorted, a_s + k, side="left")
    hit = int((n > 0).sum())
    return hit / len(a_s), hit


def pairwise(sites, overlap):
    fams = sorted(sites)
    width = {f: int(np.median(sites[f][1] - sites[f][0])) for f in fams}
    rows = []
    n = len(fams)
    for i in range(n):
        for j in range(i + 1, n):
            A, B = fams[i], fams[j]
            k = max(1, math.ceil(overlap * min(width[A], width[B])))
            fa, ha = frac_overlapping(sites[A], sites[B], k)
            fb, hb = frac_overlapping(sites[B], sites[A], k)
            rows.append({"family_a": A, "family_b": B, "width_a": width[A], "width_b": width[B],
                         "n_sites_a": len(sites[A][0]), "n_sites_b": len(sites[B][0]), "min_overlap_bp": k,
                         "frac_a_in_b": fa, "frac_b_in_a": fb, "mutual": min(fa, fb), "containment": max(fa, fb)})
        if n > 20 and (i + 1) % 20 == 0:
            log(f"pairs: {i + 1}/{n} families done")
    return fams, rows


# ------------------------------------------------------------------------------------------------ clustering
def average_linkage(fams, rows, threshold):
    """Groups from average linkage on d = 1 - mutual, merging while the closest pair of groups has mean
    mutual >= threshold. Pairs with no sites in common have mutual = 0."""
    ix = {f: i for i, f in enumerate(fams)}
    n = len(fams)
    sim = np.zeros((n, n))
    for r in rows:
        sim[ix[r["family_a"]], ix[r["family_b"]]] = sim[ix[r["family_b"]], ix[r["family_a"]]] = r["mutual"]
    groups = [[i] for i in range(n)]
    while len(groups) > 1:
        best, bi, bj = -1.0, -1, -1
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                m = sim[np.ix_(groups[i], groups[j])].mean()
                if m > best:
                    best, bi, bj = m, i, j
        if best < threshold:
            break
        groups[bi] = groups[bi] + groups[bj]
        del groups[bj]
    return [[fams[i] for i in g] for g in groups]


# ------------------------------------------------------------------------------------------------ assess
def cmd_assess(a):
    cr = os.path.abspath(a.crossrun)
    files = sorted(glob.glob(os.path.join(cr, "work", "scan", "*.peaks.sites.bed.gz")))
    if not files:
        die(f"no work/scan/*.peaks.sites.bed.gz under {cr} (step 05 writes them)")
    out = os.path.abspath(a.out or os.path.join(cr, "cooccupancy_assessment"))
    os.makedirs(out, exist_ok=True)
    cons = {}
    ff = os.path.join(cr, "work", "motifs", "family_rep_features.tsv")
    if os.path.exists(ff):
        for r in csv.DictReader(open(ff), delimiter="\t"):
            cons[r["family_id"]] = r.get("pwm_consensus", "")
    S = load_sites_bed(files)
    n_rec = S.n_records
    sites = S.finish()
    n_uniq = sum(len(v[0]) for v in sites.values())
    log(f"{len(files)} samples, {n_rec:,} site records -> {n_uniq:,} unique sites in {len(sites)} families")
    if len(sites) < 2:
        die("fewer than two families have sites; nothing to compare")
    small = [f for f in sites if len(sites[f][0]) < a.min_sites]
    if small:
        log(f"{len(small)} families with fewer than --min-sites {a.min_sites} sites are left out of the comparison: "
            + ", ".join(sorted(small)))
        for f in small:
            del sites[f]
    fams, rows = pairwise(sites, a.overlap)

    cols = ["family_a", "family_b", "consensus_a", "consensus_b", "width_a", "width_b", "n_sites_a", "n_sites_b",
            "min_overlap_bp", "frac_a_in_b", "frac_b_in_a", "mutual", "containment"]
    rows.sort(key=lambda r: -r["mutual"])
    with open(os.path.join(out, "cooccupancy_pairs.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            r["consensus_a"] = cons.get(r["family_a"], ""); r["consensus_b"] = cons.get(r["family_b"], "")
            fh.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

    mu = np.array([r["mutual"] for r in rows])
    print(f"\n{len(fams)} families, {len(rows)} pairs; overlap rule: >= {a.overlap:g} x the shorter motif\n")
    print("Distribution of 'mutual' over all pairs (most pairs should be near 0; pairs that are the same")
    print("element should pile up near 1; a good threshold sits in the empty stretch between them):")
    edges = np.linspace(0, 1, 11)
    h, _ = np.histogram(mu, bins=edges)
    for lo, hi, c in zip(edges[:-1], edges[1:], h):
        print(f"  {lo:.1f}-{hi:.1f}  {c:6d}  {'#' * min(60, int(math.ceil(60 * c / max(1, h.max()))))}")

    print("\nThreshold sweep (average linkage on 1 - mutual):")
    print("  threshold  groups  families_merged  largest_group")
    sweep = []
    for t in [round(x, 2) for x in np.arange(0.10, 0.96, 0.05)]:
        g = average_linkage(fams, rows, t)
        merged = sum(len(x) for x in g if len(x) > 1)
        big = max(len(x) for x in g)
        sweep.append((t, len(g), merged, big, g))
        print(f"  {t:9.2f}  {len(g):6d}  {merged:15d}  {big:13d}")
    with open(os.path.join(out, "cooccupancy_sweep.tsv"), "w") as fh:
        fh.write("threshold\tgroups\tfamilies_merged\tlargest_group\n")
        for t, ng, m, b, _ in sweep:
            fh.write(f"{t}\t{ng}\t{m}\t{b}\n")

    # where is the widest empty stretch among the non-trivial values? (a hint, not a decision)
    v = np.sort(mu[mu >= 0.05])
    if len(v) >= 2:
        gaps = np.diff(v)
        i = int(np.argmax(gaps))
        print(f"\nWidest empty stretch in 'mutual' above 0.05: {v[i]:.2f} to {v[i + 1]:.2f} "
              f"({i + 1} pairs below it, {len(v) - i - 1} above). A threshold inside it gives the same groups "
              f"anywhere in that range.")
    else:
        print("\nFewer than two pairs have mutual >= 0.05: these families do not share sites; nothing to merge.")

    with open(os.path.join(out, "cooccupancy_groups.tsv"), "w") as fh:
        fh.write("threshold\tgroup\tsize\tfamily_id\tconsensus\tn_sites\n")
        for t, _, _, _, g in sweep:
            for gi, x in enumerate(sorted((x for x in g if len(x) > 1), key=len, reverse=True), 1):
                for f in x:
                    fh.write(f"{t}\t{gi}\t{len(x)}\t{f}\t{cons.get(f, '')}\t{len(sites[f][0])}\n")
    for t in a.show:
        g = next((s[4] for s in sweep if abs(s[0] - t) < 1e-9), None) or average_linkage(fams, rows, t)
        multi = sorted((x for x in g if len(x) > 1), key=len, reverse=True)
        print(f"\nGroups at threshold {t:g}: {len(multi)} merged group(s)")
        for x in multi:
            print("  " + " + ".join(f"{f} {cons.get(f, '')}".strip() for f in x))
    print(f"\nTop pairs (all pairs in {os.path.join(out, 'cooccupancy_pairs.tsv')}):")
    print("  family_a  family_b  mutual  a_in_b  b_in_a  consensus_a / consensus_b")
    for r in rows[:a.top]:
        print(f"  {r['family_a']:8s}  {r['family_b']:8s}  {r['mutual']:.2f}    {r['frac_a_in_b']:.2f}    "
              f"{r['frac_b_in_a']:.2f}    {r['consensus_a']} / {r['consensus_b']}")
    print(f"\nWritten to {out}: cooccupancy_pairs.tsv, cooccupancy_sweep.tsv, cooccupancy_groups.tsv")


# ------------------------------------------------------------------------------------------------ merge (step 04)
def cmd_merge(a):
    fam_rows = list(csv.DictReader(open(a.families), delimiter="\t"))
    if not fam_rows:
        die(f"{a.families} is empty")
    cols = list(fam_rows[0].keys())
    S = load_sites_fimo(a.fimo)
    sites = S.finish()
    all_fams = sorted({r["family_id"] for r in fam_rows})
    scanned = {f: v for f, v in sites.items() if f in set(all_fams) and len(v[0]) >= a.min_sites}
    log(f"{len(all_fams)} Tomtom families; {len(scanned)} have at least {a.min_sites} sites in the peak sequences "
        f"and can be compared ({len(all_fams) - len(scanned)} cannot, and stay as they are)")
    groups = [[f] for f in all_fams]
    rows = []
    if len(scanned) >= 2:
        fams, rows = pairwise(scanned, a.overlap)
        groups = average_linkage(fams, rows, a.threshold) + [[f] for f in all_fams if f not in scanned]
    if a.pairs_out:
        pc = ["family_a", "family_b", "width_a", "width_b", "n_sites_a", "n_sites_b", "min_overlap_bp",
              "frac_a_in_b", "frac_b_in_a", "mutual", "containment"]
        with open(a.pairs_out, "w") as fh:
            fh.write("\t".join(pc) + "\n")
            for r in sorted(rows, key=lambda r: -r["mutual"]):
                fh.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in pc) + "\n")
    gid = {f: i for i, g in enumerate(groups) for f in g}

    def num(x):
        try:
            v = float(x)
            return v if math.isfinite(v) else None
        except (TypeError, ValueError):
            return None

    for r in fam_rows:                     # same ranking as R/04_cluster_families.R: SEA p, else E-value, else 1
        sp, ev = num(r.get("sea_pvalue")), num(r.get("evalue"))
        r["_rank"] = sp if sp is not None else (ev if ev is not None else 1.0)
        r["_w"] = num(r.get("width")) or 0
        r["_g"] = gid[r["family_id"]]
        r["tomtom_family"] = r["family_id"]
    by_g = {}
    for r in fam_rows:
        by_g.setdefault(r["_g"], []).append(r)
    order = sorted(by_g, key=lambda g: (-len({r["sample_id"] for r in by_g[g]}), min(r["_rank"] for r in by_g[g])))
    n_merged = 0
    for k, g in enumerate(order, 1):
        members = sorted(by_g[g], key=lambda r: (r["_rank"], -r["_w"]))
        old = sorted({r["tomtom_family"] for r in members})
        if len(old) > 1:
            n_merged += 1
            log(f"merged {' + '.join(old)} -> F{k:03d}")
        ns = len({r["sample_id"] for r in members})
        for i, r in enumerate(members):
            r["family_id"] = f"F{k:03d}"
            r["is_representative"] = "TRUE" if i == 0 else "FALSE"
            r["family_n_samples"] = ns
            r["family_size"] = len(members)
    fam_rows.sort(key=lambda r: (r["family_id"], r["is_representative"] != "TRUE", r["_rank"]))
    if "tomtom_family" not in cols:
        cols.append("tomtom_family")
    tmp = a.out + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(fam_rows)
    os.replace(tmp, a.out)
    log(f"co-occupancy stage: {len(all_fams)} Tomtom families -> {len(order)} families "
        f"({n_merged} merged group(s); mutual >= {a.threshold:g}, overlap >= {a.overlap:g} x the shorter motif)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("assess")
    s.add_argument("--crossrun", required=True, help="existing crossrun output folder (read only)")
    s.add_argument("--out", help="where to write the tables (default <crossrun>/cooccupancy_assessment)")
    s.add_argument("--overlap", type=float, default=0.5, help="sites co-occupy when they overlap by this share of the shorter motif")
    s.add_argument("--min-sites", type=int, default=20, help="families with fewer sites are not compared")
    s.add_argument("--show", type=float, nargs="*", default=[0.3, 0.5, 0.7], help="thresholds whose groups are printed")
    s.add_argument("--top", type=int, default=25)
    s = sp.add_parser("merge")
    s.add_argument("--families", required=True); s.add_argument("--fimo", nargs="+", required=True)
    s.add_argument("--threshold", type=float, required=True); s.add_argument("--overlap", type=float, default=0.5)
    s.add_argument("--min-sites", type=int, default=20)
    s.add_argument("--pairs-out"); s.add_argument("--out", required=True)
    a = ap.parse_args()
    if not 0 < a.overlap <= 1:
        die("--overlap must be above 0 and at most 1")
    if a.cmd == "merge" and not 0 < a.threshold <= 1:
        die("--threshold must be above 0 and at most 1")
    {"assess": cmd_assess, "merge": cmd_merge}[a.cmd](a)


if __name__ == "__main__":
    main()
