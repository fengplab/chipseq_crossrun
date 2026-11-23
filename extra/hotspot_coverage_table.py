#!/usr/bin/env python3
"""Supplemental table: hotspot bins with read coverage, from a finished chipseq_crossrun output folder.

Run after run_all.sh has finished (it needs step 08's bins and R/40_regions.R's table05c). Nothing in the
crossrun is changed; the table is written next to the others.

For each hotspot bin (by default exactly the rows of tables/table05c_hotspot_bins.tsv) it adds:
  mean_cov_ip         mean, over the IP libraries, of each library's mean per-base read depth in the bin,
                      per million counted reads (--norm cpm, default) or raw (--norm raw)
  mean_cov_input      the same for inputs: each IP's own input BAM (control_bam), so an input shared by several
                      IPs counts once per IP and the ratio compares like with like. An input library processed in
                      two pipeline runs is two BAMs; each IP uses the one from its own run
  log2_ip_over_input  log2((mean_cov_ip + pc) / (mean_cov_input + pc)), pc = --pseudocount (default 0.01)
  max_depth           the largest raw per-base read depth in the bin in any IP library (--max-from ip, default)
                      or any library (--max-from all), and max_depth_library, the library it came from

Depth is read depth from `samtools depth -a -b <bins> -Q <min-mapq>`: samtools' default filter (unmapped,
secondary, QC-fail and duplicate records are not counted), supplementary alignments counted, both mates counted
where they overlap. The per-million denominator counts the same records once each:
`samtools view -c -F 0xF04 -q <min-mapq>` (0xF04 = unmapped, secondary, QC-fail, duplicate, supplementary).
--min-mapq 0 (default) uses the BAMs as they are, the reads MACS used; e.g. 255 keeps only reads Bowtie2 placed
uniquely in -k mode (bins.unique_mapq).

Usage (cluster, from the bundle directory, crossrun environment active; extra/hotspot_coverage_table.sbatch):
  python3 extra/hotspot_coverage_table.py --crossrun <outdir> [--bins table05c|all] [--min-samples N]
        [--norm cpm|raw] [--pseudocount 0.01] [--min-mapq 0] [--max-from ip|all] [--threads 8] [--out PREFIX]

Writes (PREFIX default <outdir>/tables/supplemental_hotspot_bins):
  PREFIX.tsv                one row per bin: table05c's columns + the coverage columns
  PREFIX.per_library.tsv    one row per bin x library: raw mean depth, normalised mean depth, max depth
  PREFIX.columns.tsv        what every column of PREFIX.tsv means (data dictionary for the supplement)
  PREFIX.run_info.tsv       settings, samtools version, counted reads per library
"""
import argparse
import bisect
import csv
import math
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

SAMTOOLS = "samtools"


def log(*a):
    print("[hotspot_coverage]", *a, file=sys.stderr, flush=True)


def die(m):
    sys.exit(f"[hotspot_coverage] ERROR: {m}")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run(cmd):
    # errors="replace": some samtools builds print a non-UTF-8 byte in their --version banner
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, errors="replace")
    if p.returncode != 0:
        die(f"{' '.join(cmd)} failed (exit {p.returncode}): {p.stderr.strip()[:500]}")
    return p.stdout


# ------------------------------------------------------------------------------------------------ inputs
def load_bins(cr, mode, min_samples):
    if mode == "table05c":
        f = os.path.join(cr, "tables", "table05c_hotspot_bins.tsv")
        if not os.path.exists(f):
            die(f"{f} not found: run the crossrun through step 12 (R/40_regions.R) first")
        rows = read_tsv(f)
        for k in ("chrom", "start", "end", "n_samples", "total_peaks"):
            if rows and k not in rows[0]:
                die(f"{f} has no column {k}")
        log(f"bins: {len(rows)} rows of table05c_hotspot_bins (same rows and order as the main table)")
    else:
        f = os.path.join(cr, "work", "regions", "bins.tsv")
        if not os.path.exists(f):
            die(f"{f} not found: run step 08 first")
        agg = {}
        for r in read_tsv(f):
            k = (r["chrom"], int(r["start"]), int(r["end"]))
            a = agg.setdefault(k, [set(), 0])
            a[0].add(r["sample_id"]); a[1] += int(r["n_peaks"])
        rows = [{"chrom": c, "start": s, "end": e, "n_samples": len(v[0]), "total_peaks": v[1]}
                for (c, s, e), v in agg.items()]
        n_all = len(rows)
        rows = [r for r in rows if r["n_samples"] >= min_samples]
        # same order as table05c: most samples, then most peaks; ties by position
        rows.sort(key=lambda r: (-r["n_samples"], -r["total_peaks"], r["chrom"], r["start"]))
        log(f"bins: {n_all} bins with a peak in work/regions/bins.tsv -> {len(rows)} with peaks in >= {min_samples} sample(s)")
    if not rows:
        die("no bins to report")
    for r in rows:
        r["start"], r["end"] = int(r["start"]), int(r["end"])
        if r["end"] <= r["start"]:
            die(f"bin with end <= start: {r['chrom']}:{r['start']}-{r['end']}")
    return rows


def load_libraries(cr):
    """IP libraries and their matched inputs, from work/samples.tsv (sample identity from the manifest)."""
    f = os.path.join(cr, "work", "samples.tsv")
    if not os.path.exists(f):
        die(f"{f} not found: is --crossrun a crossrun output folder?")
    ips, inputs, match = {}, {}, {}
    rows = read_tsv(f)
    pairs = []                      # (sample, control_id, control BAM, run tag)
    for r in rows:
        s, bam = r.get("sample_id", ""), (r.get("bam_file") or "").strip()
        if not bam:
            log(f"WARNING: {s} has no bam_file in samples.tsv; left out of the coverage columns")
            continue
        if not os.path.exists(bam):
            die(f"{s}: bam_file {bam} not found")
        ips[s] = bam
        cbam = (r.get("control_bam") or "").strip()
        if cbam:
            if not os.path.exists(cbam):
                die(f"{s}: control_bam {cbam} not found")
            cid = (r.get("control_id") or "").strip() or os.path.basename(cbam)
            # the run the BAM came from: run_id, else the folder above bowtie2/ in an nf-core results path
            parts = cbam.split(os.sep)
            tag = (r.get("run_id") or "").strip() or (parts[parts.index("bowtie2") - 1] if "bowtie2" in parts[1:] else os.path.dirname(cbam))
            pairs.append((s, cid, cbam, tag))
        else:
            log(f"WARNING: {s} has no control_bam; it counts in mean_cov_ip but adds nothing to mean_cov_input")
    # Inputs are keyed by BAM file, not by control_id: one input library processed in two pipeline runs is two BAMs,
    # and each IP is compared with the input BAM from its own run (the one MACS used for its peaks). The display name
    # is the control_id, with the run added when one control_id has more than one BAM.
    bams_of = {}
    for _s, cid, cbam, _t in pairs:
        bams_of.setdefault(cid, set()).add(cbam)
    names = {}
    for s, cid, cbam, tag in pairs:
        name = cid if len(bams_of[cid]) == 1 else f"{cid} [{tag}]"
        if names.get(name, cbam) != cbam:          # two BAMs with the same control_id AND the same run tag
            name = f"{cid} [{cbam}]"
        names[name] = cbam
        inputs[name] = cbam
        match[s] = name
    for cid, b in sorted(bams_of.items()):
        if len(b) > 1:
            log(f"control_id {cid} has {len(b)} BAMs (one per pipeline run); each IP uses the one from its own run: "
                + "; ".join(sorted(n for n in inputs if n.startswith(cid + " ["))))
    if not ips:
        die("no IP library with a bam_file in samples.tsv")
    log(f"libraries: {len(rows)} samples in samples.tsv -> {len(ips)} IP BAMs, {len(inputs)} input BAMs "
        f"({len(match)} IPs with a matched input)")
    for b in list(ips.values()) + list(inputs.values()):
        if not (os.path.exists(b + ".bai") or os.path.exists(os.path.splitext(b)[0] + ".bai")
                or os.path.exists(b + ".csi")):
            die(f"{b} has no index (.bai/.csi); samtools depth -b needs one to read only the bins")
    return ips, inputs, match


def load_alias(cr):
    f = os.path.join(cr, "work", "annot", "chrom_alias.tsv")
    alias = {}
    if os.path.exists(f):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and not line.startswith("#"):
                    alias[p[0]] = p[1]
    else:
        log(f"WARNING: {f} not found; BAM contig names must match the bins' names exactly")
    return alias


def bam_names(bam, alias, need):
    """canonical (crossrun) contig name -> the BAM's own name, for every contig the bins use."""
    out = {}
    for line in run([SAMTOOLS, "view", "-H", bam]).splitlines():
        if line.startswith("@SQ"):
            sn = next((f[3:] for f in line.split("\t") if f.startswith("SN:")), None)
            if sn is not None:
                out.setdefault(alias.get(sn, sn), sn)
    miss = sorted(c for c in need if c not in out)
    if miss:
        die(f"{bam}: no contig matching {', '.join(miss[:5])}{' ...' if len(miss) > 5 else ''} "
            f"(contig names must agree; checked through work/annot/chrom_alias.tsv)")
    return {c: out[c] for c in need}


# ------------------------------------------------------------------------------------------------ per library
def library_coverage(job):
    """mean raw depth, max depth and position count per bin, plus the counted reads, for one BAM."""
    lib, bam, bins, names, min_mapq, threads, tmpdir = job
    bed = os.path.join(tmpdir, f"{lib}.bins.bed")
    by_chrom = {}
    with open(bed, "w") as fh:
        # BED regions: 0-based half-open, the same convention as the bins
        for i, b in enumerate(bins):
            fh.write(f"{names[b['chrom']]}\t{b['start']}\t{b['end']}\n")
            by_chrom.setdefault(names[b["chrom"]], []).append((b["start"], b["end"], i))
    starts = {c: sorted(v) for c, v in by_chrom.items()}
    keys = {c: [s for s, _, _ in v] for c, v in starts.items()}
    n = len(bins)
    tot, mx, npos = [0] * n, [0] * n, [0] * n
    seen = set()
    p = subprocess.Popen([SAMTOOLS, "depth", "-a", "-b", bed, "-Q", str(min_mapq), bam],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, errors="replace",
                         bufsize=1 << 20)
    cur_c, cur = None, None
    for line in p.stdout:
        c, pos, d = line.split("\t")
        # samtools depth prints 1-based positions; a bin is 0-based half-open [start, end), so position pos
        # (1-based) is base pos-1 (0-based) and lies in the bin when start <= pos-1 < end
        x = int(pos) - 1
        if cur is None or c != cur_c or not (cur[0] <= x < cur[1]):
            cur_c, cur = c, None
            v = starts.get(c)
            if v:
                j = bisect.bisect_right(keys[c], x) - 1
                if j >= 0 and v[j][0] <= x < v[j][1]:
                    cur = v[j]
            if cur is None:
                continue
        i = cur[2]
        d = int(d)
        tot[i] += d
        npos[i] += 1
        if d > mx[i]:
            mx[i] = d
        seen.add(i)
    err = p.stderr.read()
    if p.wait() != 0:
        die(f"samtools depth on {bam} failed: {err.strip()[:500]}")
    reads = int(run([SAMTOOLS, "view", "-c", "-F", "0xF04", "-q", str(min_mapq), "-@", str(threads), bam]).strip())
    os.remove(bed)
    return lib, tot, mx, npos, reads


# ------------------------------------------------------------------------------------------------ main
def main():
    global SAMTOOLS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--crossrun", required=True)
    ap.add_argument("--out")
    ap.add_argument("--bins", choices=["table05c", "all"], default="table05c",
                    help="table05c = the main table's rows (default); all = every bin with a peak, see --min-samples")
    ap.add_argument("--min-samples", type=int, default=1, help="--bins all only: bins with peaks in at least N samples")
    ap.add_argument("--norm", choices=["cpm", "raw"], default="cpm")
    ap.add_argument("--pseudocount", type=float, default=0.01)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--max-from", choices=["ip", "all"], default="ip")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--samtools", default="samtools")
    a = ap.parse_args()
    SAMTOOLS = a.samtools
    if a.pseudocount <= 0:
        die("--pseudocount must be above 0")
    if a.min_samples < 1:
        die("--min-samples must be at least 1")
    cr = os.path.abspath(a.crossrun)
    out = os.path.abspath(a.out or os.path.join(cr, "tables", "supplemental_hotspot_bins"))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    sv = run([SAMTOOLS, "--version"]).splitlines()[0]
    log(sv)

    bins = load_bins(cr, a.bins, a.min_samples)
    ips, inputs, match = load_libraries(cr)
    alias = load_alias(cr)
    need = sorted({b["chrom"] for b in bins})
    libs = [(f"IP:{s}", s, b, "IP") for s, b in ips.items()] + [(f"input:{c}", c, b, "input") for c, b in inputs.items()]
    tmpdir = out + ".tmp_bins"
    os.makedirs(tmpdir, exist_ok=True)
    jobs = [(key, bam, bins, bam_names(bam, alias, need), a.min_mapq, max(1, a.threads // max(1, len(libs))), tmpdir)
            for key, _name, bam, _kind in libs]
    workers = max(1, min(a.threads, len(jobs)))
    log(f"{len(bins)} bins x {len(jobs)} libraries; {workers} at a time; min MAPQ {a.min_mapq}")
    res = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for lib, tot, mx, npos, reads in ex.map(library_coverage, jobs):
            if reads <= 0:
                die(f"{lib}: no reads pass the filter (samtools view -c -F 0xF04 -q {a.min_mapq})")
            nb = sum(1 for v in tot if v > 0)
            log(f"{lib}: {reads:,} counted reads; {nb} of {len(bins)} bins with any coverage")
            short = sum(1 for i, b in enumerate(bins) if npos[i] != b["end"] - b["start"])
            if short:
                # samtools depth -a skips contigs where the library has no reads at all; those positions are depth 0
                # and the mean divides by the full bin length, so the values are right either way
                log(f"{lib}: {short} bin(s) partly or wholly on a contig with no reads in this library (depth 0 there)")
            res[lib] = (tot, mx, reads)
    os.rmdir(tmpdir)
    if not any(any(v > 0 for v in res[k][0]) for k in res):
        die("no library has any coverage in any bin: check BAMs, contig names and --min-mapq")

    def norm_mean(lib, i):
        tot, _mx, reads = res[lib]
        m = tot[i] / (bins[i]["end"] - bins[i]["start"])
        return m / (reads / 1e6) if a.norm == "cpm" else m

    unit = "mean per-base read depth per million counted reads" if a.norm == "cpm" else "mean per-base read depth"
    rows, per = [], []
    ip_keys = [f"IP:{s}" for s in ips]
    for i, b in enumerate(bins):
        ipv = [norm_mean(k, i) for k in ip_keys]
        inv = [norm_mean(f"input:{match[s]}", i) for s in ips if s in match]
        m_ip = sum(ipv) / len(ipv)
        m_in = sum(inv) / len(inv) if inv else float("nan")
        lr = math.log2((m_ip + a.pseudocount) / (m_in + a.pseudocount)) if inv else float("nan")
        pool = ip_keys if a.max_from == "ip" else list(res)
        best = max(pool, key=lambda k: res[k][1][i])
        r = {k: b[k] for k in b}
        r.update({"mean_cov_ip": m_ip, "mean_cov_input": m_in, "log2_ip_over_input": lr,
                  "max_depth": res[best][1][i],
                  "max_depth_library": best.split(":", 1)[1] if res[best][1][i] > 0 else "NA"})
        rows.append(r)
        for key, name, _bam, kind in libs:
            tot, mx, reads = res[key]
            per.append({"chrom": b["chrom"], "start": b["start"], "end": b["end"], "library": name, "library_type": kind,
                        "mean_depth_raw": tot[i] / (b["end"] - b["start"]),
                        "mean_depth_per_million_reads": tot[i] / (b["end"] - b["start"]) / (reads / 1e6),
                        "max_depth": mx[i], "counted_reads": reads})

    def fmt(v):
        if isinstance(v, float):
            return "NA" if math.isnan(v) else f"{v:.6g}"
        return str(v)

    def write(path, data, cols):
        with open(path + ".tmp", "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in data:
                fh.write("\t".join(fmt(r.get(c, "")) for c in cols) + "\n")
        os.replace(path + ".tmp", path)
        log(f"wrote {path} ({len(data)} rows)")

    base = [c for c in rows[0] if c not in ("mean_cov_ip", "mean_cov_input", "log2_ip_over_input", "max_depth", "max_depth_library")]
    cols = base + ["mean_cov_ip", "mean_cov_input", "log2_ip_over_input", "max_depth", "max_depth_library"]
    write(out + ".tsv", rows, cols)
    write(out + ".per_library.tsv", per, ["chrom", "start", "end", "library", "library_type", "mean_depth_raw",
                                          "mean_depth_per_million_reads", "max_depth", "counted_reads"])
    desc = {
        "chrom": "Chromosome (crossrun contig names)",
        "start": "Bin start, 0-based (BED convention); the bin covers start+1 to end in 1-based coordinates",
        "end": "Bin end, exclusive (BED convention)",
        "n_samples": "Number of IP samples with at least one peak in the bin",
        "total_peaks": "Merged peaks in the bin, summed over samples",
        "mean_cov_ip": f"Mean over the {len(ips)} IP libraries of the {unit} in the bin",
        "mean_cov_input": f"Mean over the IP libraries' matched inputs (an input shared by several IPs counts once per IP) of the {unit} in the bin",
        "log2_ip_over_input": f"log2((mean_cov_ip + {a.pseudocount:g}) / (mean_cov_input + {a.pseudocount:g}))",
        "max_depth": "Largest per-base read depth (raw, not normalised) in the bin in any "
                     + ("IP library" if a.max_from == "ip" else "library, IP or input"),
        "max_depth_library": "The library holding max_depth (sample_id for an IP; control_id for an input, with [run] added when one control_id has a BAM in several pipeline runs)",
    }
    with open(out + ".columns.tsv", "w") as fh:
        fh.write("column\tdescription\n")
        for c in cols:
            fh.write(f"{c}\t{desc.get(c, 'From table05c_hotspot_bins')}\n")
    log(f"wrote {out}.columns.tsv")
    with open(out + ".run_info.tsv", "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in [("crossrun", cr), ("bins", a.bins), ("min_samples", a.min_samples if a.bins == "all" else "NA"),
                     ("norm", a.norm), ("pseudocount", a.pseudocount), ("min_mapq", a.min_mapq),
                     ("max_from", a.max_from), ("samtools", sv),
                     ("depth_command", f"samtools depth -a -b <bins> -Q {a.min_mapq}"),
                     ("read_count_command", f"samtools view -c -F 0xF04 -q {a.min_mapq}")]:
            fh.write(f"{k}\t{v}\n")
        for k in res:
            fh.write(f"counted_reads:{k}\t{res[k][2]}\n")
    log(f"wrote {out}.run_info.tsv")


if __name__ == "__main__":
    main()
