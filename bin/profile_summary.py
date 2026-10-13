#!/usr/bin/env python3
"""Summarise a profiled run (CROSSRUN_PROFILE=1) into one table per step plus a resource suggestion.

Usage: profile_summary.py <outdir> [--threads N]
Reads <outdir>/logs/profile/stepNN.tsv, writes <outdir>/logs/profile/summary.tsv and prints it.

Columns: step, wall time, peak memory, peak and mean cores, CPU hours, and
  efficiency = mean cores / threads requested: near 1 means the step used what it was given; low values mean
  the step is serial (R, Python) or I/O-bound, so more threads will not speed it up.
"""
import glob, os, re, sys

NAMES = {0: "setup", 1: "prepare inputs", 2: "motif discovery", 3: "sequences", 4: "motif families", 5: "motif scans",
         6: "repeat enrichment", 7: "reads in repeats", 8: "regions", 9: "R: overview + Q1", 10: "R: Q2 repeats",
         11: "R: Q3 known motifs", 12: "R: Q4 regions", 13: "R: Q5 metadata", 14: "QC computations", 15: "R: QC",
         16: "model features", 17: "R: models", 18: "report", 19: "extra figures"}


def hms(s):
    s = int(s); return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    out = sys.argv[1]
    threads = None
    if "--threads" in sys.argv:
        threads = int(sys.argv[sys.argv.index("--threads") + 1])
    else:
        env = os.path.join(out, "work", "config.env")
        if os.path.exists(env):
            m = re.search(r"THREADS=['\"]?(\d+)", open(env).read())
            threads = int(m.group(1)) if m else None
    rows = []
    for f in sorted(glob.glob(os.path.join(out, "logs", "profile", "step*.tsv"))):
        n = int(re.search(r"step(\d+)", f).group(1))
        summ = [l for l in open(f) if l.startswith("#summary")]
        if not summ:
            continue
        d = dict(kv.split("=") for kv in summ[-1].strip().split("\t")[1:])
        d = {k: float(v) for k, v in d.items()}
        rows.append((n, d))
    if not rows:
        sys.exit(f"no profiles in {out}/logs/profile (run with CROSSRUN_PROFILE=1)")
    hdr = ["step", "name", "wall", "peak_mem_gb", "peak_cores", "mean_cores", "cpu_hours", "efficiency", "read_gb", "exit"]
    lines = ["\t".join(hdr)]
    tot_wall = tot_cpu = 0.0; peak_mem = peak_c = 0.0
    for n, d in rows:
        eff = f"{d['mean_cores'] / threads:.2f}" if threads else "NA"
        lines.append("\t".join([str(n), NAMES.get(n, ""), hms(d["wall_s"]), f"{d['peak_rss_gb']:.1f}", f"{d['peak_cores']:.1f}",
                                f"{d['mean_cores']:.1f}", f"{d['cpu_hours']:.2f}", eff, f"{d.get('read_gb', 0):.1f}", str(int(d['exit']))]))
        tot_wall += d["wall_s"]; tot_cpu += d["cpu_hours"]
        peak_mem = max(peak_mem, d["peak_rss_gb"]); peak_c = max(peak_c, d["peak_cores"])
    path = os.path.join(out, "logs", "profile", "summary.tsv")
    open(path, "w").write("\n".join(lines) + "\n")
    w = [max(len(l.split("\t")[i]) for l in lines) for i in range(len(hdr))]
    for l in lines:
        print("  ".join(c.ljust(w[i]) for i, c in enumerate(l.split("\t"))))
    print(f"\ntotal wall {hms(tot_wall)}, {tot_cpu:.1f} CPU hours, threads = {threads or '?'}")
    # suggestion: memory = peak + 25% headroom (min +2 GB); cores = what the parallel steps actually reached
    mem = max(peak_mem * 1.25, peak_mem + 2)
    print(f"suggested request for this manifest and config: --mem={int(mem + 0.999)}G  "
          f"--cpus-per-task={max(1, min(threads or 99, int(peak_c + 0.999)))}  --time={hms(tot_wall * 1.5)}")
    if threads:
        bound = [f"{n} {NAMES.get(n, '')}" for n, d in rows if d["wall_s"] >= 60 and d["mean_cores"] >= 0.75 * threads]
        if bound:
            print("steps using all their threads (more threads would shorten them, up to one per sample): " + ", ".join(bound))
    slow = sorted(rows, key=lambda r: -r[1]["wall_s"])[:3]
    print("longest steps: " + ", ".join(f"{n} {NAMES.get(n, '')} ({hms(d['wall_s'])})" for n, d in slow))
    print(f"written: {path}")


if __name__ == "__main__":
    main()
