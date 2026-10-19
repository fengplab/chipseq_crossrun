#!/usr/bin/env python3
"""Run a command and sample the resources of its whole process tree (Linux /proc; no extra packages).

Usage: profile_run.py --out stepNN.profile.tsv [--interval 5] -- command args...

Every --interval seconds it records, summed over the command and all of its descendants:
  rss_gb    memory in use: PSS (shared pages split between processes), RSS if PSS is unreadable
  cores     CPU in use (CPU seconds consumed since the last sample / seconds elapsed)
  nproc     processes alive
  read_gb / write_gb  cumulative disk I/O (when /proc/<pid>/io is readable)
A final '#summary' line gives wall time, peak RSS, peak and mean cores and total CPU hours.
The command's exit status is passed through. Sampling costs well under 1% of one core.

Why not /usr/bin/time -v: its "Maximum resident set size" is the largest SINGLE process, so it
misses memory when the step runs several samples at once (xargs -P, mclapply), which most steps do.
"""
import argparse, os, subprocess, sys, time

TCK = os.sysconf("SC_CLK_TCK")
PAGE = os.sysconf("SC_PAGE_SIZE")


def children_map():
    kids = {}
    for p in os.listdir("/proc"):
        if not p.isdigit():
            continue
        try:
            with open(f"/proc/{p}/stat") as fh:
                s = fh.read()
            ppid = int(s[s.rfind(")") + 2:].split()[1])
            kids.setdefault(ppid, []).append(int(p))
        except (OSError, ValueError, IndexError):
            pass
    return kids


def tree(root):
    kids, out, todo = children_map(), [], [root]
    while todo:
        p = todo.pop()
        out.append(p)
        todo.extend(kids.get(p, []))
    return out


def proc_stats(pid):
    """(cpu_seconds, rss_bytes, read_bytes, write_bytes) of one process, or None if it is gone."""
    try:
        with open(f"/proc/{pid}/stat") as fh:
            s = fh.read()
        f = s[s.rfind(")") + 2:].split()
        cpu = (int(f[11]) + int(f[12])) / TCK          # utime + stime
        rss = int(f[21]) * PAGE
    except (OSError, ValueError, IndexError):
        return None
    try:                                              # PSS splits shared pages among the processes sharing
        with open(f"/proc/{pid}/smaps_rollup") as fh: # them, so forked R workers (mclapply) are not counted
            for line in fh:                           # N times over; falls back to RSS if unreadable
                if line.startswith("Pss:"):
                    rss = int(line.split()[1]) * 1024
                    break
    except OSError:
        pass
    rb = wb = 0
    try:
        with open(f"/proc/{pid}/io") as fh:
            for line in fh:
                if line.startswith("read_bytes:"):
                    rb = int(line.split()[1])
                elif line.startswith("write_bytes:"):
                    wb = int(line.split()[1])
    except OSError:
        pass
    return cpu, rss, rb, wb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--interval", type=float, default=5.0)
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        sys.exit("profile_run.py: no command given")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

    t0 = time.time()
    child = subprocess.Popen(cmd)
    seen = {}                       # pid -> last cpu seconds (so CPU of finished processes is kept)
    done_cpu = 0.0                  # CPU of processes that have exited
    io_max = [0, 0]
    peak_rss = peak_cores = 0.0
    last_t, last_cpu = t0, 0.0
    rows = []
    while True:
        alive = child.poll() is None
        now = time.time()
        cur, rss, rb, wb = {}, 0, 0, 0
        for p in tree(child.pid) if alive else []:
            st = proc_stats(p)
            if st:
                cur[p] = st[0]; rss += st[1]; rb += st[2]; wb += st[3]
        for p, c in seen.items():           # processes gone since the last sample
            if p not in cur:
                done_cpu += c
        seen = cur
        cpu_total = done_cpu + sum(cur.values())
        dt = max(now - last_t, 1e-6)
        cores = max(0.0, (cpu_total - last_cpu) / dt) if rows else 0.0
        io_max = [max(io_max[0], rb), max(io_max[1], wb)]
        peak_rss, peak_cores = max(peak_rss, rss), max(peak_cores, cores)
        rows.append((round(now - t0, 1), round(rss / 2**30, 3), round(cores, 2), len(cur),
                     round(io_max[0] / 2**30, 3), round(io_max[1] / 2**30, 3)))
        last_t, last_cpu = now, cpu_total
        if not alive:
            break
        time.sleep(a.interval)
    # children reaped by the kernel: add their CPU via rusage (covers processes shorter than one interval)
    import resource
    cpu_ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu_hours = max(cpu_total, cpu_ru.ru_utime + cpu_ru.ru_stime) / 3600
    wall = time.time() - t0
    with open(a.out, "w") as fh:
        fh.write("t_s\trss_gb\tcores\tnproc\tread_gb\twrite_gb\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
        fh.write(f"#summary\twall_s={wall:.0f}\tpeak_rss_gb={peak_rss / 2**30:.2f}\tpeak_cores={peak_cores:.1f}"
                 f"\tmean_cores={cpu_hours * 3600 / max(wall, 1):.2f}\tcpu_hours={cpu_hours:.2f}"
                 f"\tread_gb={io_max[0] / 2**30:.1f}\twrite_gb={io_max[1] / 2**30:.1f}\texit={child.returncode}\n")
    sys.exit(child.returncode)


if __name__ == "__main__":
    main()
