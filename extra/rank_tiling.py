#!/usr/bin/env python3
"""PROTOTYPE, not part of the bundle: rank motif families by how strongly their sites tile.

Compares the +/-1 bp window at each multiple of the modal spacing against a window the same
width placed half a period away. The off-peak window is a local control, so the comparison is
not affected by the shape of the spacing distribution or by the size of the period -- which is
what makes periodic_frac incomparable across families.

  tiling_rungs  how many of k = 1..K have on-peak >= --min-ratio x off-peak
  tiling_ratio  median on/off ratio over those k
Usage: python3 rank_tiling.py <crossrun> [--k 6] [--min-ratio 3] [--min-spacings 30]
"""
import csv, os, sys, math
from collections import defaultdict

cr = sys.argv[1]
K = int(sys.argv[sys.argv.index("--k") + 1]) if "--k" in sys.argv else 6
MINR = float(sys.argv[sys.argv.index("--min-ratio") + 1]) if "--min-ratio" in sys.argv else 3.0
MINS = int(sys.argv[sys.argv.index("--min-spacings") + 1]) if "--min-spacings" in sys.argv else 30

def rd(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

rl = {r["family_id"]: r for r in rd(os.path.join(cr, "tables", "table02_repeat_likeness.tsv"))}
hist = defaultdict(lambda: defaultdict(float))
sdir = os.path.join(cr, "work", "scan")
for f in sorted(x for x in os.listdir(sdir) if x.endswith(".peaks.spacing.tsv")):
    for r in rd(os.path.join(sdir, f)):
        d = int(float(r["spacing"]))
        if d >= 2:
            hist[r["family_id"]][d] += float(r["n"])

def win(h, centre):
    return sum(v for d, v in h.items() if abs(d - centre) <= 1)

rows = []
for fam, h in hist.items():
    meta = rl.get(fam, {})
    tot = sum(h.values())
    try:
        mode = int(float(meta.get("spacing_mode", "nan")))
    except ValueError:
        continue
    if tot < MINS or mode < 3:
        continue
    # a long period runs out of axis: only k with the off-peak window still inside spacing_max
    # can be scored, so rungs are reported out of what was achievable, not out of K
    smax = max(h)
    k_possible = min(K, int(smax / mode - 0.5))
    rungs, ratios = 0, []
    for k in range(1, k_possible + 1):
        on = win(h, k * mode)
        off = win(h, int(round((k + 0.5) * mode)))
        if on <= 0:
            continue
        ratio = on / off if off > 0 else float("inf")
        ratios.append(ratio)
        if ratio >= MINR:
            rungs += 1
    if not ratios:
        continue
    fin = sorted(r for r in ratios if math.isfinite(r))
    med = fin[len(fin) // 2] if fin else float("inf")
    rows.append((rungs / max(k_possible, 1), rungs, k_possible, med, fam, mode, tot,
                 meta.get("mean_sites_per_hit_peak", ""), meta.get("periodic_frac", ""),
                 meta.get("spacing_mode_frac", ""), meta.get("flag", "")))

rows.sort(key=lambda r: (-r[0], -(r[3] if math.isfinite(r[3]) else 1e9)))
print(f"{'family':8s}{'rungs':>10s}{'on/off':>9s}{'mode':>6s}{'pairs':>9s}"
      f"{'sites/pk':>10s}{'per_frac':>10s}{'mode_frac':>11s}  flag")
for frac, rungs, kposs, med, fam, mode, tot, msp, pf, smf, flag in rows[:12]:
    m = "inf" if not math.isfinite(med) else f"{med:.1f}"
    print(f"{fam:8s}{f'{rungs}/{kposs}':>10s}{m:>9s}{mode:6d}{tot:9.0f}{msp:>10s}{pf:>10s}{smf:>11s}  {flag}")
print(f"\n{len(rows)} families scored (n_spacings >= {MINS}, mode >= 3), K = {K}, min ratio = {MINR}")
