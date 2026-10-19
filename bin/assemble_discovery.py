#!/usr/bin/env python3
"""
Assemble text-mode STREME + MEME + SEA output into the same files XSTREME writes, so the rest of the
bundle (02_parse_xstreme.py onwards) works unchanged:
  <out>/combined.meme     all de novo motifs (STREME then MEME), MEME motif format
  <out>/xstreme.tsv       one row per motif: SOURCE, ID, CONSENSUS, WIDTH, SITES, SEA_PVALUE, EVALUE, ...
  <out>/sea_out/sea.tsv   SEA enrichment of every de novo motif (primary vs control)

Usage: assemble_discovery.py --streme streme.txt --meme meme.txt --out DIR [--sea sea.tsv]
  Step 1 (without --sea) writes combined.meme; step 2 (with --sea) writes xstreme.tsv and sea_out/.
"""
import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memelib import read_meme, write_meme  # noqa: E402


def load(path, source):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None, []
    bg, ms = read_meme(path)
    for m in ms:
        m.source = source
    return bg, ms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--streme")
    ap.add_argument("--meme")
    ap.add_argument("--sea")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    bg1, s_m = load(a.streme, "STREME")
    bg2, m_m = load(a.meme, "MEME")
    motifs, seen = [], set()
    for m in s_m + m_m:
        base, i = m.id, 1
        while m.id in seen:                      # keep IDs unique across the two programs
            i += 1
            m.id = f"{base}_{i}"
        seen.add(m.id)
        motifs.append(m)
    comb = os.path.join(a.out, "combined.meme")
    if not a.sea:
        write_meme(comb, motifs, bg1 or bg2)
        print(f"[assemble] {len(s_m)} STREME + {len(m_m)} MEME motifs", file=sys.stderr)
        return

    sea = {}
    with open(a.sea) as fh:
        hdr = None
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = {h: i for i, h in enumerate(f)}
                continue
            sea[f[hdr["ID"]]] = f[hdr["PVALUE"]] if "PVALUE" in hdr else ""
    os.makedirs(os.path.join(a.out, "sea_out"), exist_ok=True)
    shutil.copy(a.sea, os.path.join(a.out, "sea_out", "sea.tsv"))

    def pv(m):
        try:
            return float(sea.get(m.id, "1"))
        except ValueError:
            return 1.0
    motifs.sort(key=pv)
    with open(os.path.join(a.out, "xstreme.tsv"), "w") as fh:
        fh.write("SEED_MOTIF\tRANK\tSOURCE\tID\tALT_ID\tCONSENSUS\tWIDTH\tSITES\tSEA_PVALUE\tEVALUE\tEVALUE_ACC\tSIM_SOURCE\tSIM_MOTIF\tMOTIF_URL\n")
        for r, m in enumerate(motifs, 1):
            fh.write(f"{m.id}\t{r}\t{m.source}\t{m.id}\t{m.alt}\t{m.consensus()}\t{m.width}\t"
                     f"{int(m.nsites) if m.nsites else ''}\t{sea.get(m.id, '')}\t{m.evalue or ''}\tyes\t\t\t\n")
        fh.write("# written by chipseq_crossrun assemble_discovery.py (text-mode STREME + MEME + SEA)\n")
    print(f"[assemble] xstreme.tsv with {len(motifs)} motifs", file=sys.stderr)


if __name__ == "__main__":
    main()
