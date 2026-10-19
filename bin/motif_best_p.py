#!/usr/bin/env python3
"""
Best p-value FIMO can assign to any site of each motif (a perfect match to the motif's preferred base at every
position), under the motif file's background letter frequencies:

    best_p = product over positions of the background probability of the highest-scoring base(s)

A motif whose best_p is above the FIMO threshold (motifs.fimo_thresh) can never be reported by FIMO, however
often it occurs: typically short (6-7 bp) or degenerate motifs. FIMO's pseudocount (0.1) and score rounding are
ignored; they change best_p only marginally.

Usage: motif_best_p.py family_reps.meme > family_best_p.tsv   (family_id, width, best_possible_pvalue)
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memelib import read_meme  # noqa: E402

bg, motifs = read_meme(sys.argv[1])
letters = "ACGT"
print("family_id\twidth\tbest_possible_pvalue")
for m in motifs:
    logp = 0.0
    for row in m.rows:
        sc = [math.log(max(row[i], 1e-9) / max(bg[a], 1e-9)) for i, a in enumerate(letters)]
        top = max(sc)
        logp += math.log(sum(bg[a] for i, a in enumerate(letters) if sc[i] >= top - 1e-9))
    print(f"{m.id}\t{m.width}\t{math.exp(logp):.4g}")
