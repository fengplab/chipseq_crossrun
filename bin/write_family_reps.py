#!/usr/bin/env python3
"""
Write one representative PWM per motif family (renamed to the family id) and PWM descriptors.
  work/motifs/family_reps.meme
  work/motifs/family_rep_features.tsv  family_id, rep_motif_uid, pwm_consensus, width, info_content,
                                        self_period_shift, self_period_corr, trimer_diversity
Usage: write_family_reps.py <work>
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memelib import read_meme, trimer_diversity, write_meme  # noqa: E402

work = sys.argv[1]
bg, motifs = read_meme(os.path.join(work, "motifs", "all_discovered.meme"))
by_id = {m.id: m for m in motifs}
reps = []
with open(os.path.join(work, "motifs", "families.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["is_representative"].upper() == "TRUE":
            m = by_id[r["motif_uid"]]
            m.alt = r["motif_uid"]
            reps.append((r["family_id"], m))
reps.sort()
out = []
with open(os.path.join(work, "motifs", "family_rep_features.tsv"), "w") as fh:
    fh.write("family_id\trep_motif_uid\tpwm_consensus\twidth\tinfo_content\tself_period_shift\tself_period_corr\ttrimer_diversity\n")
    for fid, m in reps:
        cons = m.consensus()
        k, c = m.self_periodicity()
        fh.write(f"{fid}\t{m.alt}\t{cons}\t{m.width}\t{m.info_content():.2f}\t{k}\t{c:.3f}\t{trimer_diversity(cons):.3f}\n")
        m.id = fid
        out.append(m)
dest = os.path.join(work, "motifs", "family_reps.meme")
write_meme(dest + ".tmp", out, bg)
if os.path.exists(dest) and open(dest).read() == open(dest + ".tmp").read():
    os.remove(dest + ".tmp")          # unchanged: keep the old file (and its timestamp) so scans are reused
else:
    os.replace(dest + ".tmp", dest)
print(f"[reps] wrote {len(out)} family representatives", file=sys.stderr)
