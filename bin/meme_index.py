#!/usr/bin/env python3
"""Index known-motif MEME files: motif_id, alt_id (TF name), consensus, width, db.
Usage: meme_index.py db1.meme [db2.meme ...] > known_motif_index.tsv"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memelib import read_meme  # noqa: E402

print("motif_id\talt_id\tconsensus\twidth\tdb")
for db in sys.argv[1:]:
    _, ms = read_meme(db)
    for m in ms:
        print(f"{m.id}\t{m.alt or m.id}\t{m.consensus()}\t{m.width}\t{os.path.basename(db)}")
