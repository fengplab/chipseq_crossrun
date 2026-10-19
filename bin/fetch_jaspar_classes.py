#!/usr/bin/env python3
"""
Optional helper: build the tf_class_table (motif_id, tf_class, tf_family) for JASPAR matrix IDs
(MA####.#) by querying the JASPAR REST API. Non-JASPAR IDs are skipped.
Usage: fetch_jaspar_classes.py known_motif_index.tsv > tf_classes.tsv
"""
import csv
import json
import sys
import time
import urllib.request

API = "https://jaspar.elixir.no/api/v1/matrix/{}/?format=json"
print("motif_id\ttf_class\ttf_family")
with open(sys.argv[1]) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        mid = r["motif_id"]
        if not mid.startswith(("MA", "PB", "PF", "PH", "UN", "CN")):
            continue
        try:
            with urllib.request.urlopen(API.format(mid), timeout=30) as resp:
                d = json.load(resp)
            cls = ";".join(d.get("class") or []) or "NA"
            fam = ";".join(d.get("family") or []) or "NA"
        except Exception as e:  # network or unknown id
            print(f"[jaspar] {mid}: {e}", file=sys.stderr)
            cls = fam = "NA"
        print(f"{mid}\t{cls}\t{fam}")
        time.sleep(0.1)
