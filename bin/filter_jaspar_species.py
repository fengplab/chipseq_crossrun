#!/usr/bin/env python3
"""
Filter a JASPAR motif file (MEME format) to the motifs annotated to given species.

JASPAR's MEME downloads carry no species, so the species come from the JASPAR REST API:
  GET https://jaspar.elixir.no/api/v1/matrix/?tax_id=<id>&collection=<c>&page_size=1000   (paged)
A motif is kept when its matrix ID (e.g. MA0079.5), or failing that its versionless base ID (MA0079), is
annotated to one of the species. Motifs whose IDs are not JASPAR-style are kept and reported.
Motif blocks are copied verbatim, so the output is byte-identical to the input apart from removed motifs.

Usage
  # online (any machine with internet access)
  filter_jaspar_species.py --meme JASPAR2024_CORE_vertebrates_non-redundant.meme --tax-id 9606 \\
      --out JASPAR2024_CORE_vertebrates_human.meme [--write-ids human_jaspar_ids.txt]
  # offline (e.g. a cluster node): reuse an ID list written by --write-ids on a machine with internet
  filter_jaspar_species.py --meme JASPAR...meme --ids-file human_jaspar_ids.txt --out JASPAR..._human.meme

Writes <out> and <out>.species_filter.tsv (motif_id, name, kept, reason).
Common NCBI taxonomy IDs: 9606 human, 10090 mouse, 10116 rat.
"""
import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request

JASPAR_ID = re.compile(r"^(MA|PB|PF|PH|PL|UN|CN|POL)\d{4}\.\d+$")


def fetch_ids(api, tax_ids, collection, release):
    ids = set()
    for tax in tax_ids:
        q = {"tax_id": tax, "page_size": 1000, "format": "json"}
        if collection:
            q["collection"] = collection
        if release:
            q["release"] = release
        url = f"{api.rstrip('/')}/matrix/?{urllib.parse.urlencode(q)}"
        n_pages = 0
        while url:
            for attempt in range(4):
                try:
                    with urllib.request.urlopen(url, timeout=60) as r:
                        d = json.load(r)
                    break
                except Exception as e:  # transient network errors: retry with backoff
                    if attempt == 3:
                        sys.exit(f"[filter_jaspar] could not query the JASPAR API ({e}).\n"
                                 "  On a machine without internet access, run this once elsewhere with --write-ids,\n"
                                 "  copy the ID file over and use --ids-file.")
                    time.sleep(2 ** attempt)
            for m in d.get("results", []):
                ids.add(m["matrix_id"])
            url = d.get("next")
            n_pages += 1
        print(f"[filter_jaspar] tax_id {tax}: {d.get('count', '?')} matrices ({n_pages} page(s))", file=sys.stderr)
    return ids


def split_meme(path):
    """Return (header_text, [(motif_id, name, block_text), ...]) preserving the original text."""
    text = open(path).read()
    parts = re.split(r"(?m)^(?=MOTIF\s)", text)
    header, blocks = parts[0], []
    for b in parts[1:]:
        f = b.split("\n", 1)[0].split()
        blocks.append((f[1], f[2] if len(f) > 2 else "", b))
    return header, blocks


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meme", required=True, help="JASPAR motifs in MEME format")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tax-id", action="append", default=[], help="NCBI taxonomy ID (repeatable); 9606 = human")
    ap.add_argument("--ids-file", help="offline: matrix IDs to keep, one per line (from --write-ids)")
    ap.add_argument("--write-ids", help="also save the fetched matrix IDs to this file for offline reuse")
    ap.add_argument("--collection", default="CORE", help="JASPAR collection to query (default CORE; '' for all)")
    ap.add_argument("--release", default="", help="JASPAR release, e.g. 2024 (default: the API's current release)")
    ap.add_argument("--api", default="https://jaspar.elixir.no/api/v1")
    a = ap.parse_args()

    if a.ids_file:
        keep = {l.strip() for l in open(a.ids_file) if l.strip() and not l.startswith("#")}
        print(f"[filter_jaspar] {len(keep)} matrix IDs from {a.ids_file}", file=sys.stderr)
    elif a.tax_id:
        keep = fetch_ids(a.api, a.tax_id, a.collection, a.release)
        if a.write_ids:
            with open(a.write_ids, "w") as fh:
                fh.write(f"# JASPAR matrix IDs for tax_id {','.join(a.tax_id)}, collection {a.collection or 'all'}\n")
                fh.write("\n".join(sorted(keep)) + "\n")
    else:
        sys.exit("give --tax-id (online) or --ids-file (offline)")
    keep_base = {i.split(".")[0] for i in keep}

    header, blocks = split_meme(a.meme)
    out, report, n_keep, n_nonjaspar = [header], [], 0, 0
    for mid, name, block in blocks:
        if not JASPAR_ID.match(mid):
            ok, why = True, "not a JASPAR ID: kept unfiltered"
            n_nonjaspar += 1
        elif mid in keep:
            ok, why = True, "matrix annotated to species"
        elif mid.split(".")[0] in keep_base:
            ok, why = True, "another version of this profile annotated to species"
        else:
            ok, why = False, "not annotated to species"
        if ok:
            out.append(block if block.endswith("\n") else block + "\n")
            n_keep += 1
        report.append(f"{mid}\t{name}\t{str(ok).upper()}\t{why}")
    if n_keep == n_nonjaspar:
        sys.exit("[filter_jaspar] no JASPAR motif matched the species list; check --tax-id / --collection "
                 "(or that the file really is JASPAR)")
    with open(a.out, "w") as fh:
        fh.write("".join(out))
    with open(a.out + ".species_filter.tsv", "w") as fh:
        fh.write("motif_id\tname\tkept\treason\n" + "\n".join(report) + "\n")
    print(f"[filter_jaspar] kept {n_keep} of {len(blocks)} motifs -> {a.out}"
          + (f" ({n_nonjaspar} non-JASPAR motifs kept unfiltered)" if n_nonjaspar else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
