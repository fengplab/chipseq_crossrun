#!/usr/bin/env python3
"""Dfam EMBL (e.g. Dfam-curated_only-1.embl.gz) -> FASTA of consensus sequences for families whose
clade (OS line) is on the human lineage, i.e. families active in human or any ancestor (MIR, L2, Alu, HSATII...).
Usage: dfam_embl_to_fasta.py Dfam-curated_only-1.embl.gz > Dfam_human_curated.fa"""
import gzip, sys
HUMAN = {"root", "cellular organisms", "Eukaryota", "Opisthokonta", "Metazoa", "Eumetazoa", "Bilateria", "Deuterostomia",
         "Chordata", "Craniata", "Vertebrata", "Gnathostomata", "Teleostomi", "Euteleostomi", "Sarcopterygii",
         "Dipnotetrapodomorpha", "Tetrapoda", "Amniota", "Mammalia", "Theria", "Eutheria", "Boreoeutheria",
         "Euarchontoglires", "Primates", "Haplorrhini", "Simiiformes", "Catarrhini", "Hominoidea", "Hominidae",
         "Homininae", "Homo", "Homo sapiens"}
op = gzip.open if sys.argv[1].endswith(".gz") else open
n_in = n_out = 0
def flush(r):
    global n_out
    if r.get("seq") and HUMAN & set(r["os"]):
        cls = r.get("type", "Unknown") + ("/" + r["sub"] if r.get("sub") else "")
        s = "".join(r["seq"]).upper()
        sys.stdout.write(f">{r['nm']}#{cls} {r.get('ac', '')} clade={','.join(r['os'])}\n")
        sys.stdout.write("\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
        n_out += 1
r, inseq = {"os": []}, False
with op(sys.argv[1], "rt") as fh:
    for line in fh:
        if line.startswith("//"):
            n_in += 1; flush(r); r, inseq = {"os": []}, False
        elif inseq:
            r["seq"].append("".join(c for c in line if c.isalpha()))
        elif line.startswith("ID"): r["ac"] = line[5:].split(";")[0].strip()
        elif line.startswith("NM"): r["nm"] = line[5:].strip()
        elif line.startswith("OS"): r["os"].append(line[5:].strip())
        elif line.startswith("CC") and "Type:" in line and "SubType" not in line: r["type"] = line.split("Type:")[1].strip()
        elif line.startswith("CC") and "SubType:" in line: r["sub"] = line.split("SubType:")[1].strip()
        elif line.startswith("SQ"): inseq, r["seq"] = True, []
print(f"[dfam] {n_out} of {n_in} families on the human lineage", file=sys.stderr)
