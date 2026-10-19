#!/usr/bin/env python3
"""
Convert a RepeatMasker or CenSat annotation to a BED4 whose 4th column is the feature label used
throughout the analysis:
  repeatmasker -> "class/family" (e.g. SINE/Alu, LINE/L1, Satellite/centr; "Simple_repeat" when equal)
  censat       -> CenSat class (active_hor, hor, dhor, mon, hsat1A, hsat2, hsat3, bsat, gsat, ct, ...)

Accepted inputs (gz ok): BED6+2 from the pipeline (class col 7, family col 8; e.g. bigBedToBed output
of genome/annotation/*.bb), UCSC rmsk.txt, RepeatMasker .out, or any BED4+ whose name encodes the class
('AluY#SINE/Alu' or T2T CenSat names such as 'hsat2_2(...)').

Usage: features_to_classbed.py --kind repeatmasker|censat --genome genome.sizes in.bed > out.bed4
"""
import argparse
import gzip
import re
import sys


def open_any(path):
    if path == "-":
        return sys.stdin
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path)


def censat_class(name):
    base = re.sub(r"\(.*$", "", name)
    base = re.sub(r"(_\d+)+$", "", base)
    return base or name


def rmsk_label(cls, fam):
    cls = cls or "Unknown"
    fam = fam or cls
    cls = cls.rstrip("?")
    fam = fam.rstrip("?")
    return cls if fam == cls else f"{cls}/{fam}"


def parse(line, kind):
    """Return (chrom, start, end, label) or None."""
    s = line.rstrip("\n")
    if not s or s.startswith(("#", "track", "browser")):
        return None
    f = s.split("\t")
    # UCSC rmsk.txt: bin swScore milliDiv milliDel milliIns genoName genoStart genoEnd genoLeft strand repName repClass repFamily ...
    if len(f) >= 13 and f[0].isdigit() and f[6].isdigit() and f[7].isdigit() and f[9] in "+-":
        return f[5], int(f[6]), int(f[7]), rmsk_label(f[11], f[12])
    ws = s.split()
    # RepeatMasker .out: score div del ins query qbegin qend (left) strand repeat class/family ...
    if len(ws) >= 11 and ws[0].isdigit() and ws[5].isdigit() and ws[6].isdigit() and ws[8] in ("+", "C"):
        cf = ws[10]
        cls, fam = (cf.split("/", 1) + [None])[:2]
        return ws[4], int(ws[5]) - 1, int(ws[6]), rmsk_label(cls, fam)
    if len(f) < 3 or not f[1].isdigit() or not f[2].isdigit():
        return None
    chrom, start, end = f[0], int(f[1]), int(f[2])
    name = f[3] if len(f) > 3 else "NA"
    if kind == "censat":
        cls = f[6] if len(f) >= 8 and f[6] not in ("", ".", "NA") and not f[6].isdigit() else censat_class(name)
        return chrom, start, end, cls
    # repeatmasker BED
    if len(f) >= 8 and f[6] not in ("", ".", "NA") and not f[6].replace(".", "").isdigit():
        return chrom, start, end, rmsk_label(f[6], f[7])
    if "#" in name:
        cf = name.split("#", 1)[1]
        cls, fam = (cf.split("/", 1) + [None])[:2]
        return chrom, start, end, rmsk_label(cls, fam)
    return chrom, start, end, name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", choices=["repeatmasker", "censat"], required=True)
    ap.add_argument("--genome", required=True, help="chrom sizes; features on other contigs are dropped")
    ap.add_argument("--alias", help="alias<TAB>canonical table; chromosome names are made canonical")
    ap.add_argument("input")
    a = ap.parse_args()
    ok = {}
    with open(a.genome) as fh:
        for line in fh:
            c, n = line.split("\t")[:2]
            ok[c] = int(n)
    alias = {}
    if a.alias:
        with open(a.alias) as fh:
            for line in fh:
                k, v = line.rstrip("\n").split("\t")[:2]
                alias[k] = v
    out = sys.stdout
    kept = dropped = renamed = 0
    unknown = {}
    with open_any(a.input) as fh:
        for line in fh:
            r = parse(line, a.kind)
            if r is None:
                continue
            chrom, start, end, lab = r
            if chrom in alias:
                renamed += chrom != alias[chrom]
                chrom = alias[chrom]
            if chrom not in ok:
                unknown[chrom] = unknown.get(chrom, 0) + 1
            if chrom not in ok or end <= start:
                dropped += 1
                continue
            end = min(end, ok[chrom])
            lab = lab.replace("\t", "_").replace(" ", "_")
            out.write(f"{chrom}\t{start}\t{end}\t{lab}\n")
            kept += 1
    print(f"[features_to_classbed] {a.kind}: kept {kept}, renamed chromosome on {renamed}, dropped {dropped}", file=sys.stderr)
    if unknown:
        top = sorted(unknown.items(), key=lambda x: -x[1])[:5]
        print(f"[features_to_classbed] {a.kind}: features on contigs not in the genome: "
              + ", ".join(f"{k} ({v})" for k, v in top), file=sys.stderr)
    if kept == 0:
        sys.exit(f"[features_to_classbed] {a.kind}: no features left after chromosome matching; "
                 f"names seen: {', '.join(list(unknown)[:5])}. Supply a chrom_alias file in config.yaml.")


if __name__ == "__main__":
    main()
