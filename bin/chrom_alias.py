#!/usr/bin/env python3
"""
Build a chromosome-name map so every input (peaks, BigBeds, BEDs, BAMs, FASTA) can be expressed in one
canonical style: UCSC 'chrN' (chr1..chr22, chrX, chrY, chrM).

Alias sources, first match wins:
  1. --alias files given by the user (any number):
       - NCBI *_assembly_report.txt (GenBank, RefSeq and sequence names -> UCSC-style-name column)
       - UCSC chromAlias.txt (hub style '# ucsc genbank refseq ...' header, or 'alias<TAB>ucsc<TAB>source')
       - any TSV with >= 2 columns; the column whose values mostly start with 'chr' is taken as canonical
  2. FASTA headers of the genome ('>NC_060925.1 Homo sapiens isolate CHM13 chromosome 1, ...')
  3. Built-in T2T-CHM13v2.0 table (RefSeq NC_0609xx.1 and GenBank CP0682xx.2 accessions)
  4. Generic rules: '1' -> 'chr1', 'X' -> 'chrX', 'MT'/'M' -> 'chrM', 'chr*' unchanged

Writes (into --outdir):
  chrom_alias.tsv   alias <TAB> canonical   (includes canonical -> canonical)
  fasta_names.tsv   canonical <TAB> name used in the genome FASTA
  genome.sizes      canonical chrom sizes of the host contigs used by the peak-based analyses: sequences present
                    in the FASTA, minus --host-exclude (e.g. chrM)
  genome_all.sizes  every contig reads can align to: all host sequences plus the viral genomes (display names)
  contigs.tsv       contig, length, group (nuclear / excluded / viral), name in the BAMs or FASTA
  chrom_naming_report.tsv   how each FASTA / chrom-sizes name was resolved
"""
import argparse
import os
import re
import sys


def t2t_chm13v2():
    m = {}
    for n in range(1, 23):
        m[f"NC_{60924 + n:06d}.1"] = f"chr{n}"
        m[f"CP{68278 - n:06d}.2"] = f"chr{n}"
    m.update({"NC_060947.1": "chrX", "CP068255.2": "chrX",
              "NC_060948.1": "chrY", "CP086569.2": "chrY",
              "NC_012920.1": "chrM", "CP068254.1": "chrM"})
    # also accept versionless accessions
    for k, v in list(m.items()):
        m.setdefault(k.split(".")[0], v)
    return m


def generic(name):
    if name.startswith("chr"):
        return name
    if re.fullmatch(r"\d{1,2}|X|Y", name):
        return "chr" + name
    if name in ("MT", "M", "chrMT"):
        return "chrM"
    return None


def from_header(desc):
    """'Homo sapiens isolate CHM13 chromosome 1, alternate assembly' -> chr1"""
    d = desc.lower()
    if "mitochondri" in d:
        return "chrM"
    if "unplaced" in d or "unlocalized" in d or "scaffold" in d:
        return None
    m = re.search(r"chromosome\s+([0-9]{1,2}|x|y)\b", d)
    return f"chr{m.group(1).upper()}" if m else None


def read_alias_file(path):
    m = {}
    with open(path) as fh:
        lines = [l.rstrip("\n\r") for l in fh if l.strip()]
    # NCBI assembly report
    hdr = next((l for l in lines if l.startswith("# Sequence-Name")), None)
    if hdr:
        cols = [c.strip() for c in hdr.lstrip("# ").split("\t")]
        idx = {c: i for i, c in enumerate(cols)}
        for l in lines:
            if l.startswith("#"):
                continue
            f = l.split("\t")
            ucsc = f[idx["UCSC-style-name"]] if "UCSC-style-name" in idx else "na"
            if ucsc in ("na", ""):
                ucsc = generic(f[idx["Sequence-Name"]]) or None
            if not ucsc:
                continue
            for c in ("Sequence-Name", "GenBank-Accn", "RefSeq-Accn"):
                if c in idx and f[idx[c]] not in ("na", ""):
                    m[f[idx[c]]] = ucsc
        return m
    # UCSC hub chromAlias.txt: '# ucsc\tgenbank\trefseq...'
    if lines[0].startswith("#") and "ucsc" in lines[0].lower():
        cols = lines[0].lstrip("# ").lower().split("\t")
        u = cols.index("ucsc")
        for l in lines[1:]:
            f = l.split("\t")
            for i, v in enumerate(f):
                if v:
                    m[v] = f[u]
        return m
    rows = [l.split("\t") for l in lines if not l.startswith("#")]
    ncol = min(len(r) for r in rows)
    if ncol < 2:
        sys.exit(f"[chrom_alias] {path}: need at least 2 tab-separated columns")
    # pick the column that looks most like UCSC names
    score = [sum(r[i].startswith("chr") for r in rows) for i in range(ncol)]
    c = score.index(max(score))
    if max(score) == 0:
        sys.exit(f"[chrom_alias] {path}: no column contains 'chr*' names")
    for r in rows:
        for i, v in enumerate(r[:ncol]):
            if v and r[c]:
                m[v] = r[c]
    return m


def fasta_headers(fasta):
    """Return {first_word: description} by seeking just before each sequence offset in the .fai."""
    out = {}
    with open(fasta + ".fai") as fai, open(fasta, "rb") as fa:
        for line in fai:
            name, _, off = line.split("\t")[:3]
            off = int(off)
            start = max(0, off - 2000)
            fa.seek(start)
            chunk = fa.read(off - start).decode(errors="replace")
            hl = chunk.rstrip("\n").split("\n")[-1]
            desc = hl[1:].split(None, 1)[1] if hl.startswith(">") and len(hl[1:].split(None, 1)) > 1 else ""
            out[name] = desc
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--sizes", required=True)
    ap.add_argument("--alias", nargs="*", default=[])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--host-exclude", nargs="*", default=[], help="host contigs left out of genome.sizes (canonical names)")
    ap.add_argument("--viral", nargs="*", default=[], help="DISPLAY=CONTIG pairs for viral genomes")
    ap.add_argument("--viral-sizes", help="TSV contig<TAB>length for the viral contigs (FASTA index or BAM header)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    user = {}
    for f in a.alias:
        user.update(read_alias_file(f))
    heads = fasta_headers(a.fasta)
    builtin = t2t_chm13v2()
    report = []

    def resolve(name, desc=""):
        for src, fn in (("alias_file", lambda n: user.get(n)),
                        ("fasta_header", lambda n: from_header(desc) if desc else None),
                        ("builtin_T2T", lambda n: builtin.get(n)),
                        ("generic", generic)):
            v = fn(name)
            if v:
                return v, src
        return name, "unchanged"

    alias = {}
    fasta_names = {}
    for n, desc in heads.items():
        c, src = resolve(n, desc)
        if c in fasta_names and fasta_names[c] != n:
            sys.exit(f"[chrom_alias] two FASTA sequences map to {c}: {fasta_names[c]} and {n}")
        fasta_names[c] = n
        alias[n] = c
        report.append(("fasta", n, c, src))
    sizes = {}
    with open(a.sizes) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) < 2 or not f[1].isdigit():
                continue
            c, src = resolve(f[0], heads.get(f[0], ""))
            alias[f[0]] = c
            sizes[c] = int(f[1])
            report.append(("chrom_sizes", f[0], c, src))
    # all known aliases of canonical names that exist in this genome
    for src_map in (builtin, user):
        for k, v in src_map.items():
            if v in fasta_names:
                alias.setdefault(k, v)
    for c in fasta_names:
        alias[c] = c
        alias.setdefault(c.replace("chr", "", 1), c)
        if c == "chrM":
            alias.setdefault("MT", c)

    missing = [c for c in fasta_names if c not in sizes]
    if missing:
        print(f"[chrom_alias] WARNING: {len(missing)} FASTA sequences absent from chrom sizes (e.g. {missing[:3]})", file=sys.stderr)
    with open(os.path.join(a.outdir, "chrom_alias.tsv"), "w") as fh:
        for k in sorted(alias):
            fh.write(f"{k}\t{alias[k]}\n")
    with open(os.path.join(a.outdir, "fasta_names.tsv"), "w") as fh:
        for c in sorted(fasta_names):
            fh.write(f"{c}\t{fasta_names[c]}\n")
    # viral genomes: display name is the canonical name; the accession (with and without version) maps to it
    viral = dict(v.split("=", 1) for v in a.viral)
    vsizes = {}
    if a.viral_sizes and os.path.exists(a.viral_sizes):
        with open(a.viral_sizes) as fh:
            for l in fh:
                f = l.rstrip("\n").split("\t")
                if len(f) >= 2 and f[1].isdigit():
                    vsizes[f[0]] = int(f[1])
    for disp, acc in viral.items():
        alias[acc] = disp
        alias.setdefault(acc.split(".")[0], disp)
        alias[disp] = disp
        if disp in fasta_names:
            sys.exit(f"[chrom_alias] viral display name {disp} is also a host sequence name")
    with open(os.path.join(a.outdir, "chrom_alias.tsv"), "w") as fh:     # rewritten with the viral names
        for k in sorted(alias):
            fh.write(f"{k}\t{alias[k]}\n")
    excl = set(a.host_exclude)
    unknown_excl = sorted(excl - set(fasta_names))
    if unknown_excl:
        print(f"[chrom_alias] note: host_exclude contigs not in the genome FASTA: {', '.join(unknown_excl)}", file=sys.stderr)
    with open(os.path.join(a.outdir, "genome.sizes"), "w") as fh:
        for c in sorted(sizes):
            if c in fasta_names and c not in excl:
                fh.write(f"{c}\t{sizes[c]}\n")
    missing_v = [d for d, acc in viral.items() if acc not in vsizes]
    if missing_v:
        print(f"[chrom_alias] WARNING: no length for viral contig(s) {', '.join(missing_v)} (not in the viral FASTA or BAM header); "
              "they are left out of the read breakdowns", file=sys.stderr)
    with open(os.path.join(a.outdir, "genome_all.sizes"), "w") as fh, open(os.path.join(a.outdir, "contigs.tsv"), "w") as ct:
        ct.write("contig\tlength\tgroup\tsource_name\n")
        for c in sorted(sizes):
            if c in fasta_names:
                fh.write(f"{c}\t{sizes[c]}\n")
                ct.write(f"{c}\t{sizes[c]}\t{'excluded' if c in excl else 'nuclear'}\t{fasta_names[c]}\n")
        for disp, acc in viral.items():
            if acc in vsizes:
                fh.write(f"{disp}\t{vsizes[acc]}\n")
                ct.write(f"{disp}\t{vsizes[acc]}\tviral\t{acc}\n")
    with open(os.path.join(a.outdir, "chrom_naming_report.tsv"), "w") as fh:
        fh.write("source\toriginal\tcanonical\tresolved_by\n")
        for r in report:
            fh.write("\t".join(r) + "\n")
    non_chr = sorted({c for c in fasta_names if not c.startswith("chr")})
    renamed = sum(1 for r in report if r[0] == "fasta" and r[1] != r[2])
    print(f"[chrom_alias] {len(fasta_names)} FASTA sequences, {renamed} renamed to chr-style"
          + (f"; {len(non_chr)} left as-is (e.g. {non_chr[:3]})" if non_chr else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
