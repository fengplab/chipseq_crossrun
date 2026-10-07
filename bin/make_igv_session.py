#!/usr/bin/env python3
"""
Build an IGV Desktop session for one ChIP-seq sample:
  - the sample's reads (BAM; IGV adds the coverage pileup track automatically)
  - the CenSat and RepeatMasker bigBed files
  - the genome (FASTA) and its annotation (GTF)
  - the peaks called on the sample (MACS3 narrowPeak / broadPeak)

Before writing anything it checks that every file uses the same chromosome notation (chr1, chr2, ... as in
UCSC/T2T) as the genome FASTA, and reports each file's notation and how many of its chromosomes match. Names
are read straight from the files: FASTA index, BAM header, the chromosome index inside bigBed files, and
column 1 of the GTF and peak files. Only the Python 3 standard library is needed.

If the notation differs between files, the script stops and explains, unless --resolve is given:
  --resolve alias      write an IGV genome definition (JSON) with a chromosome-alias file, so IGV maps the
                       names onto the FASTA's own names; no data files are copied (IGV shows the FASTA names)
  --resolve chr-fasta  also write a copy of the FASTA with chr-style headers (plus its .fai), so IGV displays
                       chr1, chr2, ...; BAM, peaks, etc. keep their names and are mapped through the alias file
Large derived files are written once, next to the originals, and reused by every later session, so sessions
reference one central copy: chr_<fasta name> (+ .fai) beside the FASTA, and sorted_<gtf name>.gz (+ .tbi), a
sorted, bgzipped, tabix-indexed GTF, beside the GTF. Only the file-name prefix changes, because IGV picks the
format from the extension. They are rebuilt if the original is newer. --shared-dir names another folder for
them; --no-index-gtf uses the GTF as given. The session folder holds only small files (session.xml,
genome.json, chrom_alias.tab, notation_check.tsv).
Aliases come from --alias-file (e.g. the crossrun bundle's work/annot/chrom_alias.tsv), the FASTA header
descriptions ("... chromosome 1, ..."), a built-in T2T-CHM13v2.0 table (RefSeq NC_0609xx / GenBank CP0682xx)
and simple rules (1 -> chr1, MT -> chrM).

Usage (explicit files):
  make_igv_session.py --out igv/HSV_24h --fasta genome.fa --gtf genes.gtf \\
      --bam S.mLb.clN.sorted.bam --peaks S_peaks.narrowPeak --censat censat.bb --repeatmasker rmsk.bb
Usage (fill everything from a crossrun run; explicit options override):
  make_igv_session.py --out igv/HSV_24h --crossrun-config config.yaml --sample HSV_THP1_WT_24h_REP1 --gtf genes.gtf
Then in IGV Desktop: File > Open Session... and choose <out>/session.xml.

Exit status: 0 = session written; 2 = notation mismatch (nothing written; rerun with --resolve); 1 = other error.
"""
import argparse
import collections
import csv
import gzip
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET

# ----------------------------------------------------------------------------------------- name readers

def _open_text(path):
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def fasta_names(fasta):
    fai = fasta + ".fai"
    if not os.path.exists(fai):
        raise SystemExit(f"FASTA index not found: {fai} (create it with: samtools faidx {fasta})")
    out = collections.OrderedDict()
    for line in open(fai):
        f = line.rstrip("\n").split("\t")
        out[f[0]] = int(f[1])
    return out


def fasta_descriptions(fasta):
    """{name: header description} read by seeking just before each sequence (cheap even for 3 Gb)."""
    out = {}
    with open(fasta + ".fai") as fai, open(fasta, "rb") as fa:
        for line in fai:
            name, _, off = line.split("\t")[:3]
            off = int(off)
            fa.seek(max(0, off - 2000))
            hl = fa.read(off - max(0, off - 2000)).decode(errors="replace").rstrip("\n").split("\n")[-1]
            parts = hl[1:].split(None, 1) if hl.startswith(">") else []
            out[name] = parts[1] if len(parts) > 1 else ""
    return out


def bam_names(bam):
    """Reference names from the BAM header (BGZF is readable with gzip)."""
    with gzip.open(bam, "rb") as fh:
        if fh.read(4) != b"BAM\x01":
            raise SystemExit(f"not a BAM file: {bam}")
        l_text = struct.unpack("<i", fh.read(4))[0]
        fh.read(l_text)
        n_ref = struct.unpack("<i", fh.read(4))[0]
        out = collections.OrderedDict()
        for _ in range(n_ref):
            l_name = struct.unpack("<i", fh.read(4))[0]
            name = fh.read(l_name).rstrip(b"\x00").decode()
            out[name] = struct.unpack("<i", fh.read(4))[0]
    return out


def bigbed_names(path):
    """Chromosome names from the B+ tree index of a bigBed (or bigWig) file."""
    with open(path, "rb") as f:
        h = f.read(64)
        if struct.unpack("<I", h[:4])[0] in (0x8789F2EB, 0x888FFC26):
            e = "<"
        elif struct.unpack(">I", h[:4])[0] in (0x8789F2EB, 0x888FFC26):
            e = ">"
        else:
            raise SystemExit(f"not a bigBed/bigWig file: {path}")
        tree = struct.unpack(e + "Q", h[8:16])[0]
        f.seek(tree)
        magic, _block, key_size, _val, _n = struct.unpack(e + "IIIIQ", f.read(24))
        if magic != 0x78CA8C91:
            raise SystemExit(f"unexpected chromosome index in {path}")
        out = collections.OrderedDict()

        def node(off):
            f.seek(off)
            is_leaf, _r, count = struct.unpack(e + "BBH", f.read(4))
            for _ in range(count):
                key = f.read(key_size).rstrip(b"\x00").decode()
                if is_leaf:
                    _cid, size = struct.unpack(e + "II", f.read(8))
                    out[key] = size
                else:
                    child = struct.unpack(e + "Q", f.read(8))[0]
                    here = f.tell()
                    node(child)
                    f.seek(here)
        node(tree + 32)
    return out


def column1_names(path):
    """Distinct values of column 1 (GTF, narrowPeak, BED); comment/track lines skipped."""
    out = collections.OrderedDict()
    with _open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c = line.split("\t", 1)[0]
            out[c] = out.get(c, 0) + 1
    return out


def is_bgzf(path):
    """True if the file is BGZF (block gzip, as written by bgzip), the only gzip variant tabix can index."""
    with open(path, "rb") as fh:
        h = fh.read(18)
    return len(h) >= 18 and h[:4] == b"\x1f\x8b\x08\x04" and h[12:14] == b"BC"


def tabix_index(path):
    """Path of an existing tabix (.tbi) or CSI (.csi) index for a bgzipped file, else None."""
    if not path.endswith(".gz") or not is_bgzf(path):
        return None
    for ext in (".tbi", ".csi"):
        if os.path.exists(path + ext):
            return path + ext
    return None


def index_names(idx):
    """Sequence names stored in a tabix .tbi or tabix-style .csi index (no need to read the data file)."""
    with gzip.open(idx, "rb") as fh:
        magic = fh.read(4)
        if magic == b"TBI\x01":
            n_ref, _fmt, _cs, _cb, _ce, _meta, _skip, l_nm = struct.unpack("<8i", fh.read(32))
            names = fh.read(l_nm)
        elif magic == b"CSI\x01":
            _min_shift, _depth, l_aux = struct.unpack("<3i", fh.read(12))
            aux = fh.read(l_aux)
            if l_aux < 28:
                return None                     # CSI without tabix metadata: names not stored
            l_nm = struct.unpack("<i", aux[24:28])[0]
            names = aux[28:28 + l_nm]
        else:
            return None
    return collections.OrderedDict((n.decode(), 0) for n in names.split(b"\x00") if n)


def gtf_names(path):
    """Chromosome names of a GTF: from its index when it has one (instant), else by reading the file."""
    idx = tabix_index(path)
    if idx:
        nm = index_names(idx)
        if nm:
            return nm
    return column1_names(path)


# ----------------------------------------------------------------------------------------- notation

def notation(name):
    if name.startswith("chr"):
        return "UCSC (chr1, chr2, ...)"
    if re.match(r"^NC_\d+(\.\d+)?$", name):
        return "RefSeq (NC_...)"
    if re.match(r"^(CP|CM)\d+(\.\d+)?$", name):
        return "GenBank (CP.../CM...)"
    if re.match(r"^(\d{1,2}|X|Y|MT|M)$", name):
        return "Ensembl/NCBI (1, 2, ..., X, MT)"
    return "other"


def t2t_table():
    m = {}
    for n in range(1, 23):
        m[f"NC_{60924 + n:06d}.1"] = f"chr{n}"
        m[f"CP{68278 - n:06d}.2"] = f"chr{n}"
    m.update({"NC_060947.1": "chrX", "CP068255.2": "chrX", "NC_060948.1": "chrY", "CP086569.2": "chrY",
              "NC_012920.1": "chrM", "CP068254.1": "chrM"})
    for k, v in list(m.items()):
        m.setdefault(k.split(".")[0], v)
    return m


def generic_chr(name):
    if name.startswith("chr"):
        return name
    if re.fullmatch(r"\d{1,2}|X|Y", name):
        return "chr" + name
    if name in ("MT", "M"):
        return "chrM"
    return None


def header_chr(desc):
    d = desc.lower()
    if "mitochondri" in d:
        return "chrM"
    if any(w in d for w in ("unplaced", "unlocalized", "scaffold")):
        return None
    m = re.search(r"chromosome\s+([0-9]{1,2}|x|y)\b", d)
    return f"chr{m.group(1).upper()}" if m else None


def build_aliases(fasta, ref, alias_file):
    """Return {fasta_name: chr_name} and {any_name: fasta_name}."""
    user = {}
    if alias_file:
        for line in open(alias_file):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] and f[1]:
                user[f[0]] = f[1]      # alias -> canonical (crossrun chrom_alias.tsv format)
    desc, t2t = fasta_descriptions(fasta), t2t_table()
    to_chr = {}
    for n in ref:
        to_chr[n] = user.get(n) or header_chr(desc.get(n, "")) or t2t.get(n) or generic_chr(n) or n
    to_fasta = {}
    for n, c in to_chr.items():
        for a in {n, c, c[3:] if c.startswith("chr") else None, "MT" if c == "chrM" else None}:
            if a:
                to_fasta.setdefault(a, n)
    for src in (user, t2t):
        for a, c in src.items():
            for n, cc in to_chr.items():
                if cc == c:
                    to_fasta.setdefault(a, n)
    return to_chr, to_fasta


# ----------------------------------------------------------------------------------------- helpers

def parse_prefixes(specs):
    """'OLD=NEW' pairs -> [(old, new)], trailing slashes removed, most specific (longest) OLD first."""
    out = []
    for x in specs:
        if "=" not in x:
            raise SystemExit(f"--rewrite-prefix needs OLD=NEW, got: {x}")
        old, new = x.split("=", 1)
        old, new = old.rstrip("/") or "/", new.rstrip("/")
        out.append((old, new))
    return sorted(out, key=lambda t: -len(t[0]))


def rewrite(path, prefixes, unmatched=None):
    """Rewrite the first matching prefix, matching whole folder names only
    (/home/u/isilon matches /home/u/isilon/x but not /home/u/isilon2/x)."""
    p = os.path.abspath(path)
    for old, new in prefixes:
        if p == old or p.startswith(old + "/"):
            return new + p[len(old):]
    if prefixes and unmatched is not None:
        unmatched.add(p)
    return p


def write_chr_fasta(fasta, to_chr, out_fa):
    """Copy FASTA with chr-style headers and write its .fai (line layout preserved)."""
    fai_rows, off = [], 0
    with open(fasta, "rb") as src, open(out_fa, "wb") as dst:
        name = seqlen = line_bases = line_width = start = None
        for raw in src:
            if raw.startswith(b">"):
                if name is not None:
                    fai_rows.append((name, seqlen, start, line_bases, line_width))
                old = raw[1:].split()[0].decode()
                name = to_chr.get(old, old)
                hdr = f">{name} {old}\n".encode()
                dst.write(hdr); off += len(hdr)
                seqlen, line_bases, line_width, start = 0, None, None, off
            else:
                dst.write(raw); off += len(raw)
                b = len(raw.rstrip(b"\r\n"))
                if line_bases is None:
                    line_bases, line_width = b, len(raw)
                seqlen += b
        if name is not None:
            fai_rows.append((name, seqlen, start, line_bases, line_width))
    with open(out_fa + ".fai", "w") as fh:
        for r in fai_rows:
            fh.write("\t".join(map(str, r)) + "\n")


def strongest_peak(peaks):
    best = None
    with _open_text(peaks) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or not f[1].isdigit():
                continue
            try:
                sig = float(f[6])
            except ValueError:
                continue
            if best is None or sig > best[0]:
                best = (sig, f[0], int(f[1]), int(f[2]))
    return best


def shared_dir_for(src, shared_dir, outdir):
    """Folder for a derived file: --shared-dir if given, else the folder of the original file, so that every
    session points at one central copy. Falls back to the session folder if that folder is not writable."""
    d = os.path.abspath(shared_dir) if shared_dir else os.path.dirname(os.path.abspath(src))
    if shared_dir:
        os.makedirs(d, exist_ok=True)
    if os.access(d, os.W_OK):
        return d
    print(f"[igv] WARNING: cannot write to {d}; the derived copy of {os.path.basename(src)} goes into the session "
          f"folder instead (use --shared-dir to name a writable central folder)", file=sys.stderr)
    return os.path.abspath(outdir)


def current(derived, src, index_ext):
    """True if the derived file and its index exist and neither is older than the original."""
    idx = derived + index_ext
    return (os.path.exists(derived) and os.path.exists(idx)
            and os.path.getmtime(derived) >= os.path.getmtime(src) and os.path.getmtime(idx) >= os.path.getmtime(src))


def adopt(legacy, dest, index_ext, src):
    """Move a copy written into a session folder by an earlier version to the shared location (saves rebuilding)."""
    if os.path.abspath(legacy) == os.path.abspath(dest) or not current(legacy, src, index_ext):
        return False
    shutil.move(legacy, dest)
    shutil.move(legacy + index_ext, dest + index_ext)
    print(f"[igv] moved {legacy} (written by an earlier version) to {dest}", file=sys.stderr)
    return True


def gtf_names_for(gtf):
    """(indexed copy's file name, legacy names) for a GTF/GFF. Only the prefix changes: IGV picks the format from
    the extension left after dropping ".gz", so the name must still end in .gtf.gz / .gff.gz / .gff3.gz."""
    name = os.path.basename(gtf)
    m = re.search(r"\.(gtf|gff3|gff)(\.gz)?$", name, flags=re.I)
    if m:
        core = name[:m.start()] + "." + m.group(1)
        legacy = [name[:m.start()] + ".sorted." + m.group(1).lower() + ".gz"]
    else:                                   # unrecognised name: the copy has to gain ".gtf" for IGV to read it
        core = name + ".gtf"
        legacy = [name + ".sorted.gtf.gz"]
    legacy.append(re.sub(r"\.gz$", "", name) + ".sorted.gz")
    return "sorted_" + core + ".gz", legacy


def index_gtf(gtf, outdir, shared_dir=None):
    """Sorted, bgzipped, tabix-indexed copy of the GTF (fast in IGV), written once beside the original as
    sorted_<name>.gtf.gz and reused by every later session. A GTF that is already bgzipped and indexed is used as-is."""
    if tabix_index(gtf):
        print(f"[igv] GTF is already indexed ({os.path.basename(tabix_index(gtf))}); using it as-is", file=sys.stderr)
        return gtf
    fname, legacy = gtf_names_for(gtf)
    d = shared_dir_for(gtf, shared_dir, outdir)
    dest = os.path.join(d, fname)
    if current(dest, gtf, ".tbi"):
        print(f"[igv] reusing the indexed GTF {dest}", file=sys.stderr)
        return dest
    for old in legacy:
        if adopt(os.path.join(outdir, old), dest, ".tbi", gtf):
            return dest
    if not (shutil.which("bgzip") and shutil.which("tabix")):
        print("[igv] bgzip/tabix not found: the GTF is used as-is (IGV reads it all at load time)", file=sys.stderr)
        return gtf
    tmp = os.path.join(d, f"tmp{os.getpid()}_{fname}")          # finished files appear under their final name only
    cat = "gzip -dc" if gtf.endswith(".gz") else "cat"
    cmd = (f"set -o pipefail; {cat} {shlexq(gtf)} | grep -v '^#' | LC_ALL=C sort -t $'\\t' -k1,1 -k4,4n | bgzip > {shlexq(tmp)} "
           f"&& tabix -p gff {shlexq(tmp)}")
    print(f"[igv] indexing the GTF (one-off): {dest}", file=sys.stderr)
    try:
        subprocess.run(["bash", "-c", cmd], check=True)
        os.replace(tmp, dest)
        os.replace(tmp + ".tbi", dest + ".tbi")
    finally:
        for f in (tmp, tmp + ".tbi"):
            if os.path.exists(f):
                os.remove(f)
    return dest


def chr_fasta(fasta, to_chr, outdir, shared_dir=None):
    """chr-named copy of the FASTA, written once beside the original as chr_<name> (+ .fai) and reused."""
    d = shared_dir_for(fasta, shared_dir, outdir)
    dest = os.path.join(d, "chr_" + os.path.basename(fasta))
    if current(dest, fasta, ".fai"):
        print(f"[igv] reusing the chr-named FASTA {dest}", file=sys.stderr)
        return dest
    if adopt(os.path.join(outdir, "genome.chr.fa"), dest, ".fai", fasta):
        return dest
    tmp = os.path.join(d, f"tmp{os.getpid()}_chr_" + os.path.basename(fasta))
    print(f"[igv] writing a chr-named copy of the FASTA (one-off): {dest}", file=sys.stderr)
    try:
        write_chr_fasta(fasta, to_chr, tmp)
        os.replace(tmp, dest)
        os.replace(tmp + ".fai", dest + ".fai")
    finally:
        for f in (tmp, tmp + ".fai"):
            if os.path.exists(f):
                os.remove(f)
    return dest


def shlexq(s):
    return "'" + s.replace("'", "'\\''") + "'"


# ----------------------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="output directory (session.xml and helper files go here)")
    ap.add_argument("--fasta"); ap.add_argument("--gtf"); ap.add_argument("--bam"); ap.add_argument("--peaks")
    ap.add_argument("--censat"); ap.add_argument("--repeatmasker")
    ap.add_argument("--crossrun-config", help="crossrun config.yaml: fills fasta, bigBeds, BAM and peaks for --sample")
    ap.add_argument("--sample", help="sample_id from the crossrun manifest (with --crossrun-config)")
    ap.add_argument("--resolve", choices=["none", "alias", "chr-fasta"], default="none",
                    help="what to do if the notation differs (default: stop and report)")
    ap.add_argument("--alias-file", help="extra aliases, 'alias<TAB>canonical' per line (e.g. crossrun work/annot/chrom_alias.tsv)")
    ap.add_argument("--locus", help="initial locus (default: the sample's strongest peak +/- 5 kb)")
    ap.add_argument("--no-index-gtf", action="store_true",
                    help="use the GTF as given instead of a sorted, bgzipped, tabix-indexed copy")
    ap.add_argument("--index-gtf", action="store_true", help=argparse.SUPPRESS)      # now the default; still accepted
    ap.add_argument("--shared-dir",
                    help="folder for the derived files (chr_<fasta>, sorted_<gtf>.gz and their indexes) "
                         "instead of the folders of the original FASTA and GTF")
    ap.add_argument("--rewrite-prefix", action="append", default=[], metavar="OLD=NEW",
                    help="rewrite path prefixes in the session (e.g. cluster path -> path on the computer running IGV)")
    ap.add_argument("--min-match", type=float, default=0.9,
                    help="a file matches the genome if at least this fraction of its chromosome names do (default 0.9)")
    a = ap.parse_args()

    # ---- fill from a crossrun run
    if a.crossrun_config:
        import yaml
        cfg = yaml.safe_load(open(a.crossrun_config))
        a.fasta = a.fasta or cfg.get("genome_fasta")
        a.censat = a.censat or cfg.get("censat")
        a.repeatmasker = a.repeatmasker or cfg.get("repeatmasker")
        work = os.path.join(os.path.abspath(cfg.get("outdir", "crossrun_out")), "work")
        if not a.alias_file and os.path.exists(os.path.join(work, "annot", "chrom_alias.tsv")):
            a.alias_file = os.path.join(work, "annot", "chrom_alias.tsv")
        if a.sample:
            rows = list(csv.DictReader(open(os.path.join(work, "samples.tsv")), delimiter="\t"))
            row = next((r for r in rows if r["sample_id"] == a.sample), None)
            if row is None:
                raise SystemExit(f"sample '{a.sample}' not in {work}/samples.tsv; available: {', '.join(r['sample_id'] for r in rows[:10])} ...")
            a.bam = a.bam or row.get("bam_file")
            a.peaks = a.peaks or row.get("peaks_file")
    tracks = [("Genome annotation (GTF)", "gtf", a.gtf), ("CenSat", "censat", a.censat),
              ("RepeatMasker", "repeatmasker", a.repeatmasker), ("Peaks", "peaks", a.peaks), ("Reads", "bam", a.bam)]
    missing = [k for _, k, p in tracks if not p] + ([] if a.fasta else ["fasta"])
    if missing:
        raise SystemExit(f"missing inputs: {', '.join('--' + m for m in missing)}")
    for _, k, p in tracks + [("Genome", "fasta", a.fasta)]:
        if not os.path.exists(p):
            raise SystemExit(f"--{k}: file not found: {p}")
    if not re.search(r"\.(gtf|gff3|gff)(\.gz)?$", a.gtf, flags=re.I):
        print(f"[igv] WARNING: {os.path.basename(a.gtf)} does not end in .gtf/.gff/.gff3 (optionally .gz); IGV chooses "
              f"the format from the file name and may not show it"
              + ("; rename it" if a.no_index_gtf else "; the indexed copy gets a name ending in .gtf.gz"),
              file=sys.stderr)
    gidx = tabix_index(a.gtf)
    if gidx:
        if os.path.getmtime(gidx) < os.path.getmtime(a.gtf):
            print(f"[igv] WARNING: {os.path.basename(gidx)} is older than the GTF; if the GTF changed, re-index it "
                  f"(tabix -p gff {a.gtf}) or IGV may show wrong features", file=sys.stderr)
    elif a.gtf.endswith(".gz") and not is_bgzf(a.gtf):
        print("[igv] note: the GTF is gzip- (not bgzip-) compressed, so it cannot be indexed in place; "
              "an indexed copy is written instead", file=sys.stderr)
    if not (os.path.exists(a.bam + ".bai") or os.path.exists(re.sub(r"\.bam$", ".bai", a.bam))):
        raise SystemExit(f"BAM index not found for {a.bam} (samtools index it); IGV needs it")

    # ---- chromosome names in every file
    ref = fasta_names(a.fasta)
    readers = {"gtf": gtf_names, "censat": bigbed_names, "repeatmasker": bigbed_names,
               "peaks": column1_names, "bam": bam_names}
    names = {"fasta": ref}
    for _, k, p in tracks:
        if k in ("censat", "repeatmasker") and not p.lower().endswith((".bb", ".bigbed")):
            names[k] = column1_names(p)      # plain BED given instead of bigBed
        else:
            names[k] = readers[k](p)
    to_chr, to_fasta = build_aliases(a.fasta, ref, a.alias_file)

    print(f"\nChromosome notation check (reference: {a.fasta})\n", file=sys.stderr)
    print(f"  {'file':<14}{'notation':<34}{'names':>7}{'exact':>8}{'via alias':>11}  examples", file=sys.stderr)
    report, ok_exact, ok_alias, all_chr = [], True, True, True
    for k, nm in names.items():
        styles = collections.Counter(notation(n) for n in nm)
        main_style = styles.most_common(1)[0][0] if styles else "none"
        exact = sum(1 for n in nm if n in ref)
        alias = sum(1 for n in nm if n in ref or n in to_fasta)
        frac_e = exact / max(1, len(nm)); frac_a = alias / max(1, len(nm))
        if k != "fasta":
            ok_exact &= frac_e >= a.min_match
            ok_alias &= frac_a >= a.min_match
        all_chr &= main_style.startswith("UCSC")
        ex = ", ".join(list(nm)[:3])
        print(f"  {k:<14}{main_style:<34}{len(nm):>7}{exact:>8}{alias:>11}  {ex}", file=sys.stderr)
        unmatched = [n for n in nm if n not in ref and n not in to_fasta]
        if unmatched and k != "fasta":
            print(f"  {'':<14}unmatched: {', '.join(unmatched[:5])}{' ...' if len(unmatched) > 5 else ''}", file=sys.stderr)
        report.append((k, main_style, len(nm), exact, alias))
    print("", file=sys.stderr)

    consistent = ok_exact and all_chr
    if consistent:
        mode = "direct"
        print("[igv] all files use chr-style names that match the genome: writing the session directly.", file=sys.stderr)
    elif a.resolve == "none":
        why = ("files use different chromosome notations" if not ok_exact
               else "all files match each other but do not use chr-style notation")
        print(f"[igv] STOP: {why}.", file=sys.stderr)
        if ok_alias:
            print("  Every file can be matched through chromosome aliases. Rerun with one of:\n"
                  "    --resolve chr-fasta  (IGV shows chr1, chr2, ...; writes a header-renamed copy of the FASTA)\n"
                  "    --resolve alias      (no copies; IGV shows the FASTA's own names)", file=sys.stderr)
        else:
            print("  Some files cannot be matched even through aliases (see 'unmatched' above); supply\n"
                  "  --alias-file with a mapping, or check that the files belong to the same assembly.", file=sys.stderr)
        sys.exit(2)
    else:
        if not ok_alias:
            raise SystemExit("[igv] some files cannot be matched even through aliases (see 'unmatched' above); "
                             "supply --alias-file or check the files belong to the same assembly")
        mode = a.resolve

    # ---- write outputs
    os.makedirs(a.out, exist_ok=True)
    prefixes = parse_prefixes(a.rewrite_prefix)
    not_rewritten = set()
    P = lambda p: rewrite(p, prefixes, not_rewritten)
    gtf = a.gtf if a.no_index_gtf else index_gtf(a.gtf, a.out, a.shared_dir)

    if mode == "direct":
        genome_ref, display = P(a.fasta), (lambda c: c)
    else:
        if mode == "chr-fasta":
            fa = chr_fasta(a.fasta, to_chr, a.out, a.shared_dir)
            display = lambda c: to_chr.get(to_fasta.get(c, c), c)
            primary = {n: to_chr[n] for n in ref}          # FASTA name -> displayed (chr) name
        else:
            fa = a.fasta
            display = lambda c: to_fasta.get(c, c)
            primary = {n: n for n in ref}
        # alias file: one line per chromosome, displayed name first, then every other known name
        names_for = collections.defaultdict(set)
        for alias_name, fname in to_fasta.items():
            names_for[fname].add(alias_name)
        with open(os.path.join(a.out, "chrom_alias.tab"), "w") as fh:
            for n in ref:
                shown = primary[n]
                others = sorted((names_for[n] | {n, to_chr[n]}) - {shown})
                fh.write("\t".join([shown] + others) + "\n")
        genome = {"id": os.path.basename(a.fasta), "name": os.path.basename(a.fasta) + " (chromosome aliases)",
                  "fastaURL": P(fa), "indexURL": P(fa + ".fai"), "aliasURL": P(os.path.join(a.out, "chrom_alias.tab"))}
        with open(os.path.join(a.out, "genome.json"), "w") as fh:
            json.dump(genome, fh, indent=2)
        genome_ref = P(os.path.join(a.out, "genome.json"))

    if a.locus:
        locus = a.locus
    else:
        best = strongest_peak(a.peaks)
        locus = (f"{display(best[1])}:{max(1, best[2] - 5000)}-{best[3] + 5000}" if best else display(next(iter(ref))))

    sess = ET.Element("Session", genome=genome_ref, locus=locus, version="8", hasGeneTrack="false", hasSequenceTrack="true")
    res = ET.SubElement(sess, "Resources")
    sample = a.sample or re.sub(r"\.mLb\.clN\.sorted\.bam$|\.bam$", "", os.path.basename(a.bam))
    labels = {"gtf": "Genes (GTF)", "censat": "CenSat", "repeatmasker": "RepeatMasker",
              "peaks": f"{sample} peaks", "bam": f"{sample} reads"}
    paths = {"gtf": gtf, "censat": a.censat, "repeatmasker": a.repeatmasker, "peaks": a.peaks, "bam": a.bam}
    for k in ("gtf", "censat", "repeatmasker", "peaks", "bam"):
        ET.SubElement(res, "Resource", name=labels[k], path=P(paths[k]))
    ET.indent(sess, space="    ")
    out = os.path.join(a.out, "session.xml")
    ET.ElementTree(sess).write(out, encoding="UTF-8", xml_declaration=True)
    with open(os.path.join(a.out, "notation_check.tsv"), "w") as fh:
        fh.write("file\tnotation\tn_names\tmatch_exact\tmatch_with_aliases\n")
        for r in report:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"[igv] wrote {out} (mode: {mode}; locus {locus}). In IGV Desktop: File > Open Session...", file=sys.stderr)
    if not_rewritten:
        print(f"[igv] WARNING: {len(not_rewritten)} path(s) in the session match no --rewrite-prefix and keep their "
              "current location; IGV on the other computer will not find them unless that location is visible there:",
              file=sys.stderr)
        for q in sorted(not_rewritten):
            print(f"    {q}", file=sys.stderr)


if __name__ == "__main__":
    main()
