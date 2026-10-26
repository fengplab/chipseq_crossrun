#!/usr/bin/env python3
"""Supplementary QC table: one row per library (each chosen IP sample and every input it uses), with sequencing quality,
alignment quality, library composition and peak quality, collected from the ChIP-seq pipeline's QC files and the
crossrun tables. No thresholds, no flags: values only.

Usage:
  python3 extra/compile_qc_table.py --crossrun <crossrun_out> --labels-file fig1_heatmap_labels.tsv [--out supp_qc_table]

  --crossrun DIR     crossrun output folder (reads work/samples.tsv and tables/)
  --labels-file TSV  sample_id, label (tab-separated, with header): the IP samples, in row order, and their labels
  --out PREFIX       writes PREFIX.tsv, PREFIX.xlsx (table + data dictionary; needs openpyxl, else skipped),
                     PREFIX_sources.tsv (the file every value came from) and PREFIX_missing.tsv (values not found)

Columns (blank = not found or not applicable; peak columns are blank for inputs):
  Sequencing (FastQC on raw read 1, Trim Galore report of read 1; summed over sequencing runs _T1, _T2, ...)
    raw_reads, read_length, pct_gc, pct_fastqc_duplicate, pct_reads_with_adapter, pct_reads_passing_trimming
  Alignment
    pct_aligned (Bowtie2 overall alignment rate), pct_properly_paired and mapped_reads (samtools flagstat before
    filtering, *.mLb.mkD*), pct_duplicates (Picard MarkDuplicates PERCENT_DUPLICATION x 100), reads_after_filtering
    (samtools flagstat of the filtered BAM, *.mLb.clN*), fingerprint_elbow_point and fingerprint_js_distance
    (deepTools plotFingerprint quality metrics)
  Composition (crossrun table00c: shares of all reads)
    pct_chrM, pct_viral, pct_unaligned
  Peaks (crossrun table00 and qc01b; IP only)
    n_peaks, median_peak_width, median_fold_enrichment, pct_peaks_q_below_0.01, pct_peaks_in_censat,
    median_peak_mappability
  Fragments (crossrun qc03)
    median_insert_outside_censat

The pipeline file names follow the nf-core/chipseq layout under each sample's results_dir and are found by pattern
(see PATTERNS). UNVERIFIED: the patterns and the text formats parsed below were written for nf-core/chipseq 2.x
outputs and checked only against synthetic test files; check PREFIX_missing.tsv and PREFIX_sources.tsv on the first
real run.
"""
import argparse, csv, glob, io, os, re, sys, zipfile
from collections import OrderedDict, defaultdict

def die(m): sys.exit(f"[qc-table] ERROR: {m}")
def say(m): print(f"[qc-table] {m}", file=sys.stderr)

# UNVERIFIED: nf-core/chipseq 2.x output layout; {lib} = library id, {al} = aligner directory (e.g. bowtie2)
PATTERNS = {
    "fastqc":      ["fastqc/**/{lib}_T*_1_fastqc.zip", "fastqc/**/{lib}_T*_fastqc.zip", "fastqc/**/{lib}_1_fastqc.zip",
                    "fastqc/**/{lib}_fastqc.zip"],
    "trimgalore":  ["trimgalore/**/{lib}_T*_1.fastq.gz_trimming_report.txt", "trimgalore/**/{lib}_T*.fastq.gz_trimming_report.txt",
                    "trimgalore/**/{lib}_1.fastq.gz_trimming_report.txt", "trimgalore/**/{lib}.fastq.gz_trimming_report.txt"],
    "bowtie2_log": ["{al}/**/{lib}_T*.bowtie2.log", "{al}/**/{lib}.bowtie2.log"],
    "picard":      ["{al}/merged_library/**/{lib}.mLb.mkD*MarkDuplicates.metrics.txt"],
    "flagstat_mkd": ["{al}/merged_library/**/{lib}.mLb.mkD*.flagstat"],
    "flagstat_cln": ["{al}/merged_library/**/{lib}.mLb.clN*.flagstat"],
    "fingerprint": ["{al}/merged_library/**/*.plotFingerprint.qcmetrics.txt"],
}

COLS = OrderedDict([
    ("sample_id", "Library id (pipeline sample id)"), ("label", "Label used in the figures"),
    ("library_type", "IP or input"), ("matched_input", "IP: its input library; input: blank"),
    ("input_for", "Input: the labelled IP samples that use it; IP: blank"),
    ("cell_type", "Cell line"), ("virus", "Virus"), ("hpi", "Hours post-infection"), ("cgas_construct", "cGAS construct"),
    ("run_id", "Sequencing run / pipeline run"),
    ("raw_reads", "Read 1 sequences before trimming (FastQC Total Sequences), summed over sequencing runs"),
    ("read_length", "Read length before trimming (FastQC Sequence length)"),
    ("pct_gc", "% GC of read 1 (FastQC), weighted by reads over sequencing runs"),
    ("pct_fastqc_duplicate", "100 - FastQC Total Deduplicated Percentage (sequence-level duplication of read 1)"),
    ("pct_reads_with_adapter", "% of read 1 with adapter trimmed (Trim Galore report)"),
    ("pct_reads_passing_trimming", "% of read 1 written after trimming (Trim Galore report)"),
    ("pct_aligned", "Bowtie2 overall alignment rate, pooled over sequencing runs"),
    ("mapped_reads", "Mapped reads before filtering (samtools flagstat, duplicate-marked merged library)"),
    ("pct_properly_paired", "% properly paired before filtering (samtools flagstat, duplicate-marked merged library)"),
    ("pct_duplicates", "Picard MarkDuplicates PERCENT_DUPLICATION x 100 (merged library)"),
    ("reads_after_filtering", "Mapped reads in the filtered BAM used for peak calling (samtools flagstat)"),
    ("fingerprint_elbow_point", "deepTools plotFingerprint elbow point (lower = more enrichment)"),
    ("fingerprint_js_distance", "deepTools plotFingerprint Jensen-Shannon distance (blank if not computed)"),
    ("pct_chrM", "% of all reads on chrM (crossrun table00c)"),
    ("pct_viral", "% of all reads on viral genomes (crossrun table00c)"),
    ("pct_unaligned", "% of all reads unaligned (crossrun table00c)"),
    ("n_peaks", "Number of peaks (crossrun table00)"),
    ("median_peak_width", "Median peak width, bp (crossrun table00)"),
    ("median_fold_enrichment", "Median MACS fold enrichment of peaks (crossrun table00)"),
    ("pct_peaks_q_below_0.01", "% of peaks with q < 0.01 (crossrun table00)"),
    ("pct_peaks_in_censat", "% of peaks in CenSat arrays (crossrun qc01b)"),
    ("median_peak_mappability", "Median k-mer mappability of peaks (crossrun qc01b)"),
    ("median_insert_outside_censat", "Median fragment length outside CenSat arrays, bp (crossrun qc03)"),
])

def read_tsv(f):
    with open(f, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def num(s):
    s = str(s).strip().replace(",", "").rstrip("%")
    try: return float(s)
    except ValueError: return None

def find(res, al, key, lib):
    hits = []
    for p in PATTERNS[key]:
        hits = sorted(glob.glob(os.path.join(glob.escape(res), p.format(lib=glob.escape(lib), al=al)), recursive=True))
        if hits: break
    return hits

def lib_name(res, al, lib):
    """The library id as the pipeline wrote it: lib itself, or lib_REP<n> (manifests may omit the replicate)."""
    for cand in [lib] + [f"{lib}_REP{i}" for i in range(1, 10)]:
        if find(res, al, "flagstat_cln", cand) or find(res, al, "fastqc", cand): return cand
    return lib

def parse_fastqc(files):
    reads = gc_w = dd_w = 0.0; length = set()
    for f in files:
        with zipfile.ZipFile(f) as z:
            name = next((n for n in z.namelist() if n.endswith("fastqc_data.txt")), None)
            if not name: continue
            kv = {}
            for line in io.TextIOWrapper(z.open(name), encoding="utf-8"):
                p = line.rstrip("\n").lstrip("#").split("\t")
                if len(p) >= 2: kv.setdefault(p[0].strip(), p[1].strip())
        n = num(kv.get("Total Sequences", ""))
        if n is None: continue
        reads += n; length.add(kv.get("Sequence length", ""))
        if num(kv.get("%GC", "")) is not None: gc_w += n * num(kv["%GC"])
        if num(kv.get("Total Deduplicated Percentage", "")) is not None: dd_w += n * num(kv["Total Deduplicated Percentage"])
    if not reads: return {}
    return {"raw_reads": int(reads), "read_length": "/".join(sorted(x for x in length if x)),
            "pct_gc": round(gc_w / reads, 1), "pct_fastqc_duplicate": round(100 - dd_w / reads, 2) if dd_w else ""}

def parse_trimgalore(files):
    tot = adp = kept = 0.0
    for f in files:
        t = open(f).read()
        m1 = re.search(r"Total reads processed:\s+([\d,]+)", t); m2 = re.search(r"Reads with adapters:\s+([\d,]+)", t)
        m3 = re.search(r"Reads written \(passing filters\):\s+([\d,]+)", t)
        if not (m1 and m2 and m3): continue
        tot += num(m1.group(1)); adp += num(m2.group(1)); kept += num(m3.group(1))
    if not tot: return {}
    return {"pct_reads_with_adapter": round(100 * adp / tot, 2), "pct_reads_passing_trimming": round(100 * kept / tot, 2)}

def parse_bowtie2(files):
    tot = ali = 0.0
    for f in files:
        t = open(f).read()
        m1 = re.search(r"^(\d+) reads; of these:", t, re.M); m2 = re.search(r"([\d.]+)% overall alignment rate", t)
        if not (m1 and m2): continue
        n = float(m1.group(1)); tot += n; ali += n * float(m2.group(1)) / 100
    return {"pct_aligned": round(100 * ali / tot, 2)} if tot else {}

def parse_flagstat(f):
    out = {}
    for line in open(f):
        m = re.match(r"(\d+) \+ (\d+) (in total|mapped|properly paired)\b", line)
        if m: out.setdefault(m.group(3), int(m.group(1)) + int(m.group(2)))
    return out

def parse_picard(f):
    lines = [l.rstrip("\n") for l in open(f)]
    for i, l in enumerate(lines):
        if l.startswith("LIBRARY\t") and i + 1 < len(lines):
            row = dict(zip(l.split("\t"), lines[i + 1].split("\t")))
            v = num(row.get("PERCENT_DUPLICATION", ""))
            return {"pct_duplicates": round(100 * v, 2)} if v is not None else {}
    return {}

def fingerprint_index(results_dirs, al):
    idx = {}
    for res in results_dirs:
        for f in find(res, al, "fingerprint", "*"):
            for r in read_tsv(f):
                lab = (r.get("Sample") or "").strip()
                if lab: idx.setdefault(lab, (r, f))
    return idx

def fp_lookup(idx, lib):
    for lab, (r, f) in idx.items():
        if lab == lib or re.match(re.escape(lib) + r"([._].*)?$", lab):
            ep, js = num(r.get("Elbow Point", "")), num(r.get("JS Distance", ""))
            return {"fingerprint_elbow_point": ep if ep is not None else "", "fingerprint_js_distance": js if js is not None else ""}, f
    return {}, None

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--crossrun", required=True); ap.add_argument("--labels-file", required=True)
    ap.add_argument("--out", default="supp_qc_table")
    a = ap.parse_args()
    sfile = os.path.join(a.crossrun, "work", "samples.tsv")
    if not os.path.exists(sfile): die(f"not found: {sfile} (run step 0 of the crossrun first)")
    S = {r["sample_id"]: r for r in read_tsv(sfile)}
    L = read_tsv(a.labels_file)
    if not L or not {"sample_id", "label"} <= set(L[0]): die(f"{a.labels_file} needs a header with columns sample_id and label")
    ids = [r["sample_id"].strip() for r in L if r["sample_id"].strip()]
    lab = {r["sample_id"].strip(): r["label"].strip() for r in L}
    if len(set(ids)) != len(ids): die("--labels-file repeats a sample_id")
    miss = [s for s in ids if s not in S]
    if miss: die("samples not in work/samples.tsv (check spelling and suffixes such as _REP1): " + ", ".join(miss))

    T = os.path.join(a.crossrun, "tables")
    def tab(n):
        f = os.path.join(T, n + ".tsv")
        if os.path.exists(f): return read_tsv(f), f
        say(f"note: {f} not found; its columns stay blank"); return [], None
    t00, f00 = tab("table00_peak_summary"); t00c, f00c = tab("table00c_contig_reads")
    q01b, fq01b = tab("qc01b_peak_mappability_by_sample"); q03, fq03 = tab("qc03_insert_size_stats")

    rows, sources, missing = [], [], []
    al_of = lambda s: S[s].get("aligner_resolved") or "bowtie2"
    fp_idx = fingerprint_index(sorted({S[s]["results_dir"] for s in ids}), "*")

    def pipeline_cols(row, res, al, lib):
        for key, fn in (("fastqc", parse_fastqc), ("trimgalore", parse_trimgalore), ("bowtie2_log", parse_bowtie2)):
            fs = find(res, al, key, lib); v = fn(fs) if fs else {}
            row.update(v)
            for c in v: sources.append((lib, c, ";".join(fs)))
            if not v: missing.append((lib, key, " | ".join(p.format(lib=lib, al=al) for p in PATTERNS[key]), res))
        for key in ("flagstat_mkd", "picard", "flagstat_cln"):
            fs = find(res, al, key, lib)
            if not fs: missing.append((lib, key, " | ".join(p.format(lib=lib, al=al) for p in PATTERNS[key]), res)); continue
            f = fs[0]
            if key == "picard": v = parse_picard(f)
            else:
                fl = parse_flagstat(f)
                if key == "flagstat_mkd":
                    v = {}
                    if "mapped" in fl: v["mapped_reads"] = fl["mapped"]
                    if "properly paired" in fl and fl.get("in total"): v["pct_properly_paired"] = round(100 * fl["properly paired"] / fl["in total"], 2)
                else: v = {"reads_after_filtering": fl["mapped"]} if "mapped" in fl else {}
            row.update(v)
            for c in v: sources.append((lib, c, f))
        v, f = fp_lookup(fp_idx, lib)
        row.update(v)
        for c in v: sources.append((lib, c, f))
        if not v: missing.append((lib, "fingerprint", "plotFingerprint.qcmetrics.txt row for this library", res))

    def composition(row, sid, libtype, lib):
        part = [r for r in t00c if r["sample_id"] == sid and r["library"] == libtype]
        if not part: return
        pct = lambda cond: round(sum(num(r["pct_reads"]) or 0 for r in part if cond(r)), 3)
        row.update({"pct_chrM": pct(lambda r: r["contig"] == "chrM"), "pct_viral": pct(lambda r: r["group"] == "viral"),
                    "pct_unaligned": pct(lambda r: r["group"] == "unaligned")})
        for c in ("pct_chrM", "pct_viral", "pct_unaligned"): sources.append((lib, c, f00c))

    def insert(row, lib, role):
        r = [x for x in q03 if x["library_id"] == lib and x["role"] == role and x["class"] == "outside_censat"]
        if r:
            row["median_insert_outside_censat"] = num(r[0]["median_insert"]); sources.append((lib, "median_insert_outside_censat", fq03))

    inputs = OrderedDict()
    for sid in ids:
        s = S[sid]; res = s["results_dir"]; al = al_of(sid); lib = s.get("pipeline_id") or sid
        ctrl = (s.get("control_id") or "").strip(); clib = lib_name(res, al, ctrl) if ctrl else ""
        row = {"sample_id": lib, "label": lab[sid], "library_type": "IP", "matched_input": clib,
               **{k: s.get(k, "") for k in ("cell_type", "virus", "hpi", "cgas_construct", "run_id")}}
        pipeline_cols(row, res, al, lib)
        composition(row, sid, "IP", lib)
        p = next((r for r in t00 if r["sample_id"] == sid), None)
        if p:
            n = num(p["n_peaks"]) or 0
            row.update({"n_peaks": int(n), "median_peak_width": num(p["width_median"]), "median_fold_enrichment": num(p["fold_enrichment_median"]),
                        "pct_peaks_q_below_0.01": round(100 * (num(p["n_q_below_0.01"]) or 0) / n, 2) if n else ""})
            for c in ("n_peaks", "median_peak_width", "median_fold_enrichment", "pct_peaks_q_below_0.01"): sources.append((lib, c, f00))
        m = next((r for r in q01b if r["sample_id"] == sid), None)
        if m:
            row.update({"pct_peaks_in_censat": round(num(m["pct_peaks_in_censat"]), 2), "median_peak_mappability": num(m["median_mappability"])})
            for c in ("pct_peaks_in_censat", "median_peak_mappability"): sources.append((lib, c, fq01b))
        insert(row, lib, "IP")
        rows.append(row)
        if clib:
            inputs.setdefault(clib, {"res": res, "al": al, "sid": sid, "users": [], "s": s})["users"].append(lab[sid])
        else:
            missing.append((lib, "control_id", "no input listed for this sample in the manifest", res))
    for clib, d in inputs.items():
        s = d["s"]
        row = {"sample_id": clib, "label": f"Input ({s.get('cell_type', '')})".replace(" ()", ""), "library_type": "input",
               "input_for": "; ".join(d["users"]), **{k: s.get(k, "") for k in ("cell_type", "run_id")}}
        pipeline_cols(row, d["res"], d["al"], clib)
        composition(row, d["sid"], "control", clib)
        insert(row, clib, "control")
        rows.append(row)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out + ".tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(COLS), delimiter="\t", extrasaction="ignore"); w.writeheader()
        for r in rows: w.writerow({k: r.get(k, "") for k in COLS})
    with open(a.out + "_sources.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["library", "column", "source_file"]); w.writerows(sources)
    with open(a.out + "_missing.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["library", "item", "looked_for", "results_dir"]); w.writerows(missing)
    try:
        import openpyxl
        from openpyxl.styles import Font
        wb = openpyxl.Workbook(); ws = wb.active; ws.title = "QC"
        ws.append(list(COLS))
        for c in ws[1]: c.font = Font(bold=True)
        for r in rows: ws.append([r.get(k, None) if r.get(k, "") != "" else None for k in COLS])
        ws.freeze_panes = "C2"
        dd = wb.create_sheet("Data dictionary"); dd.append(["column", "meaning"])
        for c in dd[1]: c.font = Font(bold=True)
        for k, v in COLS.items(): dd.append([k, v])
        wb.save(a.out + ".xlsx"); xl = a.out + ".xlsx"
    except ImportError:
        xl = None; say("openpyxl not installed: .xlsx skipped")
    filled = sum(1 for r in rows for k in COLS if r.get(k, "") != "")
    say(f"{len(ids)} IP samples + {len(inputs)} inputs -> {a.out}.tsv{' and ' + xl if xl else ''}; "
        f"{filled} of {len(rows) * len(COLS)} cells filled; {len(missing)} items not found (see {a.out}_missing.tsv)")

if __name__ == "__main__":
    main()
