#!/usr/bin/env python3
"""
Validate config + manifest, resolve every per-sample pipeline output, and write:
  <outdir>/work/config.env    shell-sourceable settings for the bash steps
  <outdir>/work/samples.tsv   manifest + resolved file paths (read by every later step)
  <outdir>/work/missing_files.tsv  report of anything that could not be found, with what was tried
  <outdir>/work/layouts.tsv        aligner/peak-type directories found and used per sample

Usage: 00_setup.py config.yaml
"""
import csv
import glob
import os
import re
import shlex
import sys

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required (pip install pyyaml)")

REQUIRED_COLS = ["sample_id", "run_id", "results_dir"]

# Spreadsheet placeholders for "no value" -> treated as blank everywhere in the manifest.
# ("none" is deliberately NOT included: it is a real level, e.g. virus = none for mock samples.)
NA_TOKENS = {"na", "n/a", "#n/a", "nan", "null", "-", "."}

ALIGNERS = {"bwa", "bowtie2", "chromap", "star"}
PEAK_TYPES = {"narrow_peak": "narrow_peak", "narrow": "narrow_peak", "narrowpeak": "narrow_peak",
              "broad_peak": "broad_peak", "broad": "broad_peak", "broadpeak": "broad_peak"}


def clean(v):
    v = (v or "").strip()
    return "" if v.lower() in NA_TOKENS else v


def read_manifest(path):
    """Read the TAB-separated manifest strictly, line by line, with clear errors.
    Tolerates a byte-order mark, Windows line endings, blank lines and trailing empty fields (stray tabs at
    the end of a line, common in spreadsheet exports). Rows with extra NON-empty fields are an error, because
    their values would not belong to any column; short rows are padded with blanks."""
    with open(path, encoding="utf-8-sig", newline="") as fh:
        lines = [l.rstrip("\r\n") for l in fh]
    lines = [(i + 1, l) for i, l in enumerate(lines) if l.strip() and not l.startswith("#")]
    if not lines:
        die("manifest is empty")
    header = [h.strip() for h in lines[0][1].split("\t")]
    while header and header[-1] == "":
        header.pop()                                   # trailing tabs on the header line
    if "" in header:
        die(f"manifest header has an empty column name (column {header.index('') + 1}); every column needs a name")
    dup = sorted({h for h in header if header.count(h) > 1})
    if dup:
        die(f"manifest header repeats column name(s): {', '.join(dup)}")
    rows, problems = [], []
    for ln, l in lines[1:]:
        f = [clean(x) for x in l.split("\t")]
        while len(f) > len(header) and f[-1] == "":
            f.pop()                                    # trailing empty fields
        if len(f) > len(header):
            extra = f[len(header):]
            problems.append(f"  line {ln} (sample '{f[0]}'): {len(f)} fields but the header has {len(header)}; "
                            f"unexpected extra value(s): {', '.join(repr(x) for x in extra[:5])}")
            continue
        f += [""] * (len(header) - len(f))
        rows.append(dict(zip(header, f)))
    if problems:
        die("manifest rows have more fields than the header has columns. A value may contain a TAB, or the "
            "header may be missing a column name:\n" + "\n".join(problems[:20]))
    return [r for r in rows if r.get("sample_id")]


def norm_aligner(v, sample, source):
    v = clean(v).lower()
    if v in ("", "auto"):
        return "auto"
    if v not in ALIGNERS:
        warn(f"{sample}: {source} aligner '{v}' is not one of {sorted(ALIGNERS)}; detecting automatically")
        return "auto"
    return v


def norm_peak_type(v, sample, source):
    v = clean(v).lower()
    if v in ("", "auto"):
        return "auto"
    if v not in PEAK_TYPES:
        warn(f"{sample}: {source} peak_type '{v}' is not narrow_peak/broad_peak; detecting automatically")
        return "auto"
    return PEAK_TYPES[v]


def die(msg):
    sys.exit(f"[setup] ERROR: {msg}")


def warn(msg):
    print(f"[setup] WARNING: {msg}", file=sys.stderr)


def first_existing(patterns):
    """First file matching any glob pattern (callers escape literal parts with glob.escape)."""
    for p in patterns:
        hits = sorted(h for h in glob.glob(p) if os.path.isfile(h))
        if hits:
            return os.path.abspath(hits[0])
    return ""


def first_dir(paths):
    for p in paths:
        if os.path.isdir(p):
            return os.path.abspath(p)
    return ""


def detect_layouts(res):
    """All (aligner, peak_type) pairs present under a pipeline --outdir."""
    out = []
    for d in sorted(glob.glob(f"{glob.escape(res)}/*/merged_library/macs3/*_peak")):
        parts = d.rstrip("/").split("/")
        if parts[-1] in ("narrow_peak", "broad_peak"):
            out.append((parts[-4], parts[-1]))
    return out


def find_lib_file(d, cid, sample, rx_suffix, what, quiet=False):
    """File in directory d named <cid><suffix> (rx_suffix: regex for the part after the id). If cid lacks _REP<n>,
    accept exactly one <cid>_REP<n><suffix>. Returns "" when nothing (or more than one candidate) matches."""
    if not cid or not os.path.isdir(d):
        return ""
    names = sorted(os.listdir(d))
    for rx, note in ((re.escape(cid) + rx_suffix, False), (re.escape(cid) + r"_REP\d+" + rx_suffix, True)):
        if note and re.search(r"_REP\d+$", cid):
            break
        hits = [n for n in names if re.fullmatch(rx, n) and os.path.isfile(os.path.join(d, n))]
        if len(hits) == 1:
            if note and not quiet:
                warn(f"{sample}: {what} '{cid}' matched {hits[0]}")
            return os.path.abspath(os.path.join(d, hits[0]))
        if len(hits) > 1:
            warn(f"{sample}: {what} '{cid}' matches several files in {d} ({', '.join(hits[:4])}); "
                 f"give the id with _REP<n> or set the file in the manifest")
            return ""
    return ""


def find_bam(ml, cid, sample):
    """<cid>.mLb.clN.sorted.bam; if cid lacks _REP<n>, accept exactly one <cid>_REP<n> BAM."""
    return find_lib_file(ml, cid, sample, r"\.mLb\.clN\.sorted\.bam", "control")


# Pipeline per-library read counts on chrM and the viral genomes (counted before host filtering), and the
# duplicate-marked alignments to the viral genomes:
#   <merged_library>/contig_reads/<id>.contig_reads.tsv
#   <merged_library>/contig_reads/viral_bam/<id>[.<anything>].bam
def find_contig_files(ml, lib_id, sample, what):
    cr = os.path.join(ml, "contig_reads")
    # the pipeline may name these files by sample (S1) rather than by replicate id (S1_REP1): try both
    ids = [lib_id] + ([re.sub(r"_REP\d+$", "", lib_id)] if re.search(r"_REP\d+$", lib_id) else [])
    tab = vb = ""
    for i in ids:
        tab = tab or find_lib_file(cr, i, sample, r"\.contig_reads\.tsv", what, quiet=True)
        vb = vb or find_lib_file(os.path.join(cr, "viral_bam"), i, sample, r"(\.[^/]*)?\.bam", what, quiet=True)
    return tab, vb


def find_peaks(res, aligner, ptype, s):
    """Exact <sample>_peaks.<ext>; if the manifest id lacks _REP<n>, accept exactly one <sample>_REP<n>."""
    ext = "narrowPeak" if ptype == "narrow_peak" else "broadPeak"
    mac = f"{res}/{aligner}/merged_library/macs3/{ptype}"
    exact = f"{mac}/{s}_peaks.{ext}"
    if os.path.isfile(exact):
        return os.path.abspath(exact), s
    if not re.search(r"_REP\d+$", s):
        hits = sorted(glob.glob(f"{glob.escape(mac)}/{glob.escape(s)}_REP[0-9]*_peaks.{ext}"))
        hits = [h for h in hits if re.fullmatch(re.escape(s) + r"_REP\d+_peaks\." + ext, os.path.basename(h))]
        if len(hits) == 1:
            return os.path.abspath(hits[0]), os.path.basename(hits[0])[: -len(f"_peaks.{ext}")]
        if len(hits) > 1:
            warn(f"{s}: several replicates match ({', '.join(os.path.basename(h) for h in hits)}); "
                 f"use the full id with _REP<n> in the manifest")
    return "", s


def main():
    if len(sys.argv) != 2:
        die("usage: 00_setup.py config.yaml")
    cfg_path = sys.argv[1]
    with open(cfg_path) as fh:
        cfg = yaml.safe_load(fh)

    outdir = os.path.abspath(cfg.get("outdir", "crossrun_out"))
    work = os.path.join(outdir, "work")
    for d in ("work", "tables", "figures", "logs"):
        os.makedirs(os.path.join(outdir, d), exist_ok=True)

    for key in ("manifest", "genome_fasta", "chrom_sizes"):
        if not cfg.get(key) or not os.path.exists(cfg[key]):
            die(f"config '{key}' missing or not found: {cfg.get(key)}")
    if not os.path.exists(cfg["genome_fasta"] + ".fai"):
        die("genome FASTA must be indexed (samtools faidx)")

    manifest_dir = os.path.dirname(os.path.abspath(cfg["manifest"]))
    rows = read_manifest(cfg["manifest"])
    if not rows:
        die("manifest is empty")
    cols = list(rows[0].keys())
    if len(cols) == 1:
        die("manifest has a single column; it must be TAB-separated (was it saved as CSV or with spaces?)")
    for c in REQUIRED_COLS:
        if c not in cols:
            die(f"manifest lacks required column '{c}'")
    ids = [r["sample_id"] for r in rows]
    dups = {i for i in ids if ids.count(i) > 1}
    if dups:
        die(f"duplicate sample_id(s): {', '.join(sorted(dups))}")

    # One genome build per analysis: peaks from different builds cannot be merged or compared
    builds = {r.get("genome_build", "").strip() for r in rows if r.get("genome_build", "").strip()}
    if len(builds) > 1:
        die(f"manifest mixes genome builds {sorted(builds)}; split into separate analyses")
    if builds and cfg.get("genome_build") and cfg["genome_build"] not in builds:
        warn(f"config genome_build '{cfg['genome_build']}' differs from manifest {builds}")

    # Covariate columns must exist
    needed = [cfg.get("batch_var"), cfg.get("series_var"), cfg.get("time_var")]
    needed += cfg.get("categorical_covariates", []) + cfg.get("numeric_covariates", [])
    needed += cfg.get("condition_vars", []) + [cfg.get("upset_group_var")]
    for c in {n for n in needed if n}:
        if c not in cols:
            die(f"covariate/grouping column '{c}' named in config is not in the manifest")

    missing = []
    resolved = []
    diag = []
    for r in rows:
        s = r["sample_id"]
        res = r["results_dir"].rstrip("/")
        if res and not os.path.isabs(res):
            res = os.path.normpath(os.path.join(manifest_dir, res))  # relative to the manifest, not the cwd
        if not os.path.isdir(res):
            die(f"{s}: results_dir does not exist: {res}")

        want_al = (norm_aligner(r.get("aligner"), s, "manifest") if r.get("aligner")
                   else norm_aligner(str(cfg.get("default_aligner", "auto")), s, "config default"))
        want_pt = (norm_peak_type(r.get("peak_type"), s, "manifest") if r.get("peak_type")
                   else norm_peak_type(str(cfg.get("default_peak_type", "auto")), s, "config default"))
        layouts = detect_layouts(res)
        if not layouts:
            die(f"{s}: no <aligner>/merged_library/macs3/<narrow_peak|broad_peak> directory under {res}; "
                f"is results_dir the pipeline --outdir?")
        cand = [(a, p) for a, p in layouts
                if want_al in ("auto", a) and want_pt in ("auto", p)]
        # prefer the layout that actually contains this sample's peaks
        with_peaks = [(a, p) for a, p in cand if find_peaks(res, a, p, s)[0]]
        pick = (with_peaks or cand or [None])[0]
        if pick is None:
            found = ", ".join(f"{a}/{p}" for a, p in layouts)
            missing.append((s, "layout", "required", f"wanted {want_al}/{want_pt}; found {found}"))
            resolved.append({**r, "results_dir": res, "peaks_file": "", "xstreme_dir": ""})
            continue
        aligner, ptype = pick
        if len(with_peaks) > 1:
            warn(f"{s}: peaks exist under several layouts {with_peaks}; using {aligner}/{ptype} "
                 f"(set aligner/peak_type columns in the manifest to choose)")
        ext = "narrowPeak" if ptype == "narrow_peak" else "broadPeak"
        ml = f"{res}/{aligner}/merged_library"
        mac = f"{ml}/macs3/{ptype}"
        peaks, pid = find_peaks(res, aligner, ptype, s)
        if r.get("peaks_file"):
            peaks, pid = r["peaks_file"], s
        # pid is the pipeline's meta.id (e.g. sample_id + _REP1 if the manifest omitted it)
        if pid != s:
            warn(f"{s}: matched pipeline id '{pid}' (manifest sample_id should include _REP<n>)")
        e = glob.escape
        paths = {
            "pipeline_id": pid,
            "peaks_file": peaks,
            "homer_file": r.get("homer_file") or first_existing([f"{e(mac)}/{e(pid)}_peaks.annotatePeaks.txt"]),
            "frip_file": r.get("frip_file") or first_existing([f"{e(mac)}/qc/{e(pid)}_peaks.FRiP_mqc.tsv",
                                                               f"{e(mac)}/qc/{e(pid)}.FRiP_mqc.tsv"]),
            "xstreme_dir": r.get("xstreme_dir") or first_dir([f"{mac}/xstreme/{pid}_xstreme"]),
            "bam_file": r.get("bam_file") or first_existing([f"{e(ml)}/{e(pid)}.mLb.clN.sorted.bam"]),
            "flagstat_file": first_existing([f"{e(ml)}/samtools_stats/{e(pid)}.mLb.clN.sorted.bam.flagstat"]),
            "feature_annot_dir": first_dir([f"{mac}/feature_annotation/{pid}"]),
        }
        ctrl = (r.get("control_id") or "").strip()
        paths["control_bam"] = find_bam(ml, ctrl, s) if ctrl else ""
        paths["contig_table"], paths["viral_bam"] = find_contig_files(ml, pid, s, "sample")
        if r.get("contig_table"):
            paths["contig_table"] = r["contig_table"]
        if r.get("viral_bam"):
            paths["viral_bam"] = r["viral_bam"]
        paths["control_contig_table"], paths["control_viral_bam"] = find_contig_files(ml, ctrl, s, "control") if ctrl else ("", "")
        paths["peak_format"] = ext
        paths["aligner_resolved"] = aligner
        paths["peak_type_resolved"] = ptype

        if not paths["peaks_file"]:
            have = sorted(os.path.basename(x) for x in glob.glob(f"{glob.escape(mac)}/*_peaks.{ext}"))
            near = [h for h in have if h.lower().startswith(s.lower()[:max(4, len(s) // 2)])]
            missing.append((s, "peaks_file", "required",
                            f"looked for {mac}/{s}_peaks.{ext}; similar: {', '.join(near[:5]) or 'none'}; "
                            f"{len(have)} peak files in that dir"))
        if not paths["xstreme_dir"]:
            xd = f"{mac}/xstreme"
            if not os.path.isdir(xd):
                note = f"{xd} does not exist (run was made without --run_xstreme?)"
            else:
                have = sorted(os.path.basename(x) for x in glob.glob(f"{glob.escape(xd)}/*_xstreme"))
                note = f"looked for {xd}/{pid}_xstreme; present: {', '.join(have[:5]) or 'none'}"
            missing.append((s, "xstreme_dir", "required", note))
        for k in ("homer_file", "frip_file", "bam_file", "flagstat_file"):
            if not paths[k]:
                missing.append((s, k, "optional", ""))
        if ctrl and not paths["control_bam"]:
            missing.append((s, "control_bam", "optional", f"control_id {ctrl}"))
        cdir = f"{ml}/contig_reads"
        if not paths["contig_table"]:
            missing.append((s, "contig_table", "optional",
                            f"looked for {cdir}/{pid}.contig_reads.tsv (chrM and viral reads are then counted from the "
                            f"filtered BAM, which has none if the pipeline removed them)"))
        if ctrl and not paths["control_contig_table"]:
            missing.append((s, "control_contig_table", "optional", f"looked for {cdir}/{ctrl}.contig_reads.tsv"))
        if paths["contig_table"] and not paths["viral_bam"]:
            missing.append((s, "viral_bam", "optional", f"looked for {cdir}/viral_bam/{pid}*.bam (viral coverage plots "
                                                        f"and the viral MAPQ table are skipped)"))
        diag.append((s, res, ",".join(f"{a}/{p}" for a, p in layouts), f"{aligner}/{ptype}"))
        resolved.append({**r, "results_dir": res, **paths})

    with open(os.path.join(work, "layouts.tsv"), "w") as fh:
        fh.write("sample_id\tresults_dir\tlayouts_found\tlayout_used\n")
        for d in diag:
            fh.write("\t".join(d) + "\n")

    extra = []
    for r in resolved:
        extra += [c for c in r.keys() if c not in cols and c not in extra]
    out_cols = cols + extra
    with open(os.path.join(work, "samples.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", extrasaction="ignore", restval="", lineterminator="\n")
        w.writeheader()
        for r in resolved:
            w.writerow(r)
    with open(os.path.join(work, "missing_files.tsv"), "w") as fh:
        fh.write("sample_id\tfile\tlevel\tdetail\n")
        for m in missing:
            fh.write("\t".join(m) + "\n")

    req_missing = [m for m in missing if m[2] == "required"]
    if req_missing:
        for m in req_missing:
            warn(f"{m[0]}: {m[1]} not found -- {m[3]}")
        warn("samples lacking peaks/XSTREME output are kept in samples.tsv but skipped by steps that need them")

    mot = cfg.get("motifs", {}) or {}
    qc = cfg.get("qc", {}) or {}
    rep = cfg.get("repeats", {}) or {}
    reg = cfg.get("regions", {}) or {}
    bns = cfg.get("bins", {}) or {}
    bmod = cfg.get("binding_models", {}) or {}
    dbs = cfg.get("known_motif_dbs") or []
    for d in dbs:
        if not os.path.exists(d):
            die(f"known motif DB not found: {d}")

    ca = cfg.get("chrom_alias") or []
    if isinstance(ca, str):
        ca = [ca]
    for f in ca:
        if not os.path.exists(f):
            die(f"chrom_alias file not found: {f}")

    ctg = cfg.get("contigs") or {}
    viral = ctg.get("viral") or {}
    if not isinstance(viral, dict):
        die("contigs.viral must map display names to contig names, e.g. {HSV-1: NC_001806.2}")
    for k, v in viral.items():
        if not re.fullmatch(r"[A-Za-z0-9._-]+", str(k)) or not re.fullmatch(r"[A-Za-z0-9._|-]+", str(v)):
            die(f"contigs.viral: '{k}: {v}' - names may contain only letters, digits, '.', '_' and '-'")
    if len(set(map(str, viral.values()))) != len(viral):
        die("contigs.viral: the same contig is listed twice")
    vfa = ctg.get("viral_fasta")
    if vfa and not os.path.exists(vfa):
        die(f"contigs.viral_fasta not found: {vfa}")
    hx = ctg.get("host_exclude", ["chrM"])
    hx = [hx] if isinstance(hx, str) else list(hx or [])
    try:
        vmq = int(ctg.get("viral_mapq", 10))
    except (TypeError, ValueError):
        die("contigs.viral_mapq must be a whole number")
    xp = cfg.get("extra_plots") or {}
    norms = xp.get("coverage_norms", ["raw", "cpm", "genome", "input"]) or []
    bad_n = [n for n in norms if n not in ("raw", "cpm", "genome", "input")]
    if bad_n:
        die(f"extra_plots.coverage_norms: unknown {bad_n} (use raw, cpm, genome, input)")
    if xp.get("coverage_bands") and not os.path.exists(xp["coverage_bands"]):
        die(f"extra_plots.coverage_bands not found: {xp['coverage_bands']}")

    # bin size for the binding models: bins.size, else the chosen size in the tuning record (tune/run_tuning.sh)
    bin_record = os.path.abspath(bns["tuning_record"]) if bns.get("tuning_record") else ""
    bin_size = bns.get("size")
    # bins.size set explicitly = bin tuning is decided: step 18 submits Level 1 and waits, so steps 19-22 follow.
    # A size taken from bins.tuning_record does not count.
    bins_size_fixed = bin_size not in (None, "")
    if bin_record and not os.path.exists(bin_record):
        die(f"bins.tuning_record not found: {bin_record}")
    if bin_size in (None, "") and bin_record:
        rec = dict(l.rstrip("\n").split("\t", 1) for l in open(bin_record) if "\t" in l)
        bin_size = rec.get("chosen_size", "")
    if bin_size not in (None, ""):
        try:
            bin_size = int(bin_size)
        except (TypeError, ValueError):
            die("bins.size must be a whole number of base pairs")
        if bin_size < 10:
            die("bins.size must be at least 10 bp")
    else:
        bin_size = ""
    if bns.get("mappability_bedgraph") and not os.path.exists(bns["mappability_bedgraph"]):
        die(f"bins.mappability_bedgraph not found: {bns['mappability_bedgraph']}")
    env = {
        "VIRAL_CONTIGS": " ".join(f"{k}={v}" for k, v in viral.items()),
        "VIRAL_FASTA": os.path.abspath(vfa) if vfa else "",
        "HOST_EXCLUDE": " ".join(map(str, hx)),
        "VIRAL_MAPQ": vmq,
        "EXTRA_ENABLED": "1" if xp.get("enabled", True) else "0",
        "EXTRA_COVERAGE": "1" if xp.get("coverage", True) else "0",
        "EXTRA_NORMS": ",".join(norms),
        "EXTRA_BANDS": os.path.abspath(xp["coverage_bands"]) if xp.get("coverage_bands") else "",
        "EXTRA_CAP": "" if xp.get("coverage_cap", 0.999) is None else xp.get("coverage_cap", 0.999),
        "EXTRA_PIES": "1" if xp.get("pies", True) else "0",
        "EXTRA_INSERT_SIZES": "1" if xp.get("insert_sizes", True) else "0",
        "CHROM_ALIAS_FILES": " ".join(os.path.abspath(f) for f in ca),
        "OUTDIR": outdir,
        "WORK": work,
        "SAMPLES": os.path.join(work, "samples.tsv"),
        "THREADS": cfg.get("threads", 4),
        "GENOME_FASTA": os.path.abspath(cfg["genome_fasta"]),
        "CHROM_SIZES": os.path.abspath(cfg["chrom_sizes"]),
        "INCLUDE_BED": os.path.abspath(cfg["include_bed"]) if cfg.get("include_bed") else "",
        "EXCLUDE_BED": os.path.abspath(cfg["exclude_bed"]) if cfg.get("exclude_bed") else "",
        "REPEATMASKER": os.path.abspath(cfg["repeatmasker"]) if cfg.get("repeatmasker") else "",
        "CENSAT": os.path.abspath(cfg["censat"]) if cfg.get("censat") else "",
        "GENES_BED": os.path.abspath(cfg["genes_bed"]) if cfg.get("genes_bed") else "",
        "KNOWN_DBS": " ".join(os.path.abspath(d) for d in dbs),
        "REPEAT_CONSENSUS_FASTA": os.path.abspath(cfg["repeat_consensus_fasta"]) if cfg.get("repeat_consensus_fasta") else "",
        "SEQ_MODE": mot.get("seq_mode", "peak"),
        "SUMMIT_FLANK": mot.get("summit_flank", 100),
        "FIMO_THRESH": mot.get("fimo_thresh", 1e-4),
        "FAMILY_EVALUE": mot.get("family_evalue", 0.05),
        "FAMILY_COOCC_MIN": "" if mot.get("family_cooccupancy_min") is None else mot.get("family_cooccupancy_min"),
        "FAMILY_COOCC_OVERLAP": mot.get("family_cooccupancy_overlap", 0.5),
        "FAMILY_COOCC_CONTAINMENT": 0 if mot.get("family_cooccupancy_containment", 0.85) is None
                                    else mot.get("family_cooccupancy_containment", 0.85),
        "FAMILY_COOCC_FLOOR": mot.get("family_cooccupancy_floor", 0.5),
        "KNOWN_EVALUE": mot.get("known_evalue", 0.05),
        "SEA_CONTROL": mot.get("sea_control", "shuffled"),
        "SPACING_MAX": mot.get("spacing_max", 400),
        "NULL_MODEL": rep.get("null_model", "uniform"),
        "GC_CANDIDATES": rep.get("gc_candidates", 10),
        "N_SHUFFLES": rep.get("n_shuffles", 100),
        "SAME_CHROM": "1" if rep.get("same_chrom", True) else "0",
        "OVERLAP_RULE": rep.get("overlap_rule", "midpoint"),
        "OVERLAP_FRACTION": rep.get("overlap_fraction", 0.5),
        "READS_IN_REPEATS": "1" if rep.get("reads_in_repeats", False) else "0",
        "MEME_PREFIX": cfg.get("meme_prefix") or "",
        "DISCOVERY": mot.get("discovery", "control_bam"),
        "DISCOVERY_ENGINE": mot.get("discovery_engine", "text"),
        "MEME_NMOTIFS": mot.get("meme_nmotifs", 5),
        "MEME_SEARCHSIZE": mot.get("meme_searchsize", 100000),
        "STREME_EXTRA_ARGS": mot.get("streme_extra_args") or "",
        "MEME_EXTRA_ARGS": mot.get("meme_extra_args") or "",
        "CONTROL_RATIO": mot.get("control_ratio", 1),
        "XSTREME_MAXW": mot.get("xstreme_maxw", 15),
        "MOTIF_MIN_WIDTH": mot.get("min_width", 6),
        "DENOVO_MAX_EVALUE": "" if mot.get("denovo_max_evalue", 0.05) is None else mot.get("denovo_max_evalue", 0.05),
        "XSTREME_EXTRA_ARGS": mot.get("xstreme_extra_args") or "",
        "QC_ENABLED": "1" if qc.get("enabled", True) else "0",
        "QC_K": qc.get("mappability_k", 150),
        "QC_MAPPABILITY_BEDGRAPH": os.path.abspath(qc["mappability_bedgraph"]) if qc.get("mappability_bedgraph") else "",
        "QC_GENMAP_INDEX": os.path.abspath(qc["genmap_index"]) if qc.get("genmap_index") else "",
        "QC_WRITE_TRACK": "1" if qc.get("write_mappability_track", True) else "0",
        "QC_INSERT_MAX": qc.get("insert_max", 1000),
        "QC_INSERT_GENOME_FRACTION": qc.get("insert_genome_fraction", 0.01),
        "BIN_SIZE": reg.get("bin_size", 100000),
        "INPUT_PEAKS_DIR": os.path.abspath(reg["input_peaks_dir"]) if reg.get("input_peaks_dir") else "",
        "MODELS_ENABLED": "1" if bmod.get("enabled", True) else "0",
        "BINS_SIZE": bin_size,
        "BINS_SIZE_FIXED": "1" if bins_size_fixed else "0",
        "BINS_RECORD": bin_record,
        "BINS_UNIQUE_MAPQ": bns.get("unique_mapq", 255),
        "BINS_INSERT_MAX": bns.get("insert_max", 2000),
        "BINS_MAPPABILITY": os.path.abspath(bns["mappability_bedgraph"]) if bns.get("mappability_bedgraph") else "",
        "BINS_DRIFT_CHECK": "1" if bns.get("drift_check", True) else "0",
        "CONFIG_YAML": os.path.abspath(cfg_path),
    }
    if cfg.get("models"):
        warn("config block 'models:' is no longer used (the regression models were replaced by the binding models: "
             "'bins:' and 'binding_models:'); it is ignored")
    try:
        mw = int(env["MOTIF_MIN_WIDTH"])
    except (TypeError, ValueError):
        die("motifs.min_width must be a whole number of base pairs")
    # STREME's own limits, from `streme --help` (MEME Suite 5.5.9): --minw must be >= 3, --maxw must be <= 30.
    # MEME has no width cap, so step 02 drops STREME above it rather than refusing the width.
    STREME_MAX_WIDTH = 30
    if mw < 3:
        die("motifs.min_width must be at least 3 bp (STREME's own minimum; shorter motifs also cannot reach "
            "motifs.fimo_thresh, so they would never be scannable)")
    try:
        xmw = int(env["XSTREME_MAXW"])
    except (TypeError, ValueError):
        die("motifs.xstreme_maxw must be a whole number of base pairs")
    if mw > xmw:
        warn(f"motifs.min_width ({mw}) is above motifs.xstreme_maxw ({xmw}): step 02 will raise the maximum to "
             f"{mw} and discovery will search that one width only. Raise xstreme_maxw to give it a range.")
    elif mw == xmw:
        warn(f"motifs.min_width and motifs.xstreme_maxw are both {mw}: discovery searches a single width, so "
             f"every motif found will be exactly {mw} bp long. Raise xstreme_maxw to give it a range.")
    # motifs.meme_searchsize -> MEME's -searchsize. From the MEME documentation: MEME samples the primary
    # sequences at random when they contain more letters than this, and 0 means use every sequence (no
    # sampling). MEME's own default is 100000, which is what this falls back to when the key is absent.
    try:
        ssz = int(env["MEME_SEARCHSIZE"])
    except (TypeError, ValueError):
        die("motifs.meme_searchsize must be a whole number of letters (0 = use all sequences, no sampling)")
    if ssz < 0:
        die("motifs.meme_searchsize must be 0 or more (0 = use all sequences, no sampling)")
    if 0 < ssz < 10000:
        warn(f"motifs.meme_searchsize ({ssz}) is far below MEME's default of 100000: discovery will search a "
             f"small random sample of each peak set, so results will depend heavily on which sample is drawn.")
    if env["DISCOVERY_ENGINE"] not in ("text", "xstreme"):
        die("motifs.discovery_engine must be text or xstreme")
    # There is no upper limit on motifs.min_width, but STREME cannot search past STREME_MAX_WIDTH, so say here
    # what step 02 will do about it rather than letting the tool refuse its own arguments per sample.
    if mw > STREME_MAX_WIDTH:
        if env["DISCOVERY_ENGINE"] == "xstreme":
            die(f"motifs.min_width ({mw}) is above STREME's maximum motif width ({STREME_MAX_WIDTH} bp), and "
                f"motifs.discovery_engine 'xstreme' runs STREME internally, so discovery cannot run. Use "
                f"motifs.discovery_engine: text, which discovers with MEME alone above {STREME_MAX_WIDTH} bp.")
        warn(f"motifs.min_width ({mw}) is above STREME's maximum motif width ({STREME_MAX_WIDTH} bp): step 02 "
             f"will skip STREME and discover with MEME alone. MEME has no width cap, but motifs.meme_nmotifs "
             f"({env['MEME_NMOTIFS']}) then caps how many motifs each sample yields, and MEME is slower at "
             f"large widths.")
    elif xmw > STREME_MAX_WIDTH:
        warn(f"motifs.xstreme_maxw ({xmw}) is above STREME's maximum motif width ({STREME_MAX_WIDTH} bp): "
             f"STREME will search {mw}-{STREME_MAX_WIDTH} bp and MEME the full {mw}-{xmw} bp.")
    if env["DISCOVERY"] not in ("control_bam", "pipeline"):
        die("motifs.discovery must be control_bam or pipeline")
    # motifs.denovo_max_evalue: one significance setting for the de novo motifs of both programs (SEA E-value)
    if env["DENOVO_MAX_EVALUE"] == "":
        warn("motifs.denovo_max_evalue is null: every de novo motif STREME and MEME report enters the motif "
             "families, significant or not.")
    else:
        try:
            dev = float(env["DENOVO_MAX_EVALUE"])
        except (TypeError, ValueError):
            die("motifs.denovo_max_evalue must be a number (SEA E-value, e.g. 0.05) or null for no gate")
        if not dev > 0:
            die("motifs.denovo_max_evalue must be above 0 (or null for no gate)")
    # motifs.family_cooccupancy_min: second family stage (site co-occupancy); null = Tomtom only
    for key, var, allow_null in (("family_cooccupancy_min", "FAMILY_COOCC_MIN", True),
                                 ("family_cooccupancy_overlap", "FAMILY_COOCC_OVERLAP", False)):
        if env[var] == "" and allow_null:
            continue
        try:
            v = float(env[var])
        except (TypeError, ValueError):
            die(f"motifs.{key} must be a number above 0 and at most 1" + (" (or null to skip the stage)" if allow_null else ""))
        if not 0 < v <= 1:
            die(f"motifs.{key} must be above 0 and at most 1" + (" (or null to skip the stage)" if allow_null else ""))
    for key, var in (("family_cooccupancy_containment", "FAMILY_COOCC_CONTAINMENT"),
                     ("family_cooccupancy_floor", "FAMILY_COOCC_FLOOR")):
        try:
            v = float(env[var])
        except (TypeError, ValueError):
            die(f"motifs.{key} must be a number between 0 and 1" + (" (or null for off)" if "containment" in key else ""))
        if not 0 <= v <= 1:
            die(f"motifs.{key} must be between 0 and 1")
    if env["QC_MAPPABILITY_BEDGRAPH"] and not os.path.exists(env["QC_MAPPABILITY_BEDGRAPH"]):
        die(f"qc.mappability_bedgraph not found: {env['QC_MAPPABILITY_BEDGRAPH']}")
    if env["NULL_MODEL"] not in ("uniform", "gc_matched"):
        die("repeats.null_model must be uniform or gc_matched")
    if env["OVERLAP_RULE"] not in ("midpoint", "fraction"):
        die("repeats.overlap_rule must be midpoint or fraction")
    with open(os.path.join(work, "config.env"), "w") as fh:
        for k, v in env.items():
            fh.write(f"export {k}={shlex.quote(str(v))}\n")

    n_ok = sum(1 for r in resolved if r["peaks_file"])
    if env["BINS_SIZE_FIXED"] == "1" and env["MODELS_ENABLED"] == "1":
        print(f"[setup] bins.size = {env['BINS_SIZE']} bp is set: bin tuning is treated as decided; step 18 submits "
              f"the Level 1 array job and waits for it (level1_mode: slurm), then the run continues to steps 19-22",
              file=sys.stderr)
    print(f"[setup] {len(resolved)} samples in manifest, {n_ok} with peaks; outputs in {outdir}")


if __name__ == "__main__":
    main()
