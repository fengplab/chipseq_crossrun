#!/usr/bin/env python3
"""
Assemble every figure into one self-contained HTML report (images embedded; no tables), with a short
description before each figure (meaning, axes, abbreviations) and, at the end, a data dictionary for
every column of every table the bundle wrote. Print to PDF from a browser for a paginated document.

Usage: build_report.py config.yaml        ->  <outdir>/report/crossrun_report.html
"""
import base64
import csv
import datetime
import glob
import html
import os
import re
import sys

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
esc = html.escape


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def header(path):
    with open(path) as fh:
        return fh.readline().rstrip("\n\r").split("\t")


def subcaption(prefix, stem, cfg):
    """Readable sub-figure caption from the file name."""
    rest = stem[len(prefix):].strip("_")
    if prefix in ("fig04", "fig07"):
        return "Centromeric/satellite (CenSat) classes" if "censat" in rest else "RepeatMasker classes"
    if prefix == "fig12":
        return "Samples grouped by '%s'" % rest.replace("upset_", "", 1)
    if prefix == "fig15":
        return "Series: %s" % rest.replace("timecourse_", "", 1)
    page = ""
    m = re.search(r"_p(\d+)$", rest)
    if m and prefix in ("fig18", "fig19"):
        page = " (part %s)" % m.group(1)
        rest = rest[:m.start()]
    if prefix == "fig18":
        return ("Motif families" + page) if page else ""
    if prefix == "fig19":
        return {"reads_class": "IP vs control reads per repeat class", "reads_peaks": "IP vs control reads in peaks",
                "peak_share": "Share of peaks per repeat class", "tandem": "Motif tandem array score",
                "density": "Motif density per kb in peaks",
                "density_class": "Motif density per kb inside each repeat class"}.get(rest.replace("covariate_effects_", "", 1), rest) + page
    if prefix == "fig20":
        return "Series: %s" % rest.replace("timepoint_contrasts_", "", 1)
    if prefix == "qc03":
        return "Sample: %s" % rest.replace("insert_", "", 1)
    return ""


CSS = """
body { font-family: -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif; color: #1f1f1d; max-width: 1150px;
       margin: 0 auto; padding: 24px 32px; line-height: 1.5; font-size: 14px; }
h1 { font-size: 26px; margin-bottom: 4px; } h2 { font-size: 20px; border-bottom: 2px solid #2a78d6; padding-bottom: 4px; margin-top: 40px; }
h3 { font-size: 16px; margin: 26px 0 6px; } h4 { font-size: 14px; margin: 18px 0 4px; color: #444; }
.meta { color: #555; } .note { background: #f3f7fc; border-left: 3px solid #2a78d6; padding: 6px 10px; margin: 8px 0; }
.empty { background: #fbf6ec; border-left: 3px solid #eda100; padding: 6px 10px; margin: 8px 0; }
figure { margin: 10px 0 20px; } figure img { max-width: 100%; height: auto; border: 1px solid #e3e2dc; }
figcaption { font-size: 13px; color: #555; margin-top: 4px; }
ul.axes { margin: 4px 0 8px 18px; padding: 0; } ul.axes li { margin: 2px 0; }
dl.abbr { display: grid; grid-template-columns: max-content 1fr; gap: 2px 14px; margin: 4px 0 8px; }
dl.abbr dt { font-weight: 600; } dl.abbr dd { margin: 0; }
table.dict { border-collapse: collapse; width: 100%; font-size: 12.5px; margin: 6px 0 18px; }
table.dict th, table.dict td { border: 1px solid #dddcd5; padding: 4px 8px; text-align: left; vertical-align: top; }
table.dict th { background: #f0efec; } table.dict td.col { font-family: Menlo, Consolas, monospace; white-space: nowrap; }
.toc a { text-decoration: none; color: #2a78d6; } .toc li { margin: 2px 0; }
.block { break-inside: avoid; }
@media print { body { max-width: none; padding: 0 10mm; font-size: 11px; }
  .figblock { break-before: page; } h2 { break-before: page; } figure img { max-height: 235mm; object-fit: contain; }
  .dictblock { break-inside: avoid; } }
"""


def main():
    cfg_path = sys.argv[1]
    cfg = yaml.safe_load(open(cfg_path))
    out = os.path.abspath(cfg.get("outdir", "crossrun_out"))
    figdir, tabdir, work = os.path.join(out, "figures"), os.path.join(out, "tables"), os.path.join(out, "work")
    spec = yaml.safe_load(open(os.path.join(ROOT, "report", "figures.yaml")))
    dictionary = read_tsv(os.path.join(ROOT, "report", "data_dictionary.tsv"))
    notes = {}
    nf = os.path.join(figdir, "figure_notes.tsv")
    if os.path.exists(nf):
        for r in read_tsv(nf):
            notes[r["figure"]] = r["note"]
    samples = read_tsv(os.path.join(work, "samples.tsv"))
    sample_ids = {s["sample_id"] for s in samples}
    path_cols = {"results_dir", "pipeline_id", "peaks_file", "homer_file", "frip_file", "xstreme_dir", "bam_file",
                 "flagstat_file", "feature_annot_dir", "control_bam", "peak_format", "aligner_resolved",
                 "peak_type_resolved", "aligner", "peak_type"}
    meta_cols = [c for c in (samples[0].keys() if samples else []) if c not in path_cols and c != "sample_id"]

    parts = []
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    rep = cfg.get("repeats", {}) or {}
    dbs = ", ".join(os.path.basename(d) for d in (cfg.get("known_motif_dbs") or [])) or "none"
    parts.append(f"<h1>ChIP-seq cross-run summary</h1><p class='meta'>Generated {esc(now)} &middot; "
                 f"{len(sample_ids)} samples &middot; genome {esc(str(cfg.get('genome_build', '')))}</p>")
    parts.append("<div class='note'><b>Settings.</b> "
                 f"Repeat enrichment null model: {esc(str(rep.get('null_model', 'uniform')))}, "
                 f"{esc(str(rep.get('n_shuffles', '')))} shuffles, overlap rule '{esc(str(rep.get('overlap_rule', 'midpoint')))}'. "
                 f"Known-motif databases: {esc(dbs)}. Motif families in the regression models: "
                 f"{esc(str((cfg.get('models') or {}).get('family_stat', 'sea_pvalue')))} <= {esc(str((cfg.get('models') or {}).get('family_cutoff')))}; "
                 f"peaks with q <= {esc(str((cfg.get('models') or {}).get('peak_qvalue')))}. Figure numbers match the figure file names; optional figures that were not "
                 "produced in this run are skipped. Tables are not included here; see the <code>tables/</code> folder and the "
                 "data dictionary at the end.</div>")

    # ---------------- figures ----------------
    toc, body, fig_no, n_figs = [], [], 0, 0
    for sec in spec["sections"]:
        sec_parts = []
        for prefix in sec["figures"]:
            f = spec["figures"][prefix]
            files = sorted(p for p in glob.glob(os.path.join(figdir, prefix + "_*.png")))
            summary = notes.get(prefix + "__summary") or next((v for k, v in notes.items() if k.startswith(prefix + "_") and k.endswith("__summary")), None)
            if not files and not summary:
                continue
            # numbering matches the file names: fig07_* = Figure 7, qc02_* = Figure QC2, pk01_* = Figure P1
            fig_no = (str(int(prefix[3:])) if prefix.startswith("fig")
                      else "P" + str(int(prefix[2:])) if prefix.startswith("pk") else "QC" + str(int(prefix[2:])))
            n_figs += 1
            anchor = f"fig{fig_no}"
            toc.append(f"<li><a href='#{anchor}'>Figure {fig_no}. {esc(f['title'])}</a></li>")
            b = [f"<div class='figblock' id='{anchor}'><h3>Figure {fig_no}. {esc(f['title'])}</h3>",
                 f"<p>{esc(' '.join(f['meaning'].split()))}</p>",
                 "<h4>How to read it</h4><ul class='axes'>" + "".join(f"<li>{esc(a)}</li>" for a in f.get("axes", [])) + "</ul>"]
            if f.get("abbreviations"):
                b.append("<h4>Abbreviations</h4><dl class='abbr'>" +
                         "".join(f"<dt>{esc(str(k))}</dt><dd>{esc(str(v))}</dd>" for k, v in f["abbreviations"].items()) + "</dl>")
            if summary:
                b.append(f"<div class='empty'>{esc(summary)}</div>")
            for i, p in enumerate(files):
                stem = os.path.splitext(os.path.basename(p))[0]
                sub = subcaption(prefix, stem, cfg)
                label = f"Figure {fig_no}{chr(97 + i) if len(files) > 1 else ''}"
                note = notes.get(stem, "")
                with open(p, "rb") as fh:
                    data = base64.b64encode(fh.read()).decode()
                b.append("<figure class='block'>"
                         + (f"<h4>{label}. {esc(sub)}</h4>" if sub else "")
                         + (f"<div class='note'>{esc(note)}</div>" if note else "")
                         + f"<img alt='{esc(label)}' src='data:image/png;base64,{data}'>"
                         + f"<figcaption>{label}{(' - ' + esc(sub)) if sub else ''}. File: figures/{esc(stem)}.pdf</figcaption></figure>")
            b.append("</div>")
            sec_parts.append("".join(b))
        if sec_parts:
            body.append(f"<h2>{esc(sec['title'])}</h2>" + "".join(sec_parts))

    # ---------------- data dictionary ----------------
    exact, generic, regex, sample_desc, table_desc = {}, {}, [], {}, {}
    for r in dictionary:
        t, c, d = r["table"], r["column"], r["description"]
        if c == "_table_":
            table_desc[t] = d
        elif c == "@sample":
            sample_desc[t] = d
        elif c == "@meta":
            generic["@meta"] = d
        elif c.startswith("re:"):
            regex.append((t, re.compile(c[3:]), d))
        elif t == "*":
            generic[c] = d
        else:
            exact[(t, c)] = d

    def describe(t, c):
        if (t, c) in exact:
            return exact[(t, c)]
        for rt, rx, d in regex:
            if rt in (t, "*"):
                m = rx.match(c)
                if m:
                    return d.replace("{1}", m.group(1) if m.groups() else c)
        if c in sample_ids and t in sample_desc:
            return f"{sample_desc[t]} (sample {c})"
        if c in generic:
            return generic[c]
        if c in meta_cols:
            return generic["@meta"].replace("{col}", c)
        return None

    missing = []
    dict_parts = ["<h2 id='dictionary'>Data dictionary</h2><p>Every column of every table written to "
                  "<code>tables/</code>. Wide matrices with one column per sample are described once for all sample columns.</p>"]
    for path in sorted(glob.glob(os.path.join(tabdir, "*.tsv"))):
        t = os.path.splitext(os.path.basename(path))[0]
        cols = header(path)
        rows_html, sample_cols_done = [], False
        for c in cols:
            if c in sample_ids and t in sample_desc:
                if sample_cols_done:
                    continue
                sample_cols_done = True
                n = sum(1 for x in cols if x in sample_ids)
                rows_html.append(f"<tr><td class='col'>&lt;sample_id&gt; ({n} columns)</td><td>{esc(sample_desc[t])}</td></tr>")
                continue
            d = describe(t, c)
            if d is None:
                missing.append(f"{t}.{c}")
                d = "(not documented)"
            rows_html.append(f"<tr><td class='col'>{esc(c)}</td><td>{esc(d)}</td></tr>")
        dict_parts.append(f"<div class='dictblock'><h3>{esc(t)}.tsv</h3><p>{esc(table_desc.get(t, ''))}</p>"
                          "<table class='dict'><tr><th>Column</th><th>Description</th></tr>" + "".join(rows_html) + "</table></div>")
    toc.append("<li><a href='#dictionary'>Data dictionary</a></li>")

    parts.append("<h2 style='break-before: avoid'>Contents</h2><ol class='toc' style='list-style:none;padding-left:0'>" + "".join(toc) + "</ol>")
    parts += body + dict_parts
    os.makedirs(os.path.join(out, "report"), exist_ok=True)
    dest = os.path.join(out, "report", "crossrun_report.html")
    with open(dest, "w") as fh:
        fh.write("<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>"
                 "<meta name='viewport' content='width=device-width, initial-scale=1'>"
                 f"<title>ChIP-seq cross-run summary</title><style>{CSS}</style></head><body>"
                 + "".join(parts) + "</body></html>")
    if missing:
        print(f"[report] WARNING: {len(missing)} table columns lack a dictionary entry: {', '.join(missing[:10])}", file=sys.stderr)
    print(f"[report] {n_figs} figures, {len(glob.glob(os.path.join(tabdir, '*.tsv')))} tables documented -> {dest}", file=sys.stderr)


if __name__ == "__main__":
    main()
