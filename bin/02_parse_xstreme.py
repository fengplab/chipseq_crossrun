#!/usr/bin/env python3
"""
Step 02: collect XSTREME results from every sample.

Writes
  tables/xstreme_motifs.tsv     one row per motif reported in xstreme.tsv (de novo and known)
  tables/xstreme_sea.tsv        SEA statistics for every motif XSTREME scored (sea_out/sea.tsv)
  tables/xstreme_run_info.tsv   per sample: --maxw used, motif counts, whether a known DB was used
  work/motifs/all_discovered.meme  every de novo motif (combined.meme) that passes motifs.min_width and
                                   motifs.denovo_max_evalue, renamed <sample>__<motif_id>

Significance gate (motifs.denovo_max_evalue, one setting for STREME and MEME alike): a de novo motif enters the
pooled set only when SEA's E-value for it, against that sample's control sequences, is at or below the setting.
SEA scores the motifs of both programs in one run with one test, so the same number means the same thing for
both; the programs' own E-values are computed differently and are not comparable with each other. A motif SEA
did not report is treated as not significant. Counts in and out are logged per sample and per program.

Usage: 02_parse_xstreme.py <work>/samples.tsv <outdir> [xstreme_root]
  xstreme_root: read <xstreme_root>/<sample>_xstreme instead of the pipeline's XSTREME directories
                (used when motif discovery is rerun with control-library sequences, step 02)
"""
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from memelib import read_meme, write_meme  # noqa: E402


def norm(h):
    return re.sub(r"[^0-9a-z]+", "_", h.strip().lower()).strip("_")


def read_tsv_skip_comments(path):
    rows, header = [], None
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                header = [norm(x) for x in f]
                continue
            rows.append(dict(zip(header, f)))
    return header or [], rows


def main():
    samples_tsv, outdir = sys.argv[1], sys.argv[2]
    xroot = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] else None
    work = os.path.join(outdir, "work")
    os.makedirs(os.path.join(work, "motifs"), exist_ok=True)
    tables = os.path.join(outdir, "tables")

    with open(samples_tsv) as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))

    motif_cols = ["seed_motif", "rank", "source", "id", "alt_id", "consensus", "width", "sites",
                  "sea_pvalue", "evalue", "evalue_acc", "sim_source", "sim_motif"]
    all_motifs, bg0 = [], None
    # motifs.min_width: de novo motifs shorter than this are left out of the pooled motif set (and so of the families).
    # Discovery itself is also run with this minimum; the filter covers motifs from the pipeline's own XSTREME.
    min_w = int(os.environ.get("MOTIF_MIN_WIDTH", "0") or 0)
    n_short = 0
    # motifs.denovo_max_evalue: empty = no significance gate
    ev_raw = (os.environ.get("DENOVO_MAX_EVALUE", "") or "").strip()
    max_ev = float(ev_raw) if ev_raw else None
    n_notsig = n_nosea = 0
    fm = open(os.path.join(tables, "xstreme_motifs.tsv"), "w")
    fm.write("\t".join(["sample_id", "motif_uid", "is_denovo", "in_combined_meme"] + motif_cols
                       + ["sea_evalue", "passes_denovo_evalue"]) + "\n")
    fs = open(os.path.join(tables, "xstreme_sea.tsv"), "w")
    fs.write("sample_id\tmotif_uid\tdb\tid\talt_id\tconsensus\ttp\ttp_pct\tfp\tfp_pct\tenr_ratio\tpvalue\tlog_pvalue\tevalue\tqvalue\n")
    fi = open(os.path.join(tables, "xstreme_run_info.tsv"), "w")
    fi.write("sample_id\tmaxw_used\tn_denovo_motifs\tn_known_rows\tknown_db_used\tstatus\tcontrol_type\tn_control_seqs\tn_primary_seqs\tdiscovery_engine"
             "\tn_denovo_reported\tn_denovo_not_significant\tdenovo_max_evalue\n")

    for s in samples:
        sid, xd = s["sample_id"], s.get("xstreme_dir", "")
        if xroot:
            xd = os.path.join(xroot, f"{sid}_xstreme")
        if not xd or not os.path.isdir(xd):
            fi.write(f"{sid}\tNA\t0\t0\tNA\tmissing\tNA\tNA\tNA\tNA\tNA\tNA\t{ev_raw or 'none'}\n")
            continue
        maxw = "NA"
        mw = os.path.join(xd, "pipeline_maxw.txt")
        if os.path.exists(mw):
            maxw = open(mw).read().split()[-1]

        # SEA E-value of every motif SEA scored for this sample (read first: the gate below needs it)
        st = os.path.join(xd, "sea_out", "sea.tsv")
        sea_ev = {}
        if os.path.exists(st):
            with open(st) as fh:
                hdr = None
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if hdr is None:
                        hdr = {h.strip(): i for i, h in enumerate(f)}
                        continue
                    if "ID" in hdr and "EVALUE" in hdr and max(hdr["ID"], hdr["EVALUE"]) < len(f):
                        try:
                            sea_ev[f[hdr["ID"]]] = float(f[hdr["EVALUE"]])
                        except ValueError:
                            pass

        combined = os.path.join(xd, "combined.meme")
        if not os.path.exists(combined):
            combined = os.path.join(xd, "xstreme.txt")
        denovo_ids, kept_ids, sig_ids = set(), set(), set()
        s_notsig = s_nosea = 0
        if os.path.exists(combined):
            bg, mots = read_meme(combined)
            bg0 = bg0 or bg
            if max_ev is not None and mots and not sea_ev:
                sys.exit(f"[02] ERROR: {sid}: motifs.denovo_max_evalue = {ev_raw} needs SEA results, but {st} is missing "
                         f"or has no ID/EVALUE columns. Rerun step 02 for this sample, or set denovo_max_evalue: null.")
            for m in mots:
                denovo_ids.add(m.id)
                if max_ev is None or (m.id in sea_ev and sea_ev[m.id] <= max_ev):
                    sig_ids.add(m.id)
                if m.width < min_w:
                    n_short += 1
                    continue
                if m.id not in sig_ids:
                    s_notsig += 1
                    s_nosea += m.id not in sea_ev
                    continue
                kept_ids.add(m.id)
                m.alt = m.id  # keep the original id as alt name
                m.id = f"{sid}__{m.id}"
                all_motifs.append(m)

        xt = os.path.join(xd, "xstreme.tsv")
        by_src = {}
        n_known = 0
        known_used = False
        if os.path.exists(xt):
            _, rows = read_tsv_skip_comments(xt)
            for r in rows:
                src = r.get("source", "")
                is_denovo = src.upper() in ("STREME", "MEME") or r.get("id", "") in denovo_ids
                if not is_denovo:
                    n_known += 1
                if r.get("sim_motif", "") or not is_denovo:
                    known_used = True
                uid = f"{sid}__{r.get('id', '')}"
                vals = [r.get(c, "") for c in motif_cols]
                rid = r.get("id", "")
                fm.write("\t".join([sid, uid, str(is_denovo).upper(), str(rid in kept_ids).upper()] + vals
                                   + [("%.4g" % sea_ev[rid]) if rid in sea_ev else "NA",
                                      (str(rid in sig_ids).upper() if is_denovo else "NA")]) + "\n")
                if is_denovo and rid in denovo_ids:
                    by_src.setdefault(src.upper() or "?", [0, 0])
                    by_src[src.upper() or "?"][0] += 1
                    by_src[src.upper() or "?"][1] += rid in kept_ids
        status = "ok" if os.path.exists(xt) else "no_xstreme_tsv"

        if os.path.exists(st):
            with open(st) as fh:
                header = None
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if header is None:
                        header = [h.strip() for h in f]
                        idx = {h: i for i, h in enumerate(header)}
                        continue
                    g = lambda k: f[idx[k]] if k in idx and idx[k] < len(f) else ""  # noqa: E731
                    fs.write("\t".join([sid, f"{sid}__{g('ID')}", g("DB"), g("ID"), g("ALT_ID"), g("CONSENSUS"),
                                        g("TP"), g("TP%"), g("FP"), g("FP%"), g("ENR_RATIO"), g("PVALUE"),
                                        g("LOG_PVALUE"), g("EVALUE"), g("QVALUE")]) + "\n")
        ci = {"control_type": "shuffled (pipeline)", "n_control_seqs": "NA", "n_primary_seqs": "NA", "engine": "pipeline XSTREME"}
        cif = os.path.join(xd, "control_info.tsv")
        if os.path.exists(cif):
            with open(cif) as fh:
                ci = next(csv.DictReader(fh, delimiter="\t"))
        fi.write(f"{sid}\t{maxw}\t{len(kept_ids)}\t{n_known}\t{str(known_used).upper()}\t{status}\t"
                 f"{ci['control_type']}\t{ci['n_control_seqs']}\t{ci['n_primary_seqs']}\t{ci.get('engine', 'xstreme')}"
                 f"\t{len(denovo_ids)}\t{s_notsig}\t{ev_raw or 'none'}\n")
        n_notsig += s_notsig; n_nosea += s_nosea
        if max_ev is not None and denovo_ids:
            per = ", ".join(f"{k} {v[1]}/{v[0]}" for k, v in sorted(by_src.items())) or "no source column"
            print(f"[02] {sid}: {len(denovo_ids)} de novo motifs reported, {len(kept_ids)} kept "
                  f"(SEA E-value <= {ev_raw}; kept/reported by program: {per})"
                  + (f"; {s_nosea} had no SEA result" if s_nosea else ""), file=sys.stderr)
            if not kept_ids:
                print(f"[02] WARNING: {sid}: no de novo motif passes motifs.denovo_max_evalue = {ev_raw}", file=sys.stderr)

    fm.close(); fs.close(); fi.close()
    # replace the pooled motif file only if its content changed, so downstream caches stay valid
    dest = os.path.join(work, "motifs", "all_discovered.meme")
    write_meme(dest + ".tmp", all_motifs, bg0)
    if os.path.exists(dest) and open(dest).read() == open(dest + ".tmp").read():
        os.remove(dest + ".tmp")
    else:
        os.replace(dest + ".tmp", dest)
    print(f"[02] {len(all_motifs)} de novo motifs from {len(samples)} samples"
          + (f"; {n_short} shorter than motifs.min_width = {min_w} bp left out" if n_short else "")
          + (f"; {n_notsig} above motifs.denovo_max_evalue = {ev_raw} left out"
             + (f" ({n_nosea} of them with no SEA result)" if n_nosea else "") if max_ev is not None else
             "; no significance gate (motifs.denovo_max_evalue: null)"), file=sys.stderr)


if __name__ == "__main__":
    main()
