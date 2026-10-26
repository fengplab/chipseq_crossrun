#!/usr/bin/env Rscript
## Read coverage along a whole chromosome (or part of one), for one or more ChIP libraries, with optional
## annotation tracks underneath: chromosome bands (ideogram) and CenSat arrays.
##
## Usage:
##   Rscript R/plot_chrom_coverage.R --chrom chr9 <tracks> [options]
##
## Tracks (at least one source; several files are comma-separated, one panel each):
##   --bam A.bam,B.bam          indexed BAM files; coverage is computed with samtools
##   --control C.bam[,D.bam]    control (input/IgG) BAM: one per --bam, or one for all. Used by --norm input,
##                              otherwise drawn in grey behind its IP (turn off with --no-control-track)
##   --bigwig X.bw,Y.bw         pre-computed coverage (needs bigWigToBedGraph or the pyBigWig Python module)
##   --bedgraph X.bg[.gz]       pre-computed coverage as bedGraph
##   --crossrun DIR --samples s1,s2
##                              take BAM, control BAM and peaks of these samples from a crossrun output folder,
##                              and (unless given) the CenSat annotation and chromosome aliases
##   --viral-bam                with --crossrun: read the pipeline's duplicate-marked viral BAMs (samples.tsv viral_bam,
##                              control_viral_bam) instead of the filtered BAMs, for a viral genome; library sizes for
##                              cpm / genome / input still come from the filtered BAMs' host nuclear reads
##   --labels a,b               panel labels (default: sample or file names)
##   --peaks P1.bed,P2.bed      peak file per track (narrowPeak/BED), drawn as ticks under the coverage
##
## Normalisation (BAM tracks; bigWig/bedGraph values are plotted as they are). Several, comma-separated, give one
## file each from a single coverage computation (e.g. --norm raw,cpm,genome,input):
##   --norm raw      mean read depth per bin (default)
##   --norm cpm      depth per million reads on the host nuclear chromosomes (compare libraries of different depth)
##   --norm genome   fold over the library's mean depth on the host nuclear chromosomes (1 = average coverage; like
##                   deepTools RPGC)
##   --norm input    log2(IP / control) per bin, each scaled to its host nuclear library size (needs --control)
##   Host nuclear chromosomes = genome.sizes of the crossrun output (no viral genomes, no chrM), or chr1-22, X, Y.
##   --mapq N        ignore reads with mapping quality below N (default 0 = keep multi-mappers' primary hits)
##
## Region and bins:
##   --region 120e6-150e6       part of the chromosome (bp; scientific notation allowed). Default: all of it
##   --bin N                    bin size in bp (default: ~4,000 bins over the region, rounded)
##   --smooth K                 running mean over K bins
##
## Annotations (BED, BED.gz or bigBed; any chromosome naming covered by the alias table):
##   --bands FILE               cytoband file: chrom start end name gieStain (UCSC/T2T cytoBand format)
##   --censat FILE              T2T CenSat annotation (raw names like hsat3_1(A1), or class BED4)
##   --chrom-alias FILE         two columns: name used in a file <TAB> name to use (e.g. NC_060933.1  chr9)
##
## Display and output:
##   --ylog                     log10(1 + value) y-axis (not for --norm input)
##   --ymax N | --cap Q         clip the y-axis at N, or at quantile Q of all bins (e.g. 0.999); clipped bins
##                              are marked with a triangle at the top
##   --same-y                   one y-axis range for all tracks
##   --title TEXT, --out FILE (.png/.pdf/.svg; default coverage_<chrom>.png), --width W, --height H (inches)
##                              {chrom} and {norm} in --out are replaced; with several --norm values and no {norm},
##                              _<norm> is added before the extension
##   --px-per-bin N             widen the figure so every bin gets at least N pixels in the coverage panels (overrides
##                              --width; 2 recommended: a 1-pixel bin usually straddles two pixels and is drawn faint).
##                              PNG is limited to 32,767 pixels wide: use a larger --bin, a --region, or .pdf/.svg.
##                              Without --bin: the smallest bin (>= 10 bp) that fits, and the figure at the full 32,767 px
##   --dpi N                    resolution of PNG output (default 200)
##   --label-pos left|center|right   where each track's sample label sits along the top of its panel (default left)
##   --cache DIR                keep the binned BAM counts here and reuse them (same BAM, contig, bins and MAPQ)
##   --table FILE               also write the binned values (TSV)
##   --threads N                parallel samtools jobs (default 4)
##
## Top peaks (narrowPeak/broadPeak files: --crossrun reads each sample's peak file from samples.tsv, or give --peaks):
##   --label-top N              number the N highest-ranked peaks in the plotted region of each track (1 = best): the
##                              numbers sit in a band under the track's peak ticks, joined to each summit by a thin line,
##                              never on the coverage; the ranked peaks are written to --top-out
##   --top-by score|signal|p|q  ranking column: 5 score (default), 7 signal value, 8 -log10 p, 9 -log10 q
##   --top-out FILE             TSV of the ranked peaks (default <out without extension>_top_peaks.tsv); input for
##                              R/plot_peak_sequences.R
##
## Examples:
##   Rscript R/plot_chrom_coverage.R --chrom chr9 --crossrun crossrun_out_human_broad \
##       --samples THP1_HSV1_8h_IP-cGAS,THP1_IAV_6h_IP-cGAS --norm genome --bands chm13v2.0_cytobands_allchrs.bed
##   Rscript R/plot_chrom_coverage.R --chrom chr1 --bam ip.bam --control input.bam --norm input \
##       --censat censat.bb --region 120e6-150e6
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
options(datatable.integer64 = "double")

## ------------------------------------------------------------------ arguments
self <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
a <- commandArgs(trailingOnly = TRUE)
if (!length(a) || any(a %in% c("-h", "--help"))) {
  cat(sub("^## ?", "", grep("^##", readLines(self), value = TRUE)), sep = "\n"); quit(status = if (length(a)) 0 else 1)
}
flags <- c("--ylog", "--same-y", "--no-control-track", "--viral-bam")
opt <- function(name, default = NULL) { i <- match(name, a); if (is.na(i)) default else a[i + 1] }
has <- function(name) name %in% a
known <- c(flags, "--chrom", "--bam", "--control", "--bigwig", "--bedgraph", "--crossrun", "--samples", "--labels", "--peaks",
           "--norm", "--mapq", "--region", "--bin", "--smooth", "--bands", "--censat",
           "--chrom-alias", "--ymax", "--cap", "--title", "--out", "--width", "--height", "--px-per-bin", "--dpi", "--label-pos", "--label-top", "--top-by", "--top-out", "--table", "--threads", "--cache")
bad <- setdiff(grep("^--", a, value = TRUE), known)
if (length(bad)) { message("ERROR: unknown option(s): ", paste(bad, collapse = ", "), " (see --help)"); quit(status = 1) }
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
csv <- function(x) if (is.null(x)) character() else strsplit(x, ",")[[1]]
die <- function(...) { message("ERROR: ", ...); quit(status = 1) }
say <- function(...) message("[coverage] ", ...)
num <- function(x) as.numeric(x)

chrom <- opt("--chrom") %||% die("--chrom is required (e.g. --chrom chr9)")
norms <- csv(opt("--norm", "raw"))
badn <- setdiff(norms, c("raw", "cpm", "genome", "input")); if (length(badn)) die("unknown --norm: ", paste(badn, collapse = ", "))
mapq <- as.integer(opt("--mapq", 0)); threads <- as.integer(opt("--threads", 4))
smooth_k <- as.integer(opt("--smooth", 1))
cr <- opt("--crossrun")
bin_root <- file.path(dirname(dirname(self)), "bin")
CHR <- chrom
tmp <- tempfile("chromcov_"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
need <- function(tool) if (!nzchar(Sys.which(tool))) die(tool, " not found on PATH")

## ------------------------------------------------------------------ chromosome aliases
alias <- data.table(from = character(), to = character())
af <- opt("--chrom-alias") %||% (if (!is.null(cr)) file.path(cr, "work", "annot", "chrom_alias.tsv") else NULL)
if (!is.null(af) && file.exists(af)) alias <- fread(af, header = FALSE, col.names = c("from", "to"), select = 1:2, colClasses = "character")
names_for <- unique(c(chrom, alias[to == chrom, from]))          # every spelling of the requested chromosome
canon <- function(x) { i <- match(x, alias$from); ifelse(is.na(i), x, alias$to[i]) }

## ------------------------------------------------------------------ tracks
tr <- data.table(kind = character(), file = character(), control = character(), label = character(), peaks = character(),
                 libfile = character(), libcontrol = character(), rawpeaks = character())
if (!is.null(cr)) {
  sm <- fread(file.path(cr, "work", "samples.tsv"), colClasses = "character")
  ss <- csv(opt("--samples")) %||% die("--crossrun needs --samples")
  miss <- setdiff(ss, sm$sample_id); if (length(miss)) die("not in ", cr, "/work/samples.tsv: ", paste(miss, collapse = ", "))
  r <- sm[match(ss, sample_id)]
  pk <- file.path(cr, "work", "peaks", paste0(ss, ".bed5"))
  if (has("--viral-bam")) {
    if (!"viral_bam" %in% names(r) || any(is.na(r$viral_bam) | !nzchar(r$viral_bam))) die("--viral-bam: no viral_bam in samples.tsv for ", paste(ss, collapse = ", "))
    vctl <- if ("control_viral_bam" %in% names(r)) r$control_viral_bam else NA_character_
    tr <- rbind(tr, data.table(kind = "bam", file = r$viral_bam, control = vctl, label = ss, peaks = NA_character_,
                               libfile = r$bam_file, libcontrol = r$control_bam %||% NA_character_, rawpeaks = NA_character_))
  } else
  tr <- rbind(tr, data.table(kind = "bam", file = r$bam_file, control = r$control_bam %||% NA_character_, label = ss,
                             peaks = ifelse(file.exists(pk), pk, NA_character_), libfile = NA_character_, libcontrol = NA_character_,
                             rawpeaks = if ("peaks_file" %in% names(r)) r$peaks_file else NA_character_))
}
add <- function(kind, files) if (length(files)) tr <<- rbind(tr, data.table(kind = kind, file = files, control = NA_character_,
                                                                           label = sub("\\.mLb\\.clN(\\.sorted)?$", "", sub("\\.(bam|bw|bigwig|bigWig|bedGraph|bedgraph|bg)(\\.gz)?$", "", basename(files))),
                                                                           peaks = NA_character_, libfile = NA_character_, libcontrol = NA_character_,
                                                                           rawpeaks = NA_character_))
add("bam", csv(opt("--bam"))); add("bigwig", csv(opt("--bigwig"))); add("bedgraph", csv(opt("--bedgraph")))
if (!nrow(tr)) die("no tracks: give --bam, --bigwig, --bedgraph or --crossrun with --samples")
ctl <- csv(opt("--control"))
if (length(ctl)) {
  ib <- which(tr$kind == "bam" & is.na(tr$control))
  if (!(length(ctl) %in% c(1, length(ib)))) die("--control needs one file, or one per --bam (", length(ib), ")")
  tr$control[ib] <- rep_len(ctl, length(ib))
}
lab <- csv(opt("--labels")); if (length(lab)) { if (length(lab) != nrow(tr)) die("--labels needs ", nrow(tr), " names"); tr$label <- lab }
pkf <- csv(opt("--peaks")); if (length(pkf)) { if (length(pkf) != nrow(tr)) die("--peaks needs one file per track (", nrow(tr), ")"); tr$peaks <- pkf; tr$rawpeaks <- pkf }
tr[control %in% c("", "NA"), control := NA_character_]
tr[libfile %in% c("", "NA"), libfile := NA_character_]; tr[libcontrol %in% c("", "NA"), libcontrol := NA_character_]
for (f in c(tr$file, na.omit(tr$control))) if (!file.exists(f)) die("file not found: ", f)
if ("input" %in% norms && any(tr$kind == "bam" & is.na(tr$control))) {
  if (length(norms) == 1) die("--norm input needs a control BAM for every BAM track")
  say("no control BAM for every track: the 'input' plot is skipped"); norms <- setdiff(norms, "input")
}
if (any(norms != "raw") && any(tr$kind != "bam")) say("note: --norm applies to BAM tracks only; bigWig/bedGraph values are plotted as they are")
## host nuclear chromosomes: the library size and genome-wide mean exclude viral genomes and chrM
gs_f <- if (!is.null(cr)) file.path(cr, "work", "annot", "genome.sizes") else NULL
host_set <- if (!is.null(gs_f) && file.exists(gs_f)) fread(gs_f, header = FALSE)$V1 else NULL
is_host <- function(x) { cx <- canon(x); if (!is.null(host_set)) cx %in% host_set else grepl("^chr([0-9]+|X|Y)$", cx) }
if (any(tr$kind == "bam")) need("samtools")

## ------------------------------------------------------------------ chromosome length, region, bins
bam_chrom <- function(bam, libbam = NA) {
  h <- fread(cmd = sprintf("samtools idxstats '%s'", bam), header = FALSE, col.names = c("chr", "len", "mapped", "unmapped"))
  if (!is.na(libbam)) {                                     # library size from another BAM (e.g. the filtered host BAM)
    l <- fread(cmd = sprintf("samtools idxstats '%s'", libbam), header = FALSE, col.names = c("chr", "len", "mapped", "unmapped"))[is_host(chr) & len > 0]
    if (!nrow(l)) die(basename(libbam), ": no host nuclear chromosomes to take the library size from")
    nm <- intersect(names_for, h$chr)
    if (!length(nm)) die(basename(bam), " has no chromosome named ", paste(names_for, collapse = " / "))
    return(list(name = nm[1], len = h[chr == nm[1], len], mapped = sum(as.numeric(l$mapped)), genome = sum(as.numeric(l$len))))
  }
  nm <- intersect(names_for, h$chr)
  if (!length(nm)) die(basename(bam), " has no chromosome named ", paste(names_for, collapse = " / "),
                       " (names in the BAM look like: ", paste(head(h$chr, 3), collapse = ", "), "; supply --chrom-alias)")
  hh <- h[is_host(chr) & len > 0]
  if (!nrow(hh)) hh <- h[len > 0]          # no recognisable host chromosomes: fall back to everything
  list(name = nm[1], len = h[chr == nm[1], len], mapped = sum(as.numeric(hh$mapped)), genome = sum(as.numeric(hh$len)))
}
info <- lapply(seq_len(nrow(tr)), function(i) if (tr$kind[i] == "bam") bam_chrom(tr$file[i], tr$libfile[i]) else NULL)
clen <- if (any(tr$kind == "bam")) info[[which(tr$kind == "bam")[1]]]$len else NA
if (is.na(clen) && !is.null(cr)) { g <- fread(file.path(cr, "work", "annot", "genome.sizes"), header = FALSE); clen <- g[V1 == chrom, V2] }
reg <- opt("--region")
if (!is.null(reg)) {
  rr <- num(strsplit(reg, "-")[[1]]); if (length(rr) != 2 || any(is.na(rr)) || rr[2] <= rr[1]) die("--region must look like 120e6-150e6")
  rs <- max(0, floor(rr[1])); re <- ceiling(rr[2]); if (length(clen) && !is.na(clen)) re <- min(re, clen)
} else {
  if (!length(clen) || is.na(clen)) die("chromosome length unknown: give a BAM track, --crossrun, or --region")
  rs <- 0; re <- clen
}
## --px-per-bin without --bin: the finest bins that still give each bin N pixels in a figure of the maximum PNG width.
## The space outside the panels is estimated here (3 in, about twice what axis text and legend need) and measured at saving.
label_pos <- tolower(opt("--label-pos", "left")); if (label_pos == "centre") label_pos <- "center"
if (!label_pos %in% c("left", "center", "right")) die("--label-pos must be left, center or right")
MAX_PX <- 32767; OVER_EST_IN <- 3
dpi <- as.numeric(opt("--dpi", 200)); if (is.na(dpi) || dpi <= 0) die("--dpi must be a positive number")
ppb <- opt("--px-per-bin")
if (!is.null(ppb)) { ppb <- as.numeric(ppb); if (is.na(ppb) || ppb <= 0) die("--px-per-bin must be a positive number") }
auto_max <- !is.null(ppb) && is.null(opt("--bin"))
nice <- function(x) { p <- 10^floor(log10(x)); m <- x / p; p * c(1, 2, 2.5, 5, 10)[which.min(abs(c(1, 2, 2.5, 5, 10) - m))] }
bin <- if (auto_max) {
  panel_px <- floor(MAX_PX - OVER_EST_IN * dpi)
  if (panel_px < 100) die(sprintf("--dpi %g leaves no room for the panels in a %d-pixel PNG; lower --dpi", dpi, MAX_PX))
  b <- max(10, ceiling(ppb * (re - rs) / panel_px))
  if (b == 10) say("bins at the 10 bp floor: the figure will be narrower than the maximum width")
  b
} else as.numeric(opt("--bin", nice(max(100, (re - rs) / 4000))))
starts <- seq(rs, re - 1, by = bin)
bins <- data.table(start = starts, end = pmin(starts + bin, re))
say(sprintf("%s:%s-%s, %d bins of %s bp", chrom, format(rs, big.mark = ",", scientific = FALSE), format(re, big.mark = ",", scientific = FALSE), nrow(bins), format(bin, big.mark = ",", scientific = FALSE)))

## ------------------------------------------------------------------ coverage
read_len <- function(bam, nm) {
  x <- suppressWarnings(system(sprintf("samtools view -F 0x904 '%s' '%s' 2>/dev/null | head -n 20000 | awk '{s += length($10)} END {if (NR) print s / NR}'", bam, nm), intern = TRUE))
  if (!length(x) || !nzchar(x)) 100 else as.numeric(x)
}
cache_dir <- opt("--cache"); if (!is.null(cache_dir)) dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
bam_cov <- function(bam, inf) {
  cf <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(substr(digest_str(paste(normalizePath(bam), inf$name, rs, re, bin, mapq)), 1, 20), ".tsv")) else NULL
  if (!is.null(cf) && file.exists(cf) && file.mtime(cf) > file.mtime(bam)) return(fread(cf))
  bf <- tempfile(tmpdir = tmp, fileext = ".bed")
  fwrite(data.table(inf$name, bins$start, bins$end), bf, sep = "\t", col.names = FALSE, scipen = 100)
  x <- fread(cmd = sprintf("samtools bedcov -c -Q %d '%s' '%s'", mapq, bf, bam), header = FALSE)
  out <- data.table(start = x$V2, depth = x$V4 / (x$V3 - x$V2), reads = x$V5)
  if (!is.null(cf)) { tf <- paste0(cf, ".", Sys.getpid()); fwrite(out, tf, sep = "\t", scipen = 100); file.rename(tf, cf) }
  out
}
digest_str <- function(x) { f <- tempfile(tmpdir = tmp); writeLines(x, f); unname(tools::md5sum(f)) }
bw_intervals <- function(f, kind) {
  out <- file.path(tmp, paste0(basename(f), ".bg"))
  if (kind == "bigwig") {
    ok <- FALSE
    for (nm in names_for) {
      if (nzchar(Sys.which("bigWigToBedGraph"))) {
        system2("bigWigToBedGraph", c(sprintf("-chrom=%s", nm), sprintf("-start=%d", as.integer(rs)), sprintf("-end=%d", as.integer(re)), f, out), stderr = FALSE)
      } else {
        py <- sprintf("import pyBigWig,sys\nb=pyBigWig.open(sys.argv[1])\nc=sys.argv[2]\nif c in b.chroms():\n e=min(int(sys.argv[4]),b.chroms()[c])\n for s,t,v in (b.intervals(c,int(sys.argv[3]),e) or []): print(c,s,t,v,sep='\\t')")
        pf <- file.path(tmp, "bw.py"); writeLines(py, pf)
        st <- system2("python3", c(pf, f, nm, as.integer(rs), as.integer(re)), stdout = out, stderr = file.path(tmp, "py.err"))
        if (st != 0) die("cannot read bigWig ", f, ": install bigWigToBedGraph (UCSC tools) or pyBigWig")
      }
      if (file.exists(out) && file.size(out) > 0) { ok <- TRUE; break }
    }
    if (!ok) return(data.table(s = numeric(), e = numeric(), v = numeric()))
    x <- fread(out, header = FALSE)
  } else {
    pat <- paste(sprintf("$1==\"%s\"", names_for), collapse = " || ")
    cat_cmd <- if (grepl("\\.gz$", f)) sprintf("gzip -dc '%s'", f) else sprintf("cat '%s'", f)
    x <- fread(cmd = sprintf("%s | awk -F'\\t' '%s'", cat_cmd, pat), header = FALSE)
    if (!nrow(x)) return(data.table(s = numeric(), e = numeric(), v = numeric()))
  }
  data.table(s = x$V2, e = x$V3, v = x$V4)[e > rs & s < re]
}
bin_intervals <- function(iv) {      # length-weighted mean of the interval values in each bin (gaps count as 0)
  b <- copy(bins)[, `:=`(bs = start, be = end)]
  if (!nrow(iv)) return(b[, .(start, value = 0)])
  setkey(b, bs, be); iv[, `:=`(s = as.numeric(s), e = as.numeric(e))]; setkey(iv, s, e)
  ov <- foverlaps(iv, b, by.x = c("s", "e"), by.y = c("bs", "be"), nomatch = 0L)
  ov[, w := pmin(e, be) - pmax(s, bs)]
  agg <- ov[, .(value = sum(v * w)), by = start]
  m <- merge(bins, agg, by = "start", all.x = TRUE)
  m[, .(start, value = fifelse(is.na(value), 0, value / (end - start)))]
}

say("computing coverage for ", nrow(tr), " track(s)")
cov_one <- function(i) {          # raw coverage of one track (and its control), computed once for every --norm
  t <- tr[i]
  if (t$kind == "bam") {
    inf <- info[[i]]; ip <- bam_cov(t$file, inf)
    res <- list(kind = "bam", label = t$label, ip = ip, inf = inf, gmean = inf$mapped * read_len(t$file, inf$name) / inf$genome)
    if (!is.na(t$control)) { ci <- bam_chrom(t$control, t$libcontrol); res$ct <- bam_cov(t$control, ci); res$ci <- ci
                             res$ct_g <- ci$mapped * read_len(t$control, ci$name) / ci$genome }
    res
  } else {
    b <- bin_intervals(bw_intervals(t$file, t$kind))
    list(kind = "file", label = t$label, b = b)
  }
}
RAW <- parallel::mclapply(seq_len(nrow(tr)), cov_one, mc.cores = max(1, min(threads, nrow(tr))))
bad_r <- vapply(RAW, function(x) inherits(x, "try-error") || is.null(x), TRUE)
if (any(bad_r)) die("coverage failed for: ", paste(tr$label[bad_r], collapse = ", "))
norm_D <- function(norm) {
  D <- rbindlist(lapply(RAW, function(r) {
    if (r$kind != "bam") return(data.table(track = r$label, start = r$b$start, value = r$b$value, control = NA_real_))
    ct <- r$ct
    if (norm == "input")
      return(data.table(track = r$label, start = r$ip$start, value = log2(((r$ip$reads + 1) / r$inf$mapped) / ((ct$reads + 1) / r$ci$mapped)), control = NA_real_))
    sc <- switch(norm, raw = c(1, 1), cpm = c(1e6 / r$inf$mapped, if (!is.null(ct)) 1e6 / r$ci$mapped else 1),
                 genome = c(1 / r$gmean, if (!is.null(ct)) 1 / r$ct_g else 1))
    data.table(track = r$label, start = r$ip$start, value = r$ip$depth * sc[1], control = if (is.null(ct)) NA_real_ else ct$depth * sc[2])
  }))
  if (smooth_k > 1) {
    rm_ <- function(x) as.numeric(stats::filter(x, rep(1 / smooth_k, smooth_k), sides = 2))
    D[, `:=`(value = rm_(value), control = if (all(is.na(control))) control else rm_(control)), by = track]
  }
  D[, track := factor(track, unique(tr$label))]
  D[, `:=`(end = pmin(start + bin, re), mid = (start + pmin(start + bin, re)) / 2)]
  D
}
out_for <- function(norm, base, ext_default = ".png") {
  o <- gsub("{chrom}", CHR, base, fixed = TRUE)
  if (grepl("{norm}", o, fixed = TRUE)) return(gsub("{norm}", norm, o, fixed = TRUE))
  if (length(norms) > 1) sub("(\\.[A-Za-z]+)$", paste0("_", norm, "\\1"), o) else o
}

## ------------------------------------------------------------------ annotations
read_any <- function(f) {
  if (grepl("\\.(bb|bigBed|bigbed)$", f)) {
    need("bigBedToBed"); out <- file.path(tmp, paste0(basename(f), ".bed"))
    for (nm in names_for) { system2("bigBedToBed", c(sprintf("-chrom=%s", nm), f, out)); if (file.size(out) > 0) break }
    return(out)
  }
  f
}
ann_classes <- function(f, kind) {   # -> chrom start end label (canonical chromosome name, region only)
  src <- read_any(f)
  sizes <- file.path(tmp, "sizes.tsv"); fwrite(data.table(chrom, max(re, 1)), sizes, sep = "\t", col.names = FALSE, scipen = 100)
  al <- file.path(tmp, "alias.tsv"); fwrite(alias, al, sep = "\t", col.names = FALSE, scipen = 100)
  grep_names <- paste(sprintf("-e '%s'", names_for), collapse = " ")
  cat_cmd <- if (grepl("\\.gz$", src)) sprintf("gzip -dc '%s'", src) else sprintf("cat '%s'", src)
  out <- file.path(tmp, paste0(kind, ".class.bed"))
  conv <- file.path(bin_root, "features_to_classbed.py")
  if (!file.exists(conv)) die("bin/features_to_classbed.py not found next to this script (keep the bundle layout)")
  st <- system(sprintf("%s | grep -w -F %s | python3 '%s' --kind %s --genome '%s' --alias '%s' - > '%s' 2> '%s.log'",
                       cat_cmd, grep_names, conv, kind, sizes, al, out, out))
  if (st != 0 || !file.exists(out) || file.size(out) == 0) { say("no ", kind, " features on ", chrom, " in ", basename(f), " (if the file uses other chromosome names, e.g. NC_060925.1, add --chrom-alias or --crossrun)"); return(NULL) }
  x <- fread(out, header = FALSE, col.names = c("chrom", "start", "end", "label"))
  x[end > rs & start < re][, `:=`(start = pmax(start, rs), end = pmin(end, re))][]
}
censat_f <- opt("--censat") %||% (if (!is.null(cr) && file.exists(file.path(cr, "work/annot/censat.class.bed"))) file.path(cr, "work/annot/censat.class.bed"))
CS <- if (!is.null(censat_f)) ann_classes(censat_f, "censat")
BD <- NULL
if (!is.null(opt("--bands"))) {
  bf <- read_any(opt("--bands"))
  BD <- fread(cmd = if (grepl("\\.gz$", bf)) sprintf("gzip -dc '%s'", bf) else sprintf("cat '%s'", bf), header = FALSE, fill = TRUE)
  if (ncol(BD) < 5) die("--bands needs 5 columns: chrom start end name gieStain")
  BD <- BD[, .(chrom = canon(V1), start = V2, end = V3, name = V4, stain = V5)][chrom == CHR & end > rs & start < re]
  if (!nrow(BD)) { say("no bands for ", chrom, " in ", opt("--bands")); BD <- NULL }
}

## ------------------------------------------------------------------ plot
u <- if (re - rs < 2e6) 1e3 else 1e6; ulab <- if (u == 1e3) "kb" else "Mb"
xlim <- c(rs, re) / u
## x-axis breaks: round numbers inside the region plus both ends of the region, always labelled. Inside breaks closer
## than 7% of the width to an end are dropped so their labels cannot collide with the end labels. End labels carry up to
## two decimals of the unit (e.g. 248.39 Mb for the end of chr1), trailing zeros dropped.
x_breaks <- function(lims) {
  inner <- pretty(lims, n = 5); inner <- inner[inner > lims[1] & inner < lims[2]]
  inner <- inner[abs(inner - lims[1]) > 0.07 * diff(lims) & abs(inner - lims[2]) > 0.07 * diff(lims)]
  c(lims[1], inner, lims[2])
}
x_labels <- function(b) {
  b <- b[!is.na(b)]; ends <- b == min(b) | b == max(b); out <- character(length(b))
  out[ends] <- formatC(round(b[ends], 2), format = "f", digits = 2, big.mark = ",", drop0trailing = TRUE)
  if (any(!ends)) out[!ends] <- format(b[!ends], big.mark = ",", trim = TRUE)
  out
}
xs <- scale_x_continuous(limits = xlim, expand = c(0, 0), breaks = x_breaks, labels = x_labels)
blank_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
## right margin wide enough for half of the right-end label, which is centred on the panel edge (axis text 8 pt)
end_half_pt <- nchar(x_labels(x_breaks(xlim))[length(x_breaks(xlim))]) * 8 * 0.6 / 2
base <- theme_bw(10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                             plot.margin = margin(2, max(8, ceiling(end_half_pt) + 2), 2, 4), legend.position = "right")
pal <- c("#2a78d6", "#d1495b", "#00798c", "#edae49", "#66a182", "#8d5a97", "#30638e", "#b5651d")
peaks <- rbindlist(lapply(seq_len(nrow(tr)), function(i) {
  f <- tr$peaks[i]; if (is.na(f) || !file.exists(f)) return(NULL)
  x <- fread(cmd = if (grepl("\\.gz$", f)) sprintf("gzip -dc '%s'", f) else sprintf("cat '%s'", f), header = FALSE, select = 1:3)
  x <- x[canon(V1) == chrom & V3 > rs & V2 < re]
  if (nrow(x)) data.table(track = tr$label[i], start = pmax(x$V2, rs), end = pmin(x$V3, re), mid = (x$V2 + x$V3) / 2) else NULL
}), fill = TRUE)
if (nrow(peaks)) peaks[, track := factor(track, unique(tr$label))]

## ---- top-N peaks per track, from the original narrowPeak/broadPeak file (needs the score and summit columns).
## Peak files are 0-based half-open (BED-like); the narrowPeak summit (column 10) is a 0-based offset from the start,
## so summit = start + offset is a 0-based position. Without a summit, the midpoint is used.
n_top <- as.integer(opt("--label-top", 0)); if (is.na(n_top) || n_top < 0) die("--label-top must be a whole number")
top_by <- opt("--top-by", "score"); top_col <- c(score = 5L, signal = 7L, p = 8L, q = 9L)[top_by]
if (is.na(top_col)) die("--top-by must be score, signal, p or q")
TOP <- data.table()
if (n_top > 0) {
  TOP <- rbindlist(lapply(seq_len(nrow(tr)), function(i) {
    f <- tr$rawpeaks[i]
    if (is.na(f) || !file.exists(f)) { say("--label-top: no peak file for ", tr$label[i], "; not labelled"); return(NULL) }
    x <- fread(cmd = if (grepl("\\.gz$", f)) sprintf("gzip -dc '%s'", f) else sprintf("cat '%s'", f), header = FALSE, fill = TRUE)
    x <- x[!grepl("^(track|browser|#)", V1)]
    if (ncol(x) < top_col) die("--label-top: ", f, " has no column ", top_col, " (", top_by, "); use --top-by with a column it has")
    nm <- function(k) if (ncol(x) >= k) suppressWarnings(as.numeric(x[[k]])) else rep(NA_real_, nrow(x))
    d <- data.table(track = tr$label[i], chrom = canon(as.character(x$V1)), start = as.numeric(x$V2), end = as.numeric(x$V3),
                    name = if (ncol(x) >= 4) as.character(x$V4) else NA_character_,
                    score = nm(5), signal = nm(7), neglog10_p = nm(8), neglog10_q = nm(9), off = nm(10))
    d[, summit := fifelse(!is.na(off) & off >= 0, start + off, floor((start + end) / 2))][, off := NULL]
    d <- d[chrom == CHR & summit >= rs & summit < re]
    d[, key := get(c(score = "score", signal = "signal", p = "neglog10_p", q = "neglog10_q")[[top_by]])]
    d <- d[!is.na(key)]
    if (!nrow(d)) { say("--label-top: no peaks of ", tr$label[i], " in the region"); return(NULL) }
    setorder(d, -key, -signal, start)                  # ties: higher signal first, then position
    d <- head(d, n_top)[, rank := seq_len(.N)][, key := NULL]
    d
  }), fill = TRUE)
  if (nrow(TOP)) {
    tf <- opt("--top-out") %||% sub("(\\.[A-Za-z]+)?$", "_top_peaks.tsv", gsub("{norm}", "", gsub("{chrom}", CHR, opt("--out", sprintf("coverage_%s.png", chrom)), fixed = TRUE), fixed = TRUE))
    dir.create(dirname(tf), showWarnings = FALSE, recursive = TRUE)
    fwrite(TOP[, .(track, rank, chrom, start, end, summit, name, score, signal, neglog10_p, neglog10_q, ranked_by = top_by)],
           tf, sep = "\t", na = "NA", quote = FALSE, scipen = 100)
    say(sprintf("top %d peaks by %s: %s (start/end 0-based half-open; summit 0-based)", n_top, top_by, tf))
    TOP[, track := factor(track, levels(peaks$track) %||% unique(tr$label))]
  }
}

for (norm in norms) {
D <- norm_D(norm)
tb <- opt("--table")
if (!is.null(tb)) { tf <- out_for(norm, tb); fwrite(D[, .(track, chrom = CHR, start, end = pmin(start + bin, re), value, control)], tf, sep = "\t", na = "NA", quote = FALSE, scipen = 100); say("wrote ", tf) }
ylab <- switch(norm, raw = "Mean read depth", cpm = "Depth per million host reads", genome = "Fold over host genome-wide mean",
               input = "log2(IP / control)")
if (all(tr$kind != "bam")) ylab <- "Signal (as in the file)" else if (any(tr$kind != "bam")) ylab <- paste(ylab, "/ file signal")
ylog <- has("--ylog") && norm != "input"
cap <- as.numeric(opt("--ymax", NA))
if (is.na(cap) && !is.null(opt("--cap"))) cap <- quantile(c(D$value, D$control), as.numeric(opt("--cap")), na.rm = TRUE)
ymax_all <- max(if (!is.na(cap)) cap else c(D$value, D$control), na.rm = TRUE)

## approximate size of one coverage panel in inches (for placing rank labels): figure size minus titles, axes and legend,
## shared among the panels by their relative heights
n_tr <- nlevels(D$track); n_st <- if (nrow(peaks)) sum(levels(D$track) %in% peaks$track) else 0
n_lab <- if (nrow(TOP)) length(intersect(levels(D$track), as.character(unique(TOP$track)))) else 0
h_est <- as.numeric(opt("--height", 1.2 + 1.5 * n_tr + 0.15 * n_st + 0.3 * n_lab + 0.6 * ((!is.null(CS) && nrow(CS) > 0) + !is.null(BD))))
panel_h_in <- max(0.3, (h_est - 1.3) / (n_tr + 0.07 * n_st + 0.22 * (!is.null(CS) && nrow(CS) > 0) + 0.2 * !is.null(BD)))
panel_w_in <- if (!is.null(ppb) && !auto_max) nrow(bins) * ppb / dpi else if (auto_max) (MAX_PX - 1) / dpi - 1.9 else as.numeric(opt("--width", 12)) - 1.9
track_plot <- function(k, i) {
  d <- D[track == k]; col <- pal[(i - 1) %% length(pal) + 1]
  ## a triangle marks a bin whose mean (the bar's height) is above the ceiling; bins are means, never single bases
  if (!is.na(cap)) d[, `:=`(clip = !is.na(value) & value > cap, value = pmin(value, cap), control = pmin(control, cap))] else d[, clip := FALSE]
  p <- ggplot(d, aes(mid / u))
  if (norm == "input") {
    p <- p + geom_col(aes(y = value, fill = value > 0), width = bin / u, show.legend = FALSE, na.rm = TRUE) +
      scale_fill_manual(values = c(`TRUE` = col, `FALSE` = "grey60")) + geom_hline(yintercept = 0, linewidth = 0.3)
  } else {
    ## one flat bar per bin at the bin's mean, so what is drawn is exactly what is compared with the ceiling
    if (any(!is.na(d$control)) && !has("--no-control-track"))
      p <- p + geom_rect(data = d[!is.na(control)], aes(xmin = start / u, xmax = end / u, ymin = 0, ymax = control), fill = "grey55", alpha = 0.55)
    p <- p + geom_rect(data = d[!is.na(value)], aes(xmin = start / u, xmax = end / u, ymin = 0, ymax = value), fill = col, alpha = 0.85)
  }
  if (any(d$clip)) p <- p + geom_point(data = d[clip == TRUE], aes(x = mid / u, y = cap), shape = 25, size = 1.1, fill = "black", inherit.aes = FALSE)
  ysc_args <- if (ylog) list(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000, 1e4, 3e4)) else list()
  tf <- if (ylog) log1p else identity; inv <- if (ylog) expm1 else identity
  ylim_fix <- NULL
  ysc <- do.call(scale_y_continuous, c(ysc_args, list(expand = if (is.null(ylim_fix)) expansion(mult = 0.05) else c(0, 0))))
  p <- p + xs + ysc + labs(y = NULL, title = NULL) + annotate("label", x = switch(label_pos, left = xlim[1], center = mean(xlim), right = xlim[2]), y = Inf, label = k,
                                      hjust = switch(label_pos, left = 0, center = 0.5, right = 1), vjust = 1, size = 3.2, fontface = "bold", fill = "white", label.size = 0, alpha = 0.85) + base + blank_x
  if (!is.null(ylim_fix)) p <- p + coord_cartesian(ylim = ylim_fix, expand = FALSE) else
  if (has("--same-y")) p <- p + coord_cartesian(ylim = c(if (norm == "input") min(D$value, na.rm = TRUE) else 0, ymax_all))
  p
}
## Rank labels (--label-top): placed along x as close to their summits as possible without touching, by a 1-D
## pool-adjacent-violators spread of label widths (inches -> bp from the panel width estimate)
rank_x <- function(tp) {
  bp_per_in <- diff(xlim) * u / panel_w_in
  tp[, wd := (nchar(as.character(rank)) * 0.078 + 0.07) * bp_per_in]            # a 9.7-pt bold label plus padding and a gap
  setorder(tp, summit)
  nn <- nrow(tp); off <- c(0, cumsum((head(tp$wd, -1) + tail(tp$wd, -1)) / 2)); blk <- list()
  for (q in seq_len(nn)) {
    blk[[length(blk) + 1]] <- c(q, q)
    repeat {
      nb <- length(blk); cur <- blk[[nb]]; ix <- cur[1]:cur[2]
      pos <- mean(tp$summit[ix] - off[ix]) + off[ix]
      if (nb == 1) break
      pv <- blk[[nb - 1]]; px <- pv[1]:pv[2]; ppos <- mean(tp$summit[px] - off[px]) + off[px]
      if (tail(ppos, 1) + (tp$wd[pv[2]] + tp$wd[cur[1]]) / 2 <= pos[1] + 1e-9) break
      blk[[nb - 1]] <- c(pv[1], cur[2]); blk[[nb]] <- NULL
    }
  }
  lx <- numeric(nn); for (bb in blk) { ix <- bb[1]:bb[2]; lx[ix] <- mean(tp$summit[ix] - off[ix]) + off[ix] }
  lim0 <- xlim[1] * u + tp$wd[1] / 2; lim1 <- xlim[2] * u - tp$wd[nn] / 2
  if (lx[1] < lim0) lx <- lx + (lim0 - lx[1]); if (lx[nn] > lim1) lx <- lx - (lx[nn] - lim1)
  if (lx[1] < lim0 - 1e-6) say("WARNING: the rank labels are wider than the panel and overlap; label fewer peaks or widen the figure")
  tp[, lx := lx][]
}
## peaks: a thin strip of their own under each track (never drawn over the coverage); each peak at its true extent,
## with a black outline so a peak narrower than a pixel still shows as a tick about 0.2 mm wide. With --label-top the
## strip grows downwards: the ticks keep their height at the top and the rank labels sit in a band below them, each
## joined to its peak's summit by a thin line that drops straight down from the tick and bends only below all ticks.
TICK_IN <- 0.07; LAB_IN <- 0.30                                                  # tick row and label band, inches
strip_h <- function(k) if (nrow(TOP) && any(TOP$track == k)) TICK_IN + LAB_IN / panel_h_in else TICK_IN   # relative to a track
peak_strip <- function(k) {
  pk <- peaks[track == k]
  lab <- nrow(TOP) && any(TOP$track == k)
  ft <- if (lab) (TICK_IN * panel_h_in) / (TICK_IN * panel_h_in + LAB_IN) else 1    # share of the strip taken by the ticks
  y0 <- 1 - ft
  p <- ggplot(pk) + geom_rect(aes(xmin = start / u, xmax = end / u, ymin = y0, ymax = 1), fill = "black", colour = "black", linewidth = 0.3, linejoin = "mitre")
  if (lab) {
    tp <- rank_x(copy(TOP[track == k]))
    bend <- y0 - 0.12 * (1 - ft); ly <- y0 - 0.30 * (1 - ft)                         # leader bend, top of the labels
    path <- rbind(tp[, .(rank, ord = 1, x = summit, yy = y0)], tp[, .(rank, ord = 2, x = summit, yy = bend)], tp[, .(rank, ord = 3, x = lx, yy = ly)])
    setorder(path, rank, ord)
    p <- p + geom_path(data = path, aes(x = x / u, y = yy, group = rank), inherit.aes = FALSE, linewidth = 0.25, colour = "grey35") +
      geom_label(data = tp, aes(x = lx / u, y = ly, label = rank), inherit.aes = FALSE, vjust = 1, size = 3.4, fontface = "bold",
                 label.size = 0, label.padding = unit(0.06, "lines"), fill = "white")
  }
  p + xs + scale_y_continuous(limits = c(0, 1), breaks = NULL, expand = c(0, 0)) + labs(y = "peaks") + base + blank_x +
    coord_cartesian(clip = "off") +
    theme(axis.title.y = element_text(angle = 0, vjust = if (lab) 1 - ft / 2 else 0.5, size = 7, colour = "black"), panel.border = element_blank(),
          plot.margin = margin(0, 8, 4, 4))
}
plots <- list(); heights <- numeric(); n_strip <- 0
for (i in seq_along(levels(D$track))) {
  k <- levels(D$track)[i]
  plots <- c(plots, list(track_plot(k, i))); heights <- c(heights, 1)
  if (nrow(peaks) && any(peaks$track == k)) { plots <- c(plots, list(peak_strip(k))); heights <- c(heights, strip_h(k)); n_strip <- n_strip + 1 }
}

## CenSat colours: distinct hues and lightness so neighbouring arrays stay apart (colour-blind-safe Okabe-Ito base)
CS_COL <- c(active_hor = "#B2182B", hor = "#E69F00", dhor = "#F0E442", mon = "#CC79A7",
            hsat1A = "#009E73", hsat1B = "#56B4E9", hsat1 = "#009E73", hsat2 = "#5B9BD5", hsat3 = "#000000",
            bsat = "#7B3294", gsat = "#8C510A", censat = "#7F7F7F", ct = "#CFCFCF", rDNA = "#B2DF8A")
if (!is.null(CS) && nrow(CS)) {
  ## a class always gets the same colour, whatever else is in the plot: known classes from CS_COL (matched without regard
  ## to case), any other class from a fixed list picked by a checksum of its name (never by the order classes appear)
  key <- tolower(names(CS_COL))
  CS[, label := { m <- match(tolower(label), key); ifelse(is.na(m), label, names(CS_COL)[m]) }]
  EXTRA <- c("#1F78B4", "#33A02C", "#FB9A99", "#FDBF6F", "#CAB2D6", "#6A3D9A", "#B15928", "#A6CEE3", "#E31A1C", "#FF7F00")
  name_col <- function(x) EXTRA[sum(utf8ToInt(x) * seq_along(utf8ToInt(x))) %% length(EXTRA) + 1]
  lv <- unique(CS$label)
  lv <- c(intersect(names(CS_COL), lv), sort(setdiff(lv, names(CS_COL))))      # legend in a fixed order
  cols <- setNames(vapply(lv, function(x) if (x %in% names(CS_COL)) CS_COL[[x]] else name_col(x), ""), lv)
  CS[, label := factor(label, lv)]
  p <- ggplot(CS) + geom_rect(aes(xmin = start / u, xmax = end / u, ymin = 0, ymax = 1, fill = label)) +
    scale_fill_manual(values = cols, name = "CenSat") + geom_rect(aes(xmin = start / u, xmax = end / u, ymin = 0, ymax = 1), fill = NA, colour = "white", linewidth = 0.1) + xs + scale_y_continuous(breaks = NULL) + labs(y = "CenSat") +
    base + blank_x + theme(axis.title.y = element_text(angle = 0, vjust = 0.5, size = 8), legend.key.size = unit(0.35, "cm"),
                           legend.text = element_text(size = 7), legend.title = element_text(size = 8)) +
    guides(fill = guide_legend(ncol = 1))
  plots <- c(plots, list(p)); heights <- c(heights, 0.22)
}
if (!is.null(BD)) {
  STAIN <- c(gneg = "#ffffff", gpos25 = "#d9d9d9", gpos33 = "#c8c8c8", gpos50 = "#a6a6a6", gpos66 = "#8c8c8c", gpos75 = "#737373",
             gpos100 = "#1f1f1f", gpos = "#525252", acen = "#c0392b", gvar = "#b8c4d6", stalk = "#7f9cc0")
  BD[, fill := STAIN[stain]]; BD[is.na(fill), fill := "#bdbdbd"]
  p <- ggplot(BD) + geom_rect(aes(xmin = pmax(start, rs) / u, xmax = pmin(end, re) / u, ymin = 0, ymax = 1), fill = BD$fill,
                              colour = "grey30", linewidth = 0.15) +
    xs + scale_y_continuous(breaks = NULL) + labs(y = "Bands") + base +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5, size = 8))
  dark <- BD$stain %in% c("gpos75", "gpos100", "gpos", "acen", "gpos66")
  if (nrow(BD) <= 60) p <- p + geom_text(aes(x = (pmax(start, rs) + pmin(end, re)) / (2 * u), y = 0.5, label = name), size = 2.2,
                                         colour = ifelse(dark, "white", "grey15"), check_overlap = TRUE)
  plots <- c(plots, list(p)); heights <- c(heights, 0.2)
}
## the last panel carries the x-axis
last <- length(plots)
plots[[last]] <- plots[[last]] + theme(axis.text.x = element_text(), axis.ticks.x = element_line(), axis.title.x = element_text()) +
  labs(x = sprintf("%s position (%s)", chrom, ulab))
sub <- c(sprintf("bins of %s bp", format(bin, big.mark = ",")), if (mapq > 0) sprintf("MAPQ >= %d", mapq),
         if (smooth_k > 1) sprintf("smoothed over %d bins", smooth_k), if (!is.na(cap)) sprintf("y clipped at %.3g (triangle = bin mean above it)", cap),
         if (any(!is.na(D$control)) && !has("--no-control-track")) "grey = control",
         if (nrow(peaks)) "black ticks under each track = peaks")
title <- opt("--title", sprintf("%s coverage", chrom))
fig <- wrap_plots(plots, ncol = 1, heights = heights, guides = "collect") +
  plot_annotation(title = title, subtitle = paste(c(ylab, sub), collapse = "; "),
                  theme = theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9)))
out <- out_for(norm, opt("--out", sprintf("coverage_%s.png", chrom)))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
w <- as.numeric(opt("--width", 12)); h <- as.numeric(opt("--height", 1.2 + 1.5 * length(levels(D$track)) + 0.15 * n_strip + 0.3 * n_lab + 0.6 * (last - length(levels(D$track)) - n_strip)))
if (!is.null(ppb)) {
  ## width outside the panels (y-axis text, labels, legend) measured from the assembled figure, then panels sized to the bins
  gt <- patchwork::patchworkGrob(fig)
  fixed <- vapply(seq_along(gt$widths), function(i) !identical(grid::unitType(gt$widths[i]), "null"), TRUE)
  png_dev <- grepl("\\.png$", out, ignore.case = TRUE)
  if (png_dev) { grDevices::png(tempfile(tmpdir = tmp, fileext = ".png"), width = 1000, height = 1000, res = dpi) } else grDevices::pdf(NULL)
  over <- sum(grid::convertWidth(gt$widths[fixed], "in", valueOnly = TRUE)); grDevices::dev.off()
  n_bins <- nrow(bins)
  w <- over + ceiling(n_bins * ppb) / dpi + 0.05
  if (auto_max && bin > 10) {
    if (over > OVER_EST_IN) say(sprintf("WARNING: axis text and legend take %.1f in, more than the %.1f in allowed for; bins get slightly under %g px", over, OVER_EST_IN, ppb))
    w <- (MAX_PX - 1) / dpi                            # full width: the panels take everything the margins leave
    if (over > OVER_EST_IN) w <- min(w, over + ceiling(n_bins * ppb) / dpi + 0.05)
  }
  say(sprintf("--px-per-bin %g: %s bins of %s bp -> figure %.1f in wide (%s px at %g dpi; %.2f px per bin)", ppb, format(n_bins, big.mark = ","),
              format(bin, big.mark = ","), w, format(floor(w * dpi), big.mark = ","), dpi, (w - over) * dpi / n_bins))
  if (png_dev && w * dpi > 32767)
    die(sprintf("a PNG this wide (%s px) is beyond the 32,767-pixel limit of the graphics device. Use --bin %s or larger, a --region, or write .pdf/.svg",
                format(ceiling(w * dpi), big.mark = ","), format(ceiling(ppb * (re - rs) / ((32767 / dpi - over) * dpi)), big.mark = ",")))
}
ggsave(out, fig, width = w, height = h, dpi = dpi, limitsize = FALSE, bg = "white")
say("wrote ", out)
}
