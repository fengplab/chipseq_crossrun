#!/usr/bin/env Rscript
## Read coverage along a whole chromosome (or part of one), for one or more ChIP libraries, with optional
## annotation tracks underneath: chromosome bands (ideogram), CenSat arrays and RepeatMasker class density.
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
##                              and (unless given) the CenSat/RepeatMasker annotation and chromosome aliases
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
##   --rmsk FILE                RepeatMasker (BED from bigBedToBed, rmsk.txt, .out, or class BED4)
##   --rmsk-classes A,B         RepeatMasker classes to show (default: the 6 covering most of the region)
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
##   --cache DIR                keep the binned BAM counts here and reuse them (same BAM, contig, bins and MAPQ)
##   --table FILE               also write the binned values (TSV)
##   --threads N                parallel samtools jobs (default 4)
##
## Examples:
##   Rscript R/plot_chrom_coverage.R --chrom chr9 --crossrun crossrun_out_human_broad \
##       --samples THP1_HSV1_8h_IP-cGAS,THP1_IAV_6h_IP-cGAS --norm genome --bands chm13v2.0_cytobands_allchrs.bed
##   Rscript R/plot_chrom_coverage.R --chrom chr1 --bam ip.bam --control input.bam --norm input \
##       --censat censat.bb --rmsk rmsk.bb --region 120e6-150e6
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
           "--norm", "--mapq", "--region", "--bin", "--smooth", "--bands", "--censat", "--rmsk", "--rmsk-classes",
           "--chrom-alias", "--ymax", "--cap", "--title", "--out", "--width", "--height", "--table", "--threads", "--cache")
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
                 libfile = character(), libcontrol = character())
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
                               libfile = r$bam_file, libcontrol = r$control_bam %||% NA_character_))
  } else
  tr <- rbind(tr, data.table(kind = "bam", file = r$bam_file, control = r$control_bam %||% NA_character_, label = ss,
                             peaks = ifelse(file.exists(pk), pk, NA_character_), libfile = NA_character_, libcontrol = NA_character_))
}
add <- function(kind, files) if (length(files)) tr <<- rbind(tr, data.table(kind = kind, file = files, control = NA_character_,
                                                                           label = sub("\\.mLb\\.clN(\\.sorted)?$", "", sub("\\.(bam|bw|bigwig|bigWig|bedGraph|bedgraph|bg)(\\.gz)?$", "", basename(files))),
                                                                           peaks = NA_character_, libfile = NA_character_, libcontrol = NA_character_))
add("bam", csv(opt("--bam"))); add("bigwig", csv(opt("--bigwig"))); add("bedgraph", csv(opt("--bedgraph")))
if (!nrow(tr)) die("no tracks: give --bam, --bigwig, --bedgraph or --crossrun with --samples")
ctl <- csv(opt("--control"))
if (length(ctl)) {
  ib <- which(tr$kind == "bam" & is.na(tr$control))
  if (!(length(ctl) %in% c(1, length(ib)))) die("--control needs one file, or one per --bam (", length(ib), ")")
  tr$control[ib] <- rep_len(ctl, length(ib))
}
lab <- csv(opt("--labels")); if (length(lab)) { if (length(lab) != nrow(tr)) die("--labels needs ", nrow(tr), " names"); tr$label <- lab }
pkf <- csv(opt("--peaks")); if (length(pkf)) { if (length(pkf) != nrow(tr)) die("--peaks needs one file per track (", nrow(tr), ")"); tr$peaks <- pkf }
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
nice <- function(x) { p <- 10^floor(log10(x)); m <- x / p; p * c(1, 2, 2.5, 5, 10)[which.min(abs(c(1, 2, 2.5, 5, 10) - m))] }
bin <- as.numeric(opt("--bin", nice(max(100, (re - rs) / 4000))))
starts <- seq(rs, re - 1, by = bin)
bins <- data.table(start = starts, end = pmin(starts + bin, re))
say(sprintf("%s:%s-%s, %d bins of %s bp", chrom, format(rs, big.mark = ","), format(re, big.mark = ","), nrow(bins), format(bin, big.mark = ",")))

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
  D[, mid := (start + pmin(start + bin, re)) / 2]
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
rmsk_f <- opt("--rmsk") %||% (if (!is.null(cr) && file.exists(file.path(cr, "work/annot/repeatmasker.class.bed"))) file.path(cr, "work/annot/repeatmasker.class.bed"))
CS <- if (!is.null(censat_f)) ann_classes(censat_f, "censat")
RM <- if (!is.null(rmsk_f)) ann_classes(rmsk_f, "repeatmasker")
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
xs <- scale_x_continuous(limits = xlim, expand = c(0, 0), labels = function(b) format(b, big.mark = ","))
blank_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
base <- theme_bw(10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
                             plot.margin = margin(2, 8, 2, 4), legend.position = "right")
pal <- c("#2a78d6", "#d1495b", "#00798c", "#edae49", "#66a182", "#8d5a97", "#30638e", "#b5651d")
peaks <- rbindlist(lapply(seq_len(nrow(tr)), function(i) {
  f <- tr$peaks[i]; if (is.na(f) || !file.exists(f)) return(NULL)
  x <- fread(cmd = if (grepl("\\.gz$", f)) sprintf("gzip -dc '%s'", f) else sprintf("cat '%s'", f), header = FALSE, select = 1:3)
  x <- x[canon(V1) == chrom & V3 > rs & V2 < re]
  if (nrow(x)) data.table(track = tr$label[i], mid = (x$V2 + x$V3) / 2) else NULL
}), fill = TRUE)
if (nrow(peaks)) peaks[, track := factor(track, unique(tr$label))]

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

track_plot <- function(k, i) {
  d <- D[track == k]; col <- pal[(i - 1) %% length(pal) + 1]
  if (!is.na(cap)) d[, `:=`(clip = value > cap, value = pmin(value, cap), control = pmin(control, cap))] else d[, clip := FALSE]
  p <- ggplot(d, aes(mid / u))
  if (norm == "input") {
    p <- p + geom_col(aes(y = value, fill = value > 0), width = bin / u, show.legend = FALSE, na.rm = TRUE) +
      scale_fill_manual(values = c(`TRUE` = col, `FALSE` = "grey60")) + geom_hline(yintercept = 0, linewidth = 0.3)
  } else {
    if (any(!is.na(d$control)) && !has("--no-control-track"))
      p <- p + geom_area(aes(y = control), fill = "grey55", alpha = 0.55, outline.type = "upper", na.rm = TRUE)
    p <- p + geom_area(aes(y = value), fill = col, alpha = 0.85, na.rm = TRUE)
  }
  if (any(d$clip)) p <- p + geom_point(data = d[clip == TRUE], aes(y = value), shape = 25, size = 1.1, fill = "black")
  pk <- if (nrow(peaks)) peaks[track == k] else peaks
  if (nrow(pk)) p <- p + geom_rug(data = pk, aes(x = mid / u), inherit.aes = FALSE, sides = "b", colour = "#e07b00",
                                  length = unit(0.06, "npc"), linewidth = 0.25)
  ysc <- if (ylog) scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000, 1e4, 3e4)) else scale_y_continuous()
  p <- p + xs + ysc + labs(y = NULL, title = NULL) + annotate("label", x = xlim[1], y = Inf, label = k, hjust = 0, vjust = 1, size = 3.2, fontface = "bold", fill = "white", label.size = 0, alpha = 0.85) + base + blank_x
  if (has("--same-y")) p <- p + coord_cartesian(ylim = c(if (norm == "input") min(D$value, na.rm = TRUE) else 0, ymax_all))
  p
}
plots <- lapply(seq_along(levels(D$track)), function(i) track_plot(levels(D$track)[i], i))
heights <- rep(1, length(plots))

CS_COL <- c(active_hor = "#c0392b", hor = "#e8847a", dhor = "#f4b183", mon = "#f9d9a8", hsat1A = "#1b9e77", hsat1B = "#7fc97f",
            hsat1 = "#1b9e77", hsat2 = "#386cb0", hsat3 = "#2c3e8f", bsat = "#984ea3", gsat = "#e7298a", censat = "#999999",
            ct = "#d9d9d9", rDNA = "#a6761d")
if (!is.null(CS) && nrow(CS)) {
  lv <- unique(CS$label); cols <- CS_COL[lv]; cols[is.na(cols)] <- grDevices::hcl.colors(sum(is.na(cols)), "Dark 3")
  names(cols) <- lv
  p <- ggplot(CS) + geom_rect(aes(xmin = start / u, xmax = end / u, ymin = 0, ymax = 1, fill = label)) +
    scale_fill_manual(values = cols, name = "CenSat") + xs + scale_y_continuous(breaks = NULL) + labs(y = "CenSat") +
    base + blank_x + theme(axis.title.y = element_text(angle = 0, vjust = 0.5, size = 8), legend.key.size = unit(0.35, "cm"),
                           legend.text = element_text(size = 7), legend.title = element_text(size = 8)) +
    guides(fill = guide_legend(ncol = 1))
  plots <- c(plots, list(p)); heights <- c(heights, 0.22)
}
if (!is.null(RM) && nrow(RM)) {
  RM[, cls := sub("/.*$", "", label)]
  want <- csv(opt("--rmsk-classes"))
  if (!length(want)) want <- RM[, .(bp = sum(end - start)), by = cls][order(-bp)][1:min(6, .N), cls]
  miss <- setdiff(want, RM$cls); if (length(miss)) say("RepeatMasker classes not on ", chrom, ": ", paste(miss, collapse = ", "))
  rm <- RM[cls %in% want]
  if (nrow(rm)) {
    b <- copy(bins)[, `:=`(bs = start, be = end)]; setkey(b, bs, be)
    rm[, `:=`(s = as.numeric(start), e = as.numeric(end))]; setkey(rm, s, e)
    ov <- foverlaps(rm[, .(s, e, cls)], b, by.x = c("s", "e"), by.y = c("bs", "be"), nomatch = 0L)
    ov[, w := pmin(e, be) - pmax(s, bs)]
    den <- ov[, .(frac = pmin(1, sum(w) / (end[1] - start[1]))), by = .(cls, start)]
    den[, cls := factor(cls, rev(intersect(want, unique(cls))))]
    p <- ggplot(den) + geom_tile(aes(x = (start + bin / 2) / u, y = cls, fill = frac), width = bin / u, height = 0.9) +
      scale_fill_gradient(low = "#f0f0f0", high = "#3b3b3b", limits = c(0, 1), name = "RepeatMasker\nfraction of bin") +
      xs + labs(y = NULL) + base + blank_x +
      theme(axis.text.y = element_text(size = 7), legend.key.height = unit(0.3, "cm"), legend.title = element_text(size = 8),
            legend.text = element_text(size = 7), panel.grid = element_blank())
    plots <- c(plots, list(p)); heights <- c(heights, 0.12 * length(want) + 0.05)
  }
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
         if (smooth_k > 1) sprintf("smoothed over %d bins", smooth_k), if (!is.na(cap)) sprintf("y clipped at %.3g (triangles)", cap),
         if (any(!is.na(D$control)) && !has("--no-control-track")) "grey = control",
         if (nrow(peaks)) "orange ticks = peaks")
title <- opt("--title", sprintf("%s coverage", chrom))
fig <- wrap_plots(plots, ncol = 1, heights = heights, guides = "collect") +
  plot_annotation(title = title, subtitle = paste(c(ylab, sub), collapse = "; "),
                  theme = theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9)))
out <- out_for(norm, opt("--out", sprintf("coverage_%s.png", chrom)))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
w <- as.numeric(opt("--width", 12)); h <- as.numeric(opt("--height", 1.2 + 1.5 * length(levels(D$track)) + 0.6 * (last - length(levels(D$track)))))
ggsave(out, fig, width = w, height = h, dpi = 200, limitsize = FALSE, bg = "white")
say("wrote ", out)
}
