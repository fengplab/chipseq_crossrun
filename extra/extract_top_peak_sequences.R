## Top peaks of one sample across the whole genome, as FASTA: a fixed window centred on each peak's summit.
## Peaks are ranked exactly as R/plot_chrom_coverage.R --label-top ranks them (same column, same tie rule, same
## summit), only over every contig instead of one chromosome.
## Usage:
##   Rscript R/extract_top_peak_sequences.R --crossrun <outdir> --sample <sample_id> [options]
##   Rscript R/extract_top_peak_sequences.R --peaks sample_peaks.narrowPeak --fasta genome.fa [--fasta-names names.tsv]
##                                          [--chrom-alias alias.tsv] [--genome-sizes genome.sizes] [options]
## Options:
##   --crossrun DIR       take the sample's MACS peak file (work/samples.tsv peaks_file), the genome FASTA
##                        (work/config.env), FASTA names (work/annot/fasta_names.tsv), chromosome aliases
##                        (work/annot/chrom_alias.tsv) and the host contigs (work/annot/genome.sizes)
##   --sample ID          the sample (with --crossrun)
##   --peaks FILE         narrowPeak/broadPeak file (.gz allowed); overrides the crossrun one
##   --fasta FILE         indexed genome FASTA (.fai next to it); overrides the crossrun one
##   --fasta-names FILE   canonical chromosome name <TAB> name in the FASTA (default: same names)
##   --chrom-alias FILE   name in the peak file <TAB> canonical name (as in plot_chrom_coverage.R)
##   --genome-sizes FILE  contigs peaks may come from (default with --crossrun: work/annot/genome.sizes, i.e. the host
##                        nuclear chromosomes crossrun analyses -- no chrM, no viral genomes)
##   --all-contigs        rank peaks on every contig of the peak file, ignoring --genome-sizes
##   --top N              number of peaks (default 10)
##   --top-by score|signal|p|q   ranking column: 5 score (default, as plot_chrom_coverage.R), 7 signal value,
##                        8 -log10 p, 9 -log10 q. Ties: higher signal value first, then lower position
##   --window N           bases per sequence, centred on the summit (default 120, as plot_peak_sequences.R)
##   --out FILE           FASTA (default <sample>_top<N>_w<window>.fa); the ranked peaks with their windows and
##                        sequences go to <out without .fa>.tsv, which R/plot_peak_sequences.R accepts as --top-peaks
## Summit: narrowPeak column 10 (0-based offset from the start), or the midpoint when there is none -- as in
## plot_chrom_coverage.R. MACS can list one peak several times, once per summit (--call-summits); like
## plot_chrom_coverage.R these are ranked as separate peaks, and overlapping windows among the top N are reported.
## Sequences are the + strand genomic window, ungapped. A window running off the end of its contig stops the script
## (give a smaller --window): the window size is fixed, never clipped.
a <- commandArgs(TRUE)
if (!length(a) || "--help" %in% a || "-h" %in% a) {
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
  h <- readLines(f); h <- h[seq_len(which(!grepl("^##", h))[1] - 1)]; cat(sub("^## ?", "", h), sep = "\n"); quit(status = 0)
}
known <- c("--crossrun", "--sample", "--peaks", "--fasta", "--fasta-names", "--chrom-alias", "--genome-sizes", "--all-contigs",
           "--top", "--top-by", "--window", "--out")
bad <- setdiff(grep("^--", a, value = TRUE), known); if (length(bad)) { message("ERROR: unknown option(s): ", paste(bad, collapse = ", ")); quit(status = 1) }
opt <- function(n, d = NULL) { i <- match(n, a); if (is.na(i)) d else if (i == length(a) || startsWith(a[i + 1], "--")) die(n, " needs a value") else a[i + 1] }
has <- function(n) n %in% a
die <- function(...) { message("ERROR: ", ...); quit(status = 1) }
say <- function(...) message("[top peak sequences] ", ...)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
suppressPackageStartupMessages(library(data.table))
options(datatable.integer64 = "double")

cr <- opt("--crossrun"); smp <- opt("--sample")
if (!is.null(cr) && is.null(smp)) die("--crossrun needs --sample")
pk_f <- opt("--peaks"); fa <- opt("--fasta")
if (!is.null(cr)) {
  sm <- fread(file.path(cr, "work", "samples.tsv"), colClasses = "character")
  if (!smp %in% sm$sample_id) die(smp, " is not in ", cr, "/work/samples.tsv")
  if (is.null(pk_f)) { pk_f <- sm[sample_id == smp, peaks_file %||% NA_character_][1]
    if (is.na(pk_f) || !nzchar(pk_f)) die(smp, " has no peaks_file in samples.tsv; give --peaks") }
  if (is.null(fa)) { ce <- file.path(cr, "work", "config.env"); if (!file.exists(ce)) die("not found: ", ce)
    l <- grep("^export GENOME_FASTA=", readLines(ce), value = TRUE); fa <- gsub("^export GENOME_FASTA=|[\"']", "", l[1]) }
}
if (is.null(pk_f)) die("give --crossrun with --sample, or --peaks")
if (!file.exists(pk_f)) die("peak file not found: ", pk_f)
if (is.null(fa) || !nzchar(fa)) die("give --fasta or --crossrun")
if (!file.exists(fa) || !file.exists(paste0(fa, ".fai"))) die("FASTA or its .fai index not found: ", fa, " (samtools faidx makes the index)")
if (!nzchar(Sys.which("samtools"))) die("samtools not on PATH")
smp <- smp %||% sub("_peaks\\.(narrow|broad)Peak(\\.gz)?$", "", basename(pk_f))

rdtab <- function(f, nm) if (!is.null(f) && file.exists(f)) fread(f, header = FALSE, select = 1:2, col.names = nm, colClasses = "character") else NULL
alias <- rdtab(opt("--chrom-alias") %||% (if (!is.null(cr)) file.path(cr, "work", "annot", "chrom_alias.tsv")), c("from", "to"))
canon <- function(x) { if (is.null(alias)) return(x); i <- match(x, alias$from); ifelse(is.na(i), x, alias$to[i]) }
FN <- rdtab(opt("--fasta-names") %||% (if (!is.null(cr)) file.path(cr, "work", "annot", "fasta_names.tsv")), c("canon", "fasta"))
gs_f <- if (has("--all-contigs")) NULL else opt("--genome-sizes") %||% (if (!is.null(cr)) file.path(cr, "work", "annot", "genome.sizes"))
GS <- if (!is.null(gs_f)) { if (!file.exists(gs_f)) die("not found: ", gs_f); fread(gs_f, header = FALSE, select = 1, col.names = "chrom", colClasses = "character") } else NULL
fai <- fread(paste0(fa, ".fai"), header = FALSE, select = 1:2, col.names = c("name", "len"), colClasses = list(character = 1))

n_top <- suppressWarnings(as.integer(opt("--top", 10))); if (is.na(n_top) || n_top < 1) die("--top must be a whole number >= 1")
W <- suppressWarnings(as.integer(opt("--window", 120))); if (is.na(W) || W < 2) die("--window must be a whole number >= 2")
top_by <- opt("--top-by", "score"); top_col <- c(score = 5L, signal = 7L, p = 8L, q = 9L)[top_by]
if (is.na(top_col)) die("--top-by must be score, signal, p or q")

## ---------------- peaks, read and ranked as in plot_chrom_coverage.R --label-top
## Peak files are 0-based half-open; the narrowPeak summit (column 10) is a 0-based offset from the start, so
## summit = start + offset is a 0-based position. Without a summit, the midpoint is used.
x <- fread(cmd = if (grepl("\\.gz$", pk_f)) sprintf("gzip -dc '%s'", pk_f) else sprintf("cat '%s'", pk_f), header = FALSE, fill = TRUE)
x <- x[!grepl("^(track|browser|#)", V1)]
if (!nrow(x)) die(pk_f, " has no peaks")
if (ncol(x) < top_col) die(pk_f, " has no column ", top_col, " (", top_by, "); use --top-by with a column it has")
nm <- function(k) if (ncol(x) >= k) suppressWarnings(as.numeric(x[[k]])) else rep(NA_real_, nrow(x))
d <- data.table(track = smp, chrom = canon(as.character(x$V1)), start = as.numeric(x$V2), end = as.numeric(x$V3),
                name = if (ncol(x) >= 4) as.character(x$V4) else NA_character_,
                score = nm(5), signal = nm(7), neglog10_p = nm(8), neglog10_q = nm(9), off = nm(10))
d[, summit := fifelse(!is.na(off) & off >= 0, start + off, floor((start + end) / 2))][, off := NULL]
n_in <- nrow(d)
if (!is.null(GS)) {
  out_c <- d[!chrom %in% GS$chrom, .N, by = chrom]
  d <- d[chrom %in% GS$chrom]
  if (nrow(out_c)) say(sprintf("%d of %d peaks on contigs outside %s left out (%s); --all-contigs keeps them", sum(out_c$N), n_in,
                               basename(gs_f), paste(sprintf("%s %d", out_c$chrom, out_c$N), collapse = ", ")))
}
d[, key := get(c(score = "score", signal = "signal", p = "neglog10_p", q = "neglog10_q")[[top_by]])]
n_na <- d[is.na(key), .N]; d <- d[!is.na(key)]
if (n_na) say(n_na, " peaks have no ", top_by, " value and cannot be ranked")
if (!nrow(d)) die("no peaks left to rank")
setorder(d, -key, -signal, chrom, start)          # ties: higher signal first, then position
if (nrow(d) < n_top) say(sprintf("only %d peaks to rank; writing all of them", nrow(d)))
TP <- head(d, n_top)[, rank := seq_len(.N)][, key := NULL]
say(sprintf("%s: %d peaks ranked by %s across %d contigs; top %d kept", smp, nrow(d), top_by, uniqueN(d$chrom), nrow(TP)))

## ---------------- fixed windows
## The window is [summit - floor(W/2), summit - floor(W/2) + W) in 0-based half-open terms; samtools faidx and the
## FASTA header use 1-based closed coordinates, (start0 + 1)-end0.
TP[, `:=`(win_start0 = summit - W %/% 2, win_end0 = summit - W %/% 2 + W)]
seqs <- vapply(seq_len(nrow(TP)), function(i) {
  ch <- TP$chrom[i]; fname <- if (!is.null(FN)) FN$fasta[match(ch, FN$canon)] else NA_character_; if (is.na(fname)) fname <- ch
  if (!fname %in% fai$name) die(ch, " (", fname, ") is not in ", fa, ".fai; give --fasta-names")
  s0 <- TP$win_start0[i]; e0 <- TP$win_end0[i]
  if (s0 < 0 || e0 > fai$len[match(fname, fai$name)])
    die(sprintf("rank %d (%s:%d) is within %d bp of its contig's end, so a %d bp window runs off it; use a smaller --window",
                TP$rank[i], ch, TP$summit[i], W %/% 2, W))
  out <- system2("samtools", c("faidx", shQuote(fa), shQuote(sprintf("%s:%d-%d", fname, s0 + 1, e0))), stdout = TRUE)
  st <- attr(out, "status"); if (!is.null(st) && st != 0) die("samtools faidx failed for rank ", TP$rank[i])
  sq <- toupper(paste(out[-1], collapse = ""))
  if (nchar(sq) != W) die(sprintf("rank %d: got %d bases from samtools, expected %d", TP$rank[i], nchar(sq), W))
  sq
}, "")
TP[, sequence := seqs]
ovl <- character()
for (i in seq_len(nrow(TP))) for (j in seq_len(nrow(TP))) if (i < j && TP$chrom[i] == TP$chrom[j] &&
    TP$win_start0[i] < TP$win_end0[j] && TP$win_start0[j] < TP$win_end0[i]) ovl <- c(ovl, paste0(TP$rank[i], "&", TP$rank[j]))
if (length(ovl)) say("windows that overlap (often two summits of one MACS peak): ranks ", paste(ovl, collapse = ", "))

out <- opt("--out", sprintf("%s_top%d_w%d.fa", smp, n_top, W))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
wrap60 <- function(z) paste(substring(z, seq(1, nchar(z), 60), pmin(seq(60, nchar(z) + 59, 60), nchar(z))), collapse = "\n")
rec <- TP[, paste0(">rank", rank, " ", chrom, ":", win_start0 + 1, "-", win_end0, " summit=", summit, " ", top_by, "=",
                   get(c(score = "score", signal = "signal", p = "neglog10_p", q = "neglog10_q")[[top_by]]), " sample=", track, "\n",
                   vapply(sequence, wrap60, ""))]
writeLines(rec, paste0(out, ".tmp")); if (!file.rename(paste0(out, ".tmp"), out)) die("could not write ", out)
tsv <- paste0(sub("\\.(fa|fasta|fna)$", "", out), ".tsv")
fwrite(TP[, .(track, rank, chrom, start, end, summit, name, score, signal, neglog10_p, neglog10_q, ranked_by = top_by,
              window_start0 = win_start0, window_end0 = win_end0, sequence)], tsv, sep = "\t", na = "NA", quote = FALSE, scipen = 100)
say(sprintf("%d sequences of %d bp -> %s (FASTA, + strand); table -> %s (peak start/end and window_start0/end0 0-based half-open)",
            nrow(TP), W, out, tsv))
