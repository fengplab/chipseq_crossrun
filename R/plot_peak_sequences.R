## Sequence grid of top peaks: the genomic sequence in a window centred on each peak's summit, one row per peak,
## each base a coloured cell with its letter (as in Figure E of the cGAS figures).
## Input: the top-peaks TSV written by R/plot_chrom_coverage.R --label-top (columns track, rank, chrom, summit).
## Usage:
##   Rscript R/plot_peak_sequences.R --top-peaks cov_top_peaks.tsv --crossrun <outdir> [options]
##   Rscript R/plot_peak_sequences.R --top-peaks cov_top_peaks.tsv --fasta genome.fa [--fasta-names names.tsv] [options]
## Options:
##   --crossrun DIR       take the genome FASTA (work/config.env) and chromosome -> FASTA name table (work/annot/fasta_names.tsv)
##   --fasta FILE         indexed genome FASTA (.fai next to it); overrides the crossrun one
##   --fasta-names FILE   two columns: canonical chromosome name, name in the FASTA (default: same names)
##   --track NAME         which track of the TSV (default: the first)
##   --window N           bases per row, centred on the summit (default 120, the width of the original Figure E grid)
##   --order rank|position   row order (default rank: 1 at the top)
##   --width PX --height PX  image size (default 3560 wide; height 30 + 50 per row, as in the original)
##   --dpi N              resolution written into the PNG (default 300)
##   --out FILE           PNG (default peak_sequences.png); the sequences are also written to <out>.tsv
## Colours sampled from the original figure: A #FBFF81 (yellow), C #00FD59 (green), G #00B0E6 (blue), T #FF93B1 (pink-red);
## N and anything else white.
a <- commandArgs(TRUE)
if (!length(a) || "--help" %in% a || "-h" %in% a) {
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
  h <- readLines(f); h <- h[seq_len(which(!grepl("^##", h))[1] - 1)]; cat(sub("^## ?", "", h), sep = "\n"); quit(status = 0)
}
known <- c("--top-peaks", "--crossrun", "--fasta", "--fasta-names", "--track", "--window", "--order", "--width", "--height", "--dpi", "--out")
bad <- setdiff(grep("^--", a, value = TRUE), known); if (length(bad)) { message("ERROR: unknown option(s): ", paste(bad, collapse = ", ")); quit(status = 1) }
opt <- function(n, d = NULL) { i <- match(n, a); if (is.na(i)) d else a[i + 1] }
die <- function(...) { message("ERROR: ", ...); quit(status = 1) }
say <- function(...) message("[peak sequences] ", ...)
suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
tp_f <- opt("--top-peaks") %||% die("--top-peaks is required (the TSV from plot_chrom_coverage.R --label-top)")
if (!file.exists(tp_f)) die("not found: ", tp_f)
TP <- fread(tp_f, colClasses = list(character = c("track", "chrom")))
need_cols <- c("track", "rank", "chrom", "summit"); if (!all(need_cols %in% names(TP))) die(tp_f, " lacks columns: ", paste(setdiff(need_cols, names(TP)), collapse = ", "))
trk <- opt("--track", TP$track[1]); if (!trk %in% TP$track) die("track '", trk, "' not in ", tp_f, " (", paste(unique(TP$track), collapse = ", "), ")")
TP <- TP[track == trk]
ord <- opt("--order", "rank"); if (!ord %in% c("rank", "position")) die("--order must be rank or position")
setorderv(TP, if (ord == "rank") "rank" else c("chrom", "summit"))

cr <- opt("--crossrun")
fa <- opt("--fasta")
if (is.null(fa) && !is.null(cr)) {
  ce <- file.path(cr, "work", "config.env"); if (!file.exists(ce)) die("not found: ", ce)
  l <- grep("^export GENOME_FASTA=", readLines(ce), value = TRUE); fa <- gsub("^export GENOME_FASTA=|[\"']", "", l[1])
}
if (is.null(fa) || !nzchar(fa)) die("give --fasta or --crossrun")
if (!file.exists(fa) || !file.exists(paste0(fa, ".fai"))) die("FASTA or its .fai index not found: ", fa)
fn_f <- opt("--fasta-names") %||% (if (!is.null(cr)) file.path(cr, "work", "annot", "fasta_names.tsv") else NULL)
FN <- if (!is.null(fn_f) && file.exists(fn_f)) fread(fn_f, header = FALSE, col.names = c("canon", "fasta"), colClasses = "character") else data.table(canon = character(), fasta = character())
fai <- fread(paste0(fa, ".fai"), header = FALSE, select = 1:2, col.names = c("name", "len"), colClasses = list(character = 1))
if (!nzchar(Sys.which("samtools"))) die("samtools not on PATH")

W <- as.integer(opt("--window", 120)); if (is.na(W) || W < 2) die("--window must be a whole number >= 2")
## Coordinates: the summit is a 0-based position (BED/narrowPeak convention). The window is [summit - floor(W/2),
## summit - floor(W/2) + W) in 0-based half-open terms; samtools faidx takes 1-based closed regions, so it is asked for
## (start0 + 1)-(start0 + W).
seqs <- rbindlist(lapply(seq_len(nrow(TP)), function(i) {
  ch <- TP$chrom[i]; fname <- FN$fasta[match(ch, FN$canon)]; if (is.na(fname)) fname <- ch
  if (!fname %in% fai$name) die(ch, " (", fname, ") is not in ", fa, ".fai; give --fasta-names")
  s0 <- TP$summit[i] - W %/% 2; e0 <- s0 + W
  if (s0 < 0 || e0 > fai$len[match(fname, fai$name)]) die("window of rank ", TP$rank[i], " runs off the end of ", ch)
  out <- system2("samtools", c("faidx", shQuote(fa), shQuote(sprintf("%s:%d-%d", fname, s0 + 1, e0))), stdout = TRUE)
  st <- attr(out, "status"); if (!is.null(st) && st != 0) die("samtools faidx failed for rank ", TP$rank[i])
  sq <- toupper(paste(out[-1], collapse = ""))
  if (nchar(sq) != W) die(sprintf("rank %d: got %d bases from samtools, expected %d", TP$rank[i], nchar(sq), W))
  data.table(rank = TP$rank[i], chrom = ch, start0 = s0, end0 = e0, summit = TP$summit[i], sequence = sq)
}))
out <- opt("--out", "peak_sequences.png")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
tsv <- paste0(sub("\\.png$", "", out), ".tsv")
fwrite(seqs, tsv, sep = "\t")
say(sprintf("%d sequences of %d bp around the summits (%s) -> %s", nrow(seqs), W, trk, tsv))

## ---------------- draw: pixel coordinates, origin top left
COL <- c(A = "#FBFF81", C = "#00FD59", G = "#00B0E6", T = "#FF93B1")
n <- nrow(seqs)
Wpx <- as.integer(opt("--width", 3560)); dpi <- as.numeric(opt("--dpi", 300))
rule_h <- 30; row_h <- 50
Hpx <- as.integer(opt("--height", rule_h + row_h * n)); row_h <- (Hpx - rule_h) / n
lab_w <- round(Wpx * 188 / 3621)                     # label column, in the original's proportion
gx0 <- lab_w; gx1 <- Wpx - 6; cw <- (gx1 - gx0) / W
if (cw < 4) say(sprintf("WARNING: cells are only %.1f px wide; letters are left out", cw))
pt_for <- function(cap_px) cap_px * 72 / (0.72 * dpi)  # font size whose capital letters are cap_px pixels tall
png(out, width = Wpx, height = Hpx, res = dpi, bg = "white")
op <- par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i"); plot.new(); plot.window(xlim = c(0, Wpx), ylim = c(Hpx, 0))
## ruler: a grey line over the grid with a tick every second base, as in the original
ry <- rule_h - 6
segments(gx0, ry, gx1, ry, col = "grey55", lwd = 4 * 96 / dpi, lend = 1)          # lwd 1 = 1/96 in: about 4 px
tk <- seq(0, W, by = 2); segments(gx0 + tk * cw, ry, gx0 + tk * cw, ry - 6, col = "grey55", lwd = 2 * 96 / dpi, lend = 1)
let_pt <- pt_for(min(0.5 * row_h, 0.85 * cw)); lab_pt <- pt_for(0.8 * row_h)
for (i in seq_len(n)) {
  y0 <- rule_h + (i - 1) * row_h; b <- strsplit(seqs$sequence[i], "")[[1]]
  fill <- COL[b]; fill[is.na(fill)] <- "white"
  rect(gx0 + (seq_len(W) - 1) * cw, y0, gx0 + seq_len(W) * cw, y0 + row_h, col = fill, border = NA)
  if (cw >= 4) text(gx0 + (seq_len(W) - 0.5) * cw, y0 + row_h / 2, b, cex = let_pt / 12, family = "sans", col = "black")
  text(lab_w - 0.12 * lab_w, y0 + row_h / 2, seqs$rank[i], adj = c(1, 0.5), font = 2, cex = lab_pt / 12, family = "sans")
}
par(op); invisible(dev.off())
say(sprintf("wrote %s (%d x %d px, %d dpi; %d rows x %d bases, cells %.1f x %.1f px)", out, Wpx, Hpx, dpi, n, W, cw, row_h))
