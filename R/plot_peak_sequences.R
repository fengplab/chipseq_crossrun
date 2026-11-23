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
##   --fasta-out FILE     also write the sequences as FASTA, one record per row in the order drawn: the genomic window
##                        on the + strand, ungapped (as extracted, before any --align). Header:
##                        >rank<R> <chrom>:<start>-<end> summit=<summit> track=<track>, coordinates 1-based closed
## Highlighting repeats (any combination):
##   --repeat-units       find each row's tandem repeat unit (the shortest period p <= --repeat-max-period, default 40,
##                        at which at least --repeat-min-identity, default 0.8, of bases equal the base p further on, with
##                        at least 3 copies in the row), box every copy (at most 20% mismatches to the unit) and print
##                        "UNIT xN" (units up to 8 bp) or "P bp xN" to the right of the row; the TSV has the full unit
##   --highlight M1,M2    outline every match of these motifs (IUPAC codes allowed, e.g. TTCGNNNNANNCGGG) on either strand.
##                        With --align, matches (and repeat copies) whose bases are split by alignment gaps are not shown;
##                        their number is logged and written to the TSV (highlight_gapped, repeat_copies_drawn)
##   --shared-kmer K      keep bases at full colour where a K-mer (either strand) occurs in at least --shared-min rows
##                        (default: half the rows, rounded up)
##   --align              align the rows with MAFFT (--auto --adjustdirection) so repeats fall into columns; reverse-
##                        complemented rows are marked * after the rank; gaps are white cells with "-"; bases that differ
##                        from their column's consensus are faded (--no-consensus-fade to keep them at full colour)
##   --fade F             strength of fading: 0 = none, 1 = white (default 0.7). With --highlight or --shared-kmer, every
##                        base that neither option marks is faded
## Colours sampled from the original figure: A #FBFF81 (yellow), C #00FD59 (green), G #00B0E6 (blue), T #FF93B1 (pink-red);
## N and anything else white.
a <- commandArgs(TRUE)
if (!length(a) || "--help" %in% a || "-h" %in% a) {
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
  h <- readLines(f); h <- h[seq_len(which(!grepl("^##", h))[1] - 1)]; cat(sub("^## ?", "", h), sep = "\n"); quit(status = 0)
}
known <- c("--top-peaks", "--crossrun", "--fasta", "--fasta-names", "--track", "--window", "--order", "--width", "--height", "--dpi", "--out", "--fasta-out",
           "--repeat-units", "--repeat-max-period", "--repeat-min-identity", "--highlight", "--shared-kmer", "--shared-min", "--align",
           "--no-consensus-fade", "--fade")
flags <- c("--repeat-units", "--align", "--no-consensus-fade")
bad <- setdiff(grep("^--", a, value = TRUE), known); if (length(bad)) { message("ERROR: unknown option(s): ", paste(bad, collapse = ", ")); quit(status = 1) }
opt <- function(n, d = NULL) { i <- match(n, a); if (is.na(i)) d else a[i + 1] }
has <- function(n) n %in% a
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
n <- nrow(seqs)
fade <- as.numeric(opt("--fade", 0.7)); if (is.na(fade) || fade < 0 || fade > 1) die("--fade must be between 0 and 1")
comp <- c(A = "T", C = "G", G = "C", T = "A", N = "N")
revcomp <- function(x) paste(rev(comp[strsplit(x, "")[[1]]]), collapse = "")

## ---------------- optional alignment (MAFFT); rows may come back reverse-complemented
seqs[, `:=`(oriented = sequence, strand = "+", aligned = sequence)]
if (has("--align")) {
  if (!nzchar(Sys.which("mafft"))) die("--align needs mafft on PATH (conda install -c bioconda mafft)")
  fin <- tempfile(fileext = ".fa"); writeLines(paste0(">r", seqs$rank, "\n", seqs$sequence), fin)
  al <- system2("mafft", c("--auto", "--adjustdirection", "--quiet", shQuote(fin)), stdout = TRUE)
  st <- attr(al, "status"); if (!is.null(st) && st != 0) die("mafft failed")
  hdr <- grep("^>", al); nm <- sub("^>", "", al[hdr])
  body <- vapply(seq_along(hdr), function(j) { e <- if (j < length(hdr)) hdr[j + 1] - 1 else length(al); toupper(paste(al[(hdr[j] + 1):e], collapse = "")) }, "")
  rc <- grepl("^_R_", nm); key <- sub("^_R_", "", nm)
  m <- match(paste0("r", seqs$rank), key); if (anyNA(m)) die("mafft output does not name every row")
  seqs[, `:=`(aligned = body[m], strand = ifelse(rc[m], "-", "+"))]
  seqs[, oriented := gsub("-", "", aligned)]
  bad_o <- seqs[nchar(oriented) != W | (strand == "+" & oriented != sequence) | (strand == "-" & oriented != vapply(sequence, revcomp, ""))]
  if (nrow(bad_o)) die("mafft changed the bases of rank(s) ", paste(bad_o$rank, collapse = ", "))
  say(sprintf("aligned %d rows with mafft: %d columns; reverse-complemented: %s", n, nchar(seqs$aligned[1]),
              if (any(seqs$strand == "-")) paste(seqs[strand == "-", rank], collapse = ", ") else "none"))
}
L <- nchar(seqs$aligned[1])
## column of each base of the oriented sequence
col_of <- lapply(seqs$aligned, function(x) which(strsplit(x, "")[[1]] != "-"))
## With --align, a match whose bases are split by alignment gaps in its row is not shown: a box drawn over it would
## cover gap cells and read as a match in the alignment. Spans are in oriented-sequence positions (start, end).
ungapped <- function(i, sp) col_of[[i]][sp[, 2]] - col_of[[i]][sp[, 1]] == sp[, 2] - sp[, 1]

## ---------------- 1. tandem repeat units
find_repeat <- function(x, maxp, minid) {
  b <- strsplit(x, "")[[1]]; w <- length(b)
  for (p in seq_len(min(maxp, w %/% 3))) {
    id <- mean(b[seq_len(w - p)] == b[(p + 1):w])
    if (id >= minid) {
      km <- vapply(seq_len(w - p + 1), function(i) paste(b[i:(i + p - 1)], collapse = ""), "")
      unit <- names(sort(table(km), decreasing = TRUE))[1]
      u <- strsplit(unit, "")[[1]]; maxmm <- floor(0.2 * p); boxes <- list(); i <- 1
      while (i <= w - p + 1) {
        if (sum(b[i:(i + p - 1)] != u) <= maxmm) { boxes[[length(boxes) + 1]] <- c(i, i + p - 1); i <- i + p } else i <- i + 1
      }
      ## label: the alphabetically first rotation of the unit, so the same repeat reads the same in every row
      rot <- vapply(0:(p - 1), function(k) paste(u[c(seq_len(p - k) + k, seq_len(k))], collapse = ""), "")
      if (length(boxes) >= 3) return(list(unit = min(rot), period = p, identity = id, boxes = do.call(rbind, boxes)))
    }
  }
  NULL
}
REP <- vector("list", n)
if (has("--repeat-units")) {
  maxp <- as.integer(opt("--repeat-max-period", 40)); minid <- as.numeric(opt("--repeat-min-identity", 0.8))
  REP <- lapply(seqs$oriented, find_repeat, maxp, minid)
  seqs[, `:=`(repeat_unit = vapply(REP, function(r) if (is.null(r)) NA_character_ else r$unit, ""),
              repeat_period = vapply(REP, function(r) if (is.null(r)) NA_integer_ else r$period, 1L),
              repeat_copies = vapply(REP, function(r) if (is.null(r)) NA_integer_ else nrow(r$boxes), 1L))]
  say("tandem repeat units: ", paste(sprintf("%s=%s", seqs$rank, ifelse(is.na(seqs$repeat_unit), "none", paste0(seqs$repeat_unit, "x", seqs$repeat_copies))), collapse = "  "))
  for (i in seq_len(n)) if (!is.null(REP[[i]])) { k <- ungapped(i, REP[[i]]$boxes); REP[[i]]$boxes <- REP[[i]]$boxes[k, , drop = FALSE] }
  seqs[, repeat_copies_drawn := vapply(REP, function(r) if (is.null(r)) NA_integer_ else nrow(r$boxes), 1L)]
  if (has("--align") && any(seqs$repeat_copies_drawn < seqs$repeat_copies, na.rm = TRUE))
    say(sprintf("repeat copies split by alignment gaps, not boxed: %d (the label counts the boxed copies; the TSV has both)",
                sum(seqs$repeat_copies - seqs$repeat_copies_drawn, na.rm = TRUE)))
}

## ---------------- 2. motif matches (IUPAC, both strands)
IUPAC <- c(A = "A", C = "C", G = "G", T = "T", R = "[AG]", Y = "[CT]", S = "[CG]", W = "[AT]", K = "[GT]", M = "[AC]",
           B = "[CGT]", D = "[AGT]", H = "[ACT]", V = "[ACG]", N = "[ACGT]")
HL <- vector("list", n); hl_mask <- lapply(seqs$oriented, function(x) rep(FALSE, nchar(x)))
if (!is.null(opt("--highlight"))) {
  mot <- toupper(trimws(strsplit(opt("--highlight"), ",")[[1]])); mot <- mot[nzchar(mot)]
  badm <- mot[grepl("[^ACGTRYSWKMBDHVN]", mot)]; if (length(badm)) die("--highlight: not IUPAC DNA: ", paste(badm, collapse = ", "))
  rc_iupac <- function(m) { cm <- c(A = "T", C = "G", G = "C", T = "A", R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
                                     B = "V", D = "H", H = "D", V = "B", N = "N"); paste(rev(cm[strsplit(m, "")[[1]]]), collapse = "") }
  pats <- unique(c(mot, vapply(mot, rc_iupac, "")))
  seqs[, highlight_gapped := 0L]
  for (i in seq_len(n)) {
    x <- seqs$oriented[i]; spans <- NULL
    for (pt in pats) {
      rx <- paste0("(?=", paste(IUPAC[strsplit(pt, "")[[1]]], collapse = ""), ")")
      st <- gregexpr(rx, x, perl = TRUE)[[1]]
      if (st[1] > 0) spans <- rbind(spans, cbind(as.integer(st), as.integer(st) + nchar(pt) - 1))
    }
    gapped <- 0L
    if (!is.null(spans)) {
      spans <- unique(spans); k <- ungapped(i, spans); gapped <- sum(!k); spans <- spans[k, , drop = FALSE]
      if (nrow(spans)) { HL[[i]] <- spans; for (r in seq_len(nrow(spans))) hl_mask[[i]][spans[r, 1]:spans[r, 2]] <- TRUE }
    }
    seqs[i, highlight_gapped := gapped]
  }
  seqs[, highlight_matches := vapply(HL, function(h) if (is.null(h)) 0L else nrow(h), 1L)]
  say(sprintf("--highlight %s: %d matches in %d of %d rows%s", paste(mot, collapse = ","), sum(seqs$highlight_matches),
              sum(seqs$highlight_matches > 0), n,
              if (sum(seqs$highlight_gapped)) sprintf("; %d more split by alignment gaps, not shown (column highlight_gapped)", sum(seqs$highlight_gapped)) else ""))
}

## ---------------- 3. shared k-mers
sh_mask <- lapply(seqs$oriented, function(x) rep(FALSE, nchar(x)))
if (!is.null(opt("--shared-kmer"))) {
  K <- as.integer(opt("--shared-kmer")); if (is.na(K) || K < 3 || K > W) die("--shared-kmer must be between 3 and the window")
  smin <- as.integer(opt("--shared-min", ceiling(n / 2)))
  canon_k <- function(k) { r <- revcomp(k); ifelse(k < r, k, r) }
  kms <- lapply(seqs$oriented, function(x) vapply(seq_len(nchar(x) - K + 1), function(i) canon_k(substr(x, i, i + K - 1)), ""))
  cnt <- table(unlist(lapply(kms, unique)))
  common <- names(cnt)[cnt >= smin]
  for (i in seq_len(n)) for (j in which(kms[[i]] %in% common)) sh_mask[[i]][j:(j + K - 1)] <- TRUE
  seqs[, shared_bases := vapply(sh_mask, sum, 1L)]
  say(sprintf("--shared-kmer %d: %d %d-mers occur in >= %d rows; %.0f%% of bases kept at full colour", K, length(common), K, smin,
              100 * sum(seqs$shared_bases) / (n * W)))
}

## ---------------- fading: per row and column
emph_opts <- !is.null(opt("--highlight")) || !is.null(opt("--shared-kmer"))
cons_fade <- has("--align") && !has("--no-consensus-fade")
G <- do.call(rbind, strsplit(seqs$aligned, ""))                        # rows x columns, "-" for gaps
cons <- apply(G, 2, function(v) { v <- v[v != "-"]; if (!length(v)) "-" else names(which.max(table(v))) })
faded <- matrix(FALSE, n, L)
for (i in seqs[, seq_len(.N)]) {
  cc <- col_of[[i]]
  if (emph_opts) faded[i, cc] <- !(hl_mask[[i]] | sh_mask[[i]])
  if (cons_fade) faded[i, cc] <- faded[i, cc] | (G[i, cc] != cons[cc])
}
tsv <- paste0(sub("\\.png$", "", out), ".tsv")
fwrite(seqs, tsv, sep = "\t")
say(sprintf("%d sequences of %d bp around the summits (%s) -> %s", n, W, trk, tsv))
fa_out <- opt("--fasta-out")
if (!is.null(fa_out)) {
  if (is.na(fa_out) || startsWith(fa_out, "--")) die("--fasta-out needs a file name")
  dir.create(dirname(fa_out), showWarnings = FALSE, recursive = TRUE)
  ## the window is 0-based half-open [start0, end0); the header gives it 1-based closed: (start0 + 1)-end0, the same
  ## region samtools faidx was asked for. The record is the + strand genomic sequence, as in the 'sequence' column.
  rec <- seqs[, paste0(">rank", rank, " ", chrom, ":", start0 + 1, "-", end0, " summit=", summit, " track=", trk, "\n",
                       vapply(sequence, function(z) paste(substring(z, seq(1, nchar(z), 60), pmin(seq(60, nchar(z) + 59, 60), nchar(z))), collapse = "\n"), ""))]
  tmp <- paste0(fa_out, ".tmp"); writeLines(rec, tmp)
  if (!file.rename(tmp, fa_out)) die("could not write ", fa_out)
  say(sprintf("%d sequences -> %s (FASTA, + strand, ungapped)", n, fa_out))
}

## ---------------- draw: pixel coordinates, origin top left
COL <- c(A = "#FBFF81", C = "#00FD59", G = "#00B0E6", T = "#FF93B1")
mix_white <- function(cl, f) { m <- grDevices::col2rgb(cl) / 255; grDevices::rgb(t(m * (1 - f) + f)) }
Wpx <- as.integer(opt("--width", 3560)); dpi <- as.numeric(opt("--dpi", 300))
rule_h <- 30; row_h <- 50
Hpx <- as.integer(opt("--height", rule_h + row_h * n)); row_h <- (Hpx - rule_h) / n
lab_w <- round(Wpx * 188 / 3621)                     # label column, in the original's proportion
rep_w <- if (has("--repeat-units")) round(Wpx * 0.085) else 0   # right-hand column for "UNIT xN"
gx0 <- lab_w; gx1 <- Wpx - 6 - rep_w; cw <- (gx1 - gx0) / L
if (cw < 4) say(sprintf("WARNING: cells are only %.1f px wide; letters are left out", cw))
pt_for <- function(cap_px) cap_px * 72 / (0.72 * dpi)  # font size whose capital letters are cap_px pixels tall
px_lwd <- function(px) px * 96 / dpi                   # line width in px -> lwd (1 lwd = 1/96 in)
png(out, width = Wpx, height = Hpx, res = dpi, bg = "white")
op <- par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i"); plot.new(); plot.window(xlim = c(0, Wpx), ylim = c(Hpx, 0))
## ruler: a grey line over the grid with a tick every second column, as in the original
ry <- rule_h - 6
segments(gx0, ry, gx1, ry, col = "grey55", lwd = px_lwd(4), lend = 1)
tk <- seq(0, L, by = 2); segments(gx0 + tk * cw, ry, gx0 + tk * cw, ry - 6, col = "grey55", lwd = px_lwd(2), lend = 1)
let_pt <- pt_for(min(0.5 * row_h, 0.85 * cw)); lab_pt <- pt_for(0.8 * row_h)
for (i in seq_len(n)) {
  y0 <- rule_h + (i - 1) * row_h; b <- G[i, ]
  fill <- COL[b]; fill[is.na(fill)] <- "white"
  fill[faded[i, ] & b != "-"] <- mix_white(fill[faded[i, ] & b != "-"], fade)
  xl <- gx0 + (seq_len(L) - 1) * cw
  rect(xl, y0, xl + cw, y0 + row_h, col = fill, border = NA)
  if (cw >= 4) text(xl + cw / 2, y0 + row_h / 2, b, cex = let_pt / 12, family = "sans",
                    col = ifelse(faded[i, ], grDevices::grey(0.15 + 0.6 * fade), "black"))
  ## spans (in oriented-sequence positions) -> column rectangles
  box <- function(sp, border, lw, inset) for (r in seq_len(nrow(sp))) {
    c1 <- col_of[[i]][sp[r, 1]]; c2 <- col_of[[i]][sp[r, 2]]
    rect(gx0 + (c1 - 1) * cw + inset, y0 + inset, gx0 + c2 * cw - inset, y0 + row_h - inset, border = border, lwd = px_lwd(lw), col = NA)
  }
  if (!is.null(REP[[i]])) box(REP[[i]]$boxes, "grey15", 2, 2)
  if (!is.null(HL[[i]])) box(HL[[i]], "black", 5, 0)
  lbl <- paste0(seqs$rank[i], if (seqs$strand[i] == "-") "*" else "")
  text(lab_w - 0.12 * lab_w, y0 + row_h / 2, lbl, adj = c(1, 0.5), font = 2, cex = lab_pt / 12, family = "sans")
  if (rep_w > 0 && !is.null(REP[[i]]))
    text(gx1 + 12, y0 + row_h / 2,
         if (REP[[i]]$period <= 8) sprintf("%s x%d", REP[[i]]$unit, nrow(REP[[i]]$boxes)) else sprintf("%d bp x%d", REP[[i]]$period, nrow(REP[[i]]$boxes)),
         adj = c(0, 0.5), cex = pt_for(min(0.42 * row_h, rep_w / 12 * 1.1)) / 12, family = "sans")
}
par(op); invisible(dev.off())
say(sprintf("wrote %s (%d x %d px, %d dpi; %d rows x %d columns, cells %.1f x %.1f px)%s", out, Wpx, Hpx, dpi, n, L, cw, row_h,
            if (any(seqs$strand == "-")) "; * = reverse complement" else ""))
