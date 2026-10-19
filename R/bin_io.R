## Readers for the bin table written by bin/bin_table.py, and the signal measures shared by the bin-size tuning
## (R/tune_bins.R) and the drift check (R/16_drift_check.R). Source after utils.R or on its own (needs data.table).
suppressPackageStartupMessages(library(data.table))

## one column of one chromosome: <dir>/<chrom>/<col>.<type>, raw little-endian
read_col <- function(dir, chrom, col, n) {
  f <- Sys.glob(file.path(dir, chrom, paste0(col, ".*")))
  f <- f[!grepl("/\\.", f)]
  if (!length(f)) return(NULL)
  type <- sub("^.*\\.", "", f[1])
  switch(type,
         float32 = readBin(f[1], "double", n, size = 4, endian = "little"),
         int16 = readBin(f[1], "integer", n, size = 2, signed = TRUE, endian = "little"),
         int32 = readBin(f[1], "integer", n, size = 4, endian = "little"),
         stop("unknown column type: ", f[1]))
}

bin_manifest <- function(dir) {
  m <- fread(file.path(dir, "manifest.tsv"), colClasses = list(character = "chrom"))
  m[, bin_size := as.integer(bin_size)][]
}

## annotation columns for some chromosomes (all if chroms = NULL)
load_annot <- function(dir, chroms = NULL, cols = c("class", "array_id", "edge_dist", "repeat_frac", "mappability", "gc", "n_frac")) {
  m <- bin_manifest(dir); if (!is.null(chroms)) m <- m[chrom %in% chroms]
  cl <- fread(file.path(dir, "classes.tsv"))
  rbindlist(lapply(seq_len(nrow(m)), function(i) {
    n <- m$n_bins[i]; d <- data.table(chrom = m$chrom[i], bin = seq_len(n) - 1L,
                                      bin_len = c(rep(m$bin_size[i], n - 1), m$length[i] - (n - 1) * m$bin_size[i]))
    for (cc in cols) { v <- read_col(dir, m$chrom[i], cc, n); if (!is.null(v)) set(d, j = cc, value = v) }
    d
  }))[, class_name := cl$class[match(class, cl$code)]][, class_source := cl$source[match(class, cl$code)]][]
}

## signal columns of one library (dir = <bins>/w<size>/<library name>)
load_signal <- function(dir, chroms = NULL, cols = c("all", "uniq", "h1", "h2", "fl")) {
  m <- bin_manifest(dir); if (!is.null(chroms)) m <- m[chrom %in% chroms]
  out <- rbindlist(lapply(seq_len(nrow(m)), function(i) {
    n <- m$n_bins[i]; d <- data.table(chrom = m$chrom[i], bin = seq_len(n) - 1L)
    for (cc in cols) { v <- read_col(dir, m$chrom[i], cc, n); if (!is.null(v)) set(d, j = cc, value = v) }
    d
  }))
  lib <- jsonlite_read(file.path(dir, "library.json"))
  attr(out, "library") <- lib
  out
}
jsonlite_read <- function(f) {                       # library.json is flat: avoid a jsonlite dependency
  x <- paste(readLines(f, warn = FALSE), collapse = " ")
  kv <- regmatches(x, gregexpr('"[^"]+"\\s*:\\s*("[^"]*"|\\[[^]]*\\]|[-0-9.eE]+|true|false)', x))[[1]]
  setNames(lapply(kv, function(p) { v <- trimws(sub('^"[^"]+"\\s*:\\s*', "", p))
    if (grepl('^"', v)) gsub('"', "", v) else if (v %in% c("true", "false")) v == "true" else if (grepl("^\\[", v)) v else as.numeric(v) }),
    sub('^"([^"]+)".*$', "\\1", kv))
}

## exact aggregation of per-bin means from a fine size to k times that size (bins are aligned at 0)
aggregate_bins <- function(d, k, value_cols, len_col = "bin_len") {
  d <- copy(d)[, g := bin %/% k]
  d[, c(lapply(.SD, function(v) sum(v * get(len_col)) / sum(get(len_col))), list(bin_len = sum(get(len_col)))),
    by = .(chrom, bin = g), .SDcols = value_cols]
}

## binding: log2 IP/input with both scaled to reads per million fragments; eps = 1st percentile of non-zero input
binding <- function(ip, input, n_ip, n_input) {
  a <- ip / (n_ip / 1e6); b <- input / (n_input / 1e6)
  eps <- as.numeric(quantile(b[b > 0], 0.01, na.rm = TRUE)); if (!is.finite(eps) || eps <= 0) eps <- 1e-3
  log2((a + eps) / (b + eps))
}

## signal retained (split-half covariance) and reproducibility (split-half correlation) over a set of bins
split_half <- function(y1, y2) {
  ok <- is.finite(y1) & is.finite(y2)
  if (sum(ok) < 50) return(list(cov = NA_real_, cor = NA_real_, n = sum(ok)))
  list(cov = cov(y1[ok], y2[ok]), cor = cor(y1[ok], y2[ok]), n = sum(ok))
}

## Split-half signal at several bin sizes for one IP/control pair, from tables at a finer base size.
## annot: load_annot() at the base size; ip, ctl: load_signal() at the base size (columns h1, h2).
## Returns one row per size and group (satellite = CenSat bins, other): split-half covariance and correlation.
retention_curves <- function(annot, ip, ctl, sizes, base) {
  if (any(sizes %% base != 0)) stop("every size must be a multiple of the base size ", base)
  A <- annot[, .(chrom, bin, bin_len, sat = as.numeric(!is.na(class_source) & class_source == "censat"),
                 nf = if ("n_frac" %in% names(annot)) n_frac else 0)]
  same_bins(A, ip); same_bins(A, ctl)                  # same chromosomes and bins, same order: bind columns directly
  d <- cbind(A, ip[, .(ip1 = h1, ip2 = h2)], ctl[, .(in1 = h1, in2 = h2)])
  Li <- attr(ip, "library"); Lc <- attr(ctl, "library")
  rbindlist(lapply(sizes, function(w) {
    a <- if (w == base) copy(d) else aggregate_bins(d, w %/% base, c("sat", "nf", "ip1", "ip2", "in1", "in2"))
    keep <- a$nf < 0.5 & (a$in1 + a$in2) > 0
    y1 <- binding(a$ip1, a$in1, Li$h1, Lc$h1); y2 <- binding(a$ip2, a$in2, Li$h2, Lc$h2)
    rbindlist(lapply(c("satellite", "other"), function(g) {
      sel <- keep & ((a$sat >= 0.5) == (g == "satellite"))
      sh <- split_half(y1[sel], y2[sel])
      data.table(size = w, group = g, cov = sh$cov, cor = sh$cor, n_bins = sh$n)
    }))
  }))
}

same_bins <- function(a, b) {
  if (nrow(a) != nrow(b) || !identical(a$chrom, b$chrom) || !identical(a$bin, b$bin))
    stop("bin tables do not line up (different chromosomes or bin size): rebuild them with the same --bin and --chroms")
  invisible(TRUE)
}

## IP/control pairs from a libraries.tsv (sample_id, library, name, dir)
lib_pairs <- function(libs) {
  ip <- libs[library == "IP", .(sample_id, ip_dir = dir)]
  ct <- libs[library == "control", .(sample_id, ctl_dir = dir)]
  merge(ip, ct, by = "sample_id")
}
