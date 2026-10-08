## Shared helpers for the R analysis scripts.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

## Counts above 2^31 (e.g. genome-wide base counts) would otherwise be read as integer64, which needs the
## bit64 package; without it the values are silently garbled in arithmetic. Doubles are exact up to 2^53.
options(datatable.integer64 = "double")

args_config <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) < 1) stop("usage: Rscript <script>.R config.yaml")
  cfg <- yaml::read_yaml(a[1])
  cfg$outdir <- normalizePath(cfg$outdir %||% "crossrun_out", mustWork = FALSE)
  cfg$work <- file.path(cfg$outdir, "work")
  cfg$tables <- file.path(cfg$outdir, "tables")
  cfg$figures <- file.path(cfg$outdir, "figures")
  dir.create(cfg$tables, showWarnings = FALSE, recursive = TRUE)
  dir.create(cfg$figures, showWarnings = FALSE, recursive = TRUE)
  cfg
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

msg <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)

rd <- function(path, ...) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  fread(path, sep = "\t", quote = "", na.strings = c("", "NA"), ...)
}

wt <- function(dt, cfg, name) {
  f <- file.path(cfg$tables, paste0(name, ".tsv"))
  fwrite(dt, f, sep = "\t", na = "NA", quote = FALSE)
  msg("wrote ", f)
  invisible(f)
}

save_plot <- function(p, cfg, name, width = 7, height = 5) {
  for (ext in c("pdf", "png")) {
    f <- file.path(cfg$figures, paste0(name, ".", ext))
    if (ext == "png") ggsave(f, p, width = width, height = height, dpi = 150, bg = "white")
    else ggsave(f, p, width = width, height = height, device = cairo_pdf_or_pdf())
  }
  msg("figure ", name)
}

cairo_pdf_or_pdf <- function() if (capabilities("cairo")) cairo_pdf else "pdf"

## pheatmap to pdf + png
save_pheatmap <- function(fun, cfg, name, width = 8, height = 6) {
  grDevices::pdf(file.path(cfg$figures, paste0(name, ".pdf")), width = width, height = height)
  fun(); grDevices::dev.off()
  grDevices::png(file.path(cfg$figures, paste0(name, ".png")), width = width, height = height,
                 units = "in", res = 150)
  fun(); grDevices::dev.off()
  msg("figure ", name)
}

theme_x <- function(base = 10) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

samples_dt <- function(cfg) {
  s <- rd(file.path(cfg$work, "samples.tsv"), colClasses = "character")
  if (is.null(s)) stop("run 00_setup.py first")
  for (v in unlist(cfg$numeric_covariates)) if (v %in% names(s)) s[[v]] <- suppressWarnings(as.numeric(s[[v]]))
  if (!is.null(cfg$time_var) && cfg$time_var %in% names(s)) s[[cfg$time_var]] <- suppressWarnings(as.numeric(s[[cfg$time_var]]))
  s
}

## annotation data.frame for pheatmap (rows = sample ids)
sample_annotation <- function(cfg, s, ids) {
  vars <- unique(c(unlist(cfg$categorical_covariates), unlist(cfg$numeric_covariates)))  # run/batch ID deliberately excluded
  vars <- vars[vars %in% names(s)]
  a <- as.data.frame(s[match(ids, s$sample_id), ..vars])
  rownames(a) <- ids
  for (v in vars) if (!is.numeric(a[[v]])) a[[v]] <- factor(ifelse(is.na(a[[v]]) | a[[v]] == "", "NA", a[[v]]))
  a
}

bh <- function(p) p.adjust(p, method = "BH")

## short sample label for plots
short_id <- function(x) sub("_REP[0-9]+$", "", x)

## read a pooled tomtom text output
read_tomtom <- function(path) {
  if (!file.exists(path)) return(NULL)
  l <- readLines(path)
  l <- l[!grepl("^#", l) & nzchar(l)]
  if (length(l) < 2) return(NULL)
  d <- fread(text = paste(l, collapse = "\n"), sep = "\t")
  setnames(d, tolower(gsub("[^A-Za-z0-9]+", "_", names(d))))
  setnames(d, c("query_id", "target_id"), c("query", "target"), skip_absent = TRUE)
  d[, e_value := as.numeric(e_value)][]
}

## ---------------------------------------------------------------------------------------------
## Human-readable labels
## ---------------------------------------------------------------------------------------------
CENSAT_NAMES <- c(
  active_hor = "Alpha-satellite HOR, active (kinetochore site)",
  hor        = "Alpha-satellite HOR, inactive",
  dhor       = "Alpha-satellite HOR, divergent",
  mon        = "Alpha-satellite, monomeric",
  hsat1A     = "Human satellite 1A (HSat1A)",
  hsat1B     = "Human satellite 1B (HSat1B)",
  hsat1      = "Human satellite 1 (HSat1)",
  hsat2      = "Human satellite 2 (HSat2)",
  hsat3      = "Human satellite 3 (HSat3)",
  bsat       = "Beta-satellite (BSat)",
  gsat       = "Gamma-satellite (GSat)",
  censat     = "Other centromeric satellite",
  ct         = "Centromeric transition region",
  rDNA       = "rDNA array",
  ANY        = "Any centromeric/satellite region")

## CenSat class code -> readable name (unknown codes are returned unchanged)
pretty_censat <- function(x) {
  out <- unname(CENSAT_NAMES[x])
  miss <- is.na(out)
  ci <- match(tolower(x[miss]), tolower(names(CENSAT_NAMES)))
  out[miss] <- ifelse(is.na(ci), x[miss], unname(CENSAT_NAMES[ci]))
  out
}

pretty_repeat <- function(x) ifelse(x == "ANY", "Any RepeatMasker repeat", x)

pretty_label <- function(feature_set, label) {
  feature_set <- rep_len(feature_set, length(label))
  ifelse(feature_set == "censat", pretty_censat(label), pretty_repeat(label))
}

HOMER_NAMES <- c(`promoter-TSS` = "promoter/TSS", TTS = "transcription termination site (TTS)",
                 exon = "exon", intron = "intron", Intergenic = "intergenic", `non-coding` = "non-coding RNA",
                 `3' UTR` = "3' UTR", `5' UTR` = "5' UTR")
pretty_homer <- function(x) { y <- unname(HOMER_NAMES[x]); ifelse(is.na(y), x, y) }

## ---------------------------------------------------------------------------------------------
## Figure notes: one line per figure explaining what subset is shown (read by the report builder)
## ---------------------------------------------------------------------------------------------
fig_note <- function(cfg, name, text) {
  f <- file.path(cfg$figures, "figure_notes.tsv")
  old <- if (file.exists(f)) fread(f, sep = "\t", quote = "") else data.table(figure = character(), note = character())
  old <- old[figure != name]
  fwrite(rbind(old, data.table(figure = name, note = gsub("[\t\n]", " ", text))), f, sep = "\t", quote = FALSE)
}

## remove figures a script produced on a previous run (e.g. per-series files that no longer apply)
clear_figures <- function(cfg, pattern) {
  unlink(list.files(cfg$figures, pattern = pattern, full.names = TRUE))
}

## wrap long axis/facet labels (vectorised, base R)
stringr_wrap <- function(x, width) vapply(x, function(z) paste(strwrap(z, width), collapse = "\n"), "", USE.NAMES = FALSE)

## q/p-values for figure labels: values below double precision print as "< 2.2e-16" rather than 0
fmt_q <- function(q) ifelse(is.na(q), "NA", ifelse(q < 2.2e-16, "< 2.2e-16", formatC(q, digits = 2, format = "g")))

## ---------------------------------------------------------------------------------------------
## Heatmap rows: for every column (effect / covariate level / sample), keep every row that is significant in that
## column, or the column's `min_n` rows with the smallest p-values when fewer are significant; show the union.
## d: long table with one row per cell; row, col, p: column names; sig: name of a logical column.
## Returns the selected row keys ordered by their best p-value.
## ---------------------------------------------------------------------------------------------
select_heatmap_rows <- function(d, row, col, p, sig, min_n = 5) {
  x <- data.frame(r = as.character(d[[row]]), k = as.character(d[[col]]), p = as.numeric(d[[p]]), s = d[[sig]] %in% TRUE,
                  stringsAsFactors = FALSE)
  x <- x[!is.na(x$p), , drop = FALSE]
  if (!nrow(x)) return(character())
  keep <- unique(unlist(lapply(split(x, x$k), function(y) {
    y <- y[order(y$p), , drop = FALSE]
    s <- unique(y$r[y$s])
    if (length(s) >= min_n) s else head(unique(y$r), min_n)
  })))
  best <- tapply(x$p[x$r %in% keep], x$r[x$r %in% keep], min)
  names(sort(best))
}

## split a vector of heatmap rows into pages of at most max_rows (one figure file per page)
paginate <- function(rows, max_rows = 60) {
  if (!length(rows)) return(list())
  unname(split(rows, ceiling(seq_along(rows) / max_rows)))
}

## file name for page i of n (no suffix when there is only one page)
page_name <- function(base, i, n) if (n > 1) sprintf("%s_p%0*d", base, nchar(n), i) else base   # zero-padded so files sort in order
