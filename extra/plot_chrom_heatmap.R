## Peaks-per-chromosome heatmap for a figure panel, drawn at its printed size so the text comes out at a fixed point
## size. Same data and colours as the report's pk02 figure: one row per sample, one column per chromosome, tile
## colour = log2(share of the sample's peaks / share of the genome), tile text = number of peaks.
## Usage: Rscript extra/plot_chrom_heatmap.R <crossrun_out>/tables/table00b_peaks_by_chromosome.tsv [options]
##   --out FILE              .png, .pdf or .svg (default chrom_heatmap.pdf)
##   --width W               figure width in inches = the width it will be placed at (default 5.75)
##   --height H              figure height in inches (default: fitted to the number of samples)
##   --font PT               all text except the tile counts, in points (default 9)
##   --count-font PT|auto    tile counts in points (default auto: the largest size up to --count-max that fits the
##                           widest count in a tile, measured with the chosen font; 0 = no counts). With a fixed size
##                           that does not fit, a warning gives the size that would
##   --count-max PT          largest count size tried by auto (default 6); counts that need less than 4 pt are left out
##   --family NAME           font family for all text (default Arial, as in the other panels; "sans" = R's default)
##   --x-angle 0|90|auto     chromosome names upright or turned 90 degrees (default auto: upright when they fit a tile)
##   --labels-file TSV       columns sample_id, label (tab-separated, with that header). Also chooses the rows: only the
##                           samples listed are drawn, top to bottom in the file's order; the samples left out are
##                           listed in the log. A sample_id not in the table, a repeated sample_id or an empty label
##                           stops the script
##   --samples a,b,...       which samples and in what order, top to bottom (default: the --labels-file order, or as in
##                           the table). With --labels-file, --samples overrides its selection and order
##   --labels a=Label,...    row labels on the command line (default: sample_id); added to or replacing --labels-file
##                           labels. Samples left without a label are listed, since a raw sample ID widens the label column
##   --chroms chr1,chr2,...  which chromosomes and in what order (default: all, in genome order 1-22, X, Y, M)
##   --lim L                 colour scale from -L to L (default: as in the report, from the data, 2 to 6)
##   --legend bottom|right   colour key position (default bottom: saves width)
##   --dpi N                 PNG resolution (default 600)
## Tiles with no peaks are grey (log2 ratio undefined) and carry no count.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
a <- commandArgs(TRUE)
if (!length(a) || a[1] %in% c("-h", "--help")) {
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)); h <- readLines(f)
  cat(sub("^## ?", "", h[seq_len(which(!grepl("^##", h))[1] - 1)]), sep = "\n"); quit(status = 0)
}
opt <- function(n, d = NULL) { i <- match(n, a); if (is.na(i)) d else a[i + 1] }
die <- function(...) { message("ERROR: ", ...); quit(status = 1) }
known <- c("--out", "--width", "--height", "--font", "--count-font", "--count-max", "--family", "--x-angle", "--samples",
           "--labels", "--labels-file", "--chroms", "--lim", "--legend", "--dpi")
bad <- setdiff(grep("^--", a, value = TRUE), known); if (length(bad)) die("unknown option(s): ", paste(bad, collapse = ", "))
if (!file.exists(a[1])) die("not found: ", a[1])
D <- fread(a[1])
need <- c("sample_id", "chrom", "n_peaks", "log2_density_ratio")
if (!all(need %in% names(D))) die(a[1], " needs columns ", paste(need, collapse = ", "))
D[, `:=`(sample_id = as.character(sample_id), chrom = as.character(chrom))]
message(sprintf("[heatmap] read %d rows: %d samples x %d chromosomes", nrow(D), uniqueN(D$sample_id), uniqueN(D$chrom)))

L <- NULL
if (!is.null(opt("--labels-file"))) {
  if (!file.exists(opt("--labels-file"))) die("not found: ", opt("--labels-file"))
  L <- fread(opt("--labels-file"), sep = "\t", colClasses = "character")
  if (!all(c("sample_id", "label") %in% names(L))) die("--labels-file needs a header with columns sample_id and label (tab-separated)")
  L[, `:=`(sample_id = trimws(sample_id), label = trimws(label))]
  L <- L[nzchar(sample_id)]
  if (!nrow(L)) die("--labels-file lists no samples")
  dup <- unique(L$sample_id[duplicated(L$sample_id)]); if (length(dup)) die("--labels-file repeats: ", paste(dup, collapse = ", "))
  emp <- L$sample_id[is.na(L$label) | !nzchar(L$label)]; if (length(emp)) die("--labels-file has no label for: ", paste(emp, collapse = ", "))
  miss <- setdiff(L$sample_id, D$sample_id)
  if (length(miss)) die("--labels-file samples not in the table (check spelling and suffixes such as _REP1): ", paste(miss, collapse = ", "))
}
sm <- if (!is.null(opt("--samples"))) strsplit(opt("--samples"), ",")[[1]] else if (!is.null(L)) L$sample_id else unique(D$sample_id)
miss <- setdiff(sm, D$sample_id); if (length(miss)) die("samples not in the table: ", paste(miss, collapse = ", "))
if (anyDuplicated(sm)) die("--samples repeats: ", paste(unique(sm[duplicated(sm)]), collapse = ", "))
left <- setdiff(unique(D$sample_id), sm)
message(sprintf("[heatmap] drawing %d of %d samples (rows from %s)", length(sm), uniqueN(D$sample_id),
                if (!is.null(opt("--samples"))) "--samples" else if (!is.null(L)) "--labels-file" else "the table"))
if (length(left)) message("[heatmap] left out: ", paste(left, collapse = ", "))
message("[heatmap] rows, top to bottom: ", paste(sm, collapse = ", "))
rank_of <- function(x) { b <- sub("^chr", "", x); n <- suppressWarnings(as.numeric(b))
  ifelse(!is.na(n), n, ifelse(b == "X", 100, ifelse(b == "Y", 101, ifelse(b %in% c("M", "MT"), 102, 200)))) }
ch <- if (!is.null(opt("--chroms"))) strsplit(opt("--chroms"), ",")[[1]] else { u <- unique(D$chrom); u[order(rank_of(u), u)] }
miss <- setdiff(ch, D$chrom); if (length(miss)) die("chromosomes not in the table: ", paste(miss, collapse = ", "))
D <- D[sample_id %in% sm & chrom %in% ch]
if (nrow(D) != length(sm) * length(ch)) die(sprintf("expected %d sample x chromosome rows, found %d (duplicates or gaps in the table)",
                                                    length(sm) * length(ch), nrow(D)))

lab <- setNames(sm, sm)
if (!is.null(L)) { hit <- intersect(L$sample_id, sm); lab[hit] <- L$label[match(hit, L$sample_id)] }
if (!is.null(opt("--labels"))) { kv <- strsplit(strsplit(opt("--labels"), ",")[[1]], "=")
  k <- sapply(kv, `[`, 1); v <- sapply(kv, `[`, 2)
  miss <- setdiff(k, sm); if (length(miss)) die("--labels for samples not shown: ", paste(miss, collapse = ", "))
  lab[k] <- v }
unl <- names(lab)[lab == names(lab)]
if (length(unl)) message("[heatmap] no label given for (shown as the sample ID): ", paste(unl, collapse = ", "))

fs <- as.numeric(opt("--font", 9)); fam <- opt("--family", "Arial")
cf_opt <- opt("--count-font", "auto"); cmax <- as.numeric(opt("--count-max", 6)); xa_opt <- opt("--x-angle", "auto")
if (!xa_opt %in% c("0", "90", "auto")) die("--x-angle must be 0, 90 or auto")
if (cf_opt != "auto" && is.na(suppressWarnings(as.numeric(cf_opt)))) die("--count-font must be a number or auto")
W <- as.numeric(opt("--width", 5.75)); pos <- opt("--legend", "bottom")
if (!pos %in% c("bottom", "right")) die("--legend must be bottom or right")
lim <- if (!is.null(opt("--lim"))) as.numeric(opt("--lim")) else
  min(6, max(2, ceiling(max(abs(D$log2_density_ratio), na.rm = TRUE) * 2) / 2))     # as in R/05_peaks.R
lh <- fs * 1.25 / 72
H <- as.numeric(opt("--height", length(sm) * max(lh * 1.15, 0.16) + 3.4 * lh + (if (pos == "bottom") 3.2 * lh else 0)))
D[, `:=`(chrom = factor(chrom, ch), sample_id = factor(sample_id, rev(sm)))]
out <- opt("--out", "chrom_heatmap.pdf"); ext <- tolower(tools::file_ext(out))
if (!ext %in% c("png", "pdf", "svg")) die("--out must end in .png, .pdf or .svg")

## check that the font family exists (fontconfig, where available: Linux, and macOS with fontconfig installed)
if (fam != "sans" && nzchar(Sys.which("fc-match"))) {
  m <- system2("fc-match", c("-f", "'%{family}'", shQuote(fam)), stdout = TRUE)
  mm <- paste(m, collapse = " ")
  if (tolower(fam) == "arial" && grepl("Liberation Sans|Arimo", mm))
    message(sprintf("[heatmap] note: Arial not installed; using %s, which has the same character widths", mm))
  else if (!grepl(fam, mm, ignore.case = TRUE))
    message(sprintf("[heatmap] WARNING: font '%s' not found; fontconfig substitutes %s, which has different widths. Install it or pass --family", fam, paste(m, collapse = " ")))
}

make_plot <- function(xa, cf) {
  p <- ggplot(D, aes(chrom, sample_id)) +
    geom_tile(aes(fill = pmax(pmin(log2_density_ratio, lim), -lim)), colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0, limits = c(-lim, lim),
                         na.value = "grey85", name = "log2(share of peaks / share of genome)",
                         breaks = c(-lim, 0, lim), labels = function(x) format(x, drop0trailing = TRUE, trim = TRUE)) +
    scale_y_discrete(labels = lab, expand = c(0, 0)) + scale_x_discrete(labels = function(x) sub("^chr", "", x), expand = c(0, 0)) +
    labs(x = "Chromosome", y = NULL) +
    theme_minimal(base_size = fs, base_family = fam) +
    theme(text = element_text(size = fs, colour = "black", family = fam), axis.text = element_text(size = fs, colour = "black"),
          axis.text.x = if (xa == 90) element_text(angle = 90, vjust = 0.5, hjust = 1) else element_text(),
          panel.grid = element_blank(), plot.margin = margin(2, 2, 2, 2), legend.position = pos,
          legend.title = element_text(size = fs), legend.text = element_text(size = fs),
          legend.key.height = unit(if (pos == "bottom") fs * 0.8 else fs * 1.2, "pt"),
          legend.key.width = unit(if (pos == "bottom") fs * 1.6 else fs * 0.8, "pt"),   # a colour bar is 5 keys long
          legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(fs * 0.4, "pt")) +
    guides(fill = guide_colourbar(title.position = if (pos == "bottom") "left" else "top", title.vjust = 0.9))
  if (cf > 0) p <- p + geom_text(aes(label = ifelse(n_peaks > 0, n_peaks, "")), size = cf / .pt, family = fam)
  p
}

## measure with the real font on a device of the final size: tile size from the laid-out plot, text widths from the font
use_cairo <- capabilities("cairo")
mdev <- tempfile(fileext = ".pdf")
if (use_cairo) grDevices::cairo_pdf(mdev, width = W, height = H, family = fam) else grDevices::pdf(mdev, width = W, height = H)
tile <- function(p) {
  g <- ggplotGrob(p); pc <- unique(g$layout[grepl("^panel", g$layout$name), c("l", "t")])
  aw <- sum(sapply(seq_along(g$widths)[-pc$l], function(i) grid::convertWidth(g$widths[i], "in", valueOnly = TRUE)))
  ah <- sum(sapply(seq_along(g$heights)[-pc$t], function(i) grid::convertHeight(g$heights[i], "in", valueOnly = TRUE)))
  c(w = (W - aw) / length(ch), h = (H - ah) / length(sm))
}
txt_w <- function(s, pt) max(sapply(s, function(x) grid::convertWidth(grid::stringWidth(x), "in", valueOnly = TRUE) *
                                       pt / grid::get.gpar("fontsize")$fontsize))
grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = fam, fontsize = 10)))
t0 <- tile(make_plot(0, 0))
chr_w <- txt_w(sub("^chr", "", ch), fs)
xa <- if (xa_opt == "auto") (if (chr_w <= 0.8 * t0["w"]) 0 else 90) else as.numeric(xa_opt)   # 0.8: keep a gap between neighbours
t1 <- tile(make_plot(xa, 0))
if (t1["w"] < 2 / 72) { grid::popViewport(); invisible(dev.off()); unlink(mdev)
  die(sprintf("the row labels and key leave %.1f pt per chromosome at --width %g; shorten the labels or widen the figure", t1["w"] * 72, W)) }
cnt <- as.character(D$n_peaks[D$n_peaks > 0]); cnt_w1 <- txt_w(cnt[which.max(nchar(cnt))], 1)   # widest count at 1 pt
fit <- floor(min(0.85 * t1["w"] / cnt_w1, 0.75 * t1["h"] * 72) * 4) / 4                    # largest size that fits, 0.25 pt steps
grid::popViewport(); invisible(dev.off()); unlink(mdev)
cf <- if (cf_opt == "auto") { if (fit >= 4) min(cmax, fit) else 0 } else as.numeric(cf_opt)
if (cf_opt == "auto" && fit < 4) message(sprintf("[heatmap] WARNING: counts would need %.2f pt to fit a tile; left out. Shorten the labels or widen the figure", fit))
if (cf_opt != "auto" && cf > fit) message(sprintf("[heatmap] WARNING: counts at %g pt are wider than a tile; %.2f pt would fit (or use --count-font auto)", cf, fit))
if (xa == 0 && chr_w > 0.8 * t1["w"]) message("[heatmap] WARNING: chromosome names are wider than a tile; use --x-angle 90 or auto")

p <- make_plot(xa, cf)
dev <- switch(ext, png = "png", pdf = if (use_cairo) grDevices::cairo_pdf else "pdf", svg = if (use_cairo) grDevices::svg else "svg")
ggsave(out, p, width = W, height = H, dpi = as.numeric(opt("--dpi", 600)), bg = "white", device = dev)
message(sprintf("[heatmap] %d samples x %d chromosomes -> %s (%.2f x %.2f in; %s %g pt; tiles %.1f x %.1f pt; chromosome names %s; counts %s)",
                length(sm), length(ch), out, W, H, fam, fs, t1["w"] * 72, t1["h"] * 72, if (xa == 0) "upright" else "turned 90",
                if (cf > 0) paste(cf, "pt") else "off"))
