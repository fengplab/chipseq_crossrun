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
##   --no-totals             leave out the totals (default: a Total column with each sample's peaks, a Total row with
##                           each chromosome's peaks over the samples drawn, and the grand total where they meet; totals
##                           cover only the samples and chromosomes drawn, in bold on grey tiles outside the colour scale)
##   --total-font PT|auto    totals in points (default auto: the same size as the counts). The Total column is as wide
##                           as the widest sample total needs; chromosome totals too wide for their tile are turned
##                           90 degrees and the Total row is made taller (the figure grows unless --height is given)
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
           "--labels", "--labels-file", "--chroms", "--lim", "--legend", "--dpi", "--no-totals", "--total-font")
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
tot_on <- !"--no-totals" %in% a; tf_opt <- opt("--total-font", "auto")
if (!xa_opt %in% c("0", "90", "auto")) die("--x-angle must be 0, 90 or auto")
for (o in c("--count-font", "--total-font")) { v <- opt(o, "auto")
  if (v != "auto" && is.na(suppressWarnings(as.numeric(v)))) die(o, " must be a number or auto") }
W <- as.numeric(opt("--width", 5.75)); pos <- opt("--legend", "bottom")
if (!pos %in% c("bottom", "right")) die("--legend must be bottom or right")
lim <- if (!is.null(opt("--lim"))) as.numeric(opt("--lim")) else
  min(6, max(2, ceiling(max(abs(D$log2_density_ratio), na.rm = TRUE) * 2) / 2))     # as in R/05_peaks.R
lh <- fs * 1.25 / 72
nr <- length(sm) + (if (tot_on) 1.3 else 0)                                          # tile rows, incl. the Total row and its gap
H <- as.numeric(opt("--height", nr * max(lh * 1.15, 0.16) + 3.4 * lh + (if (pos == "bottom") 3.2 * lh else 0)))
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

## layout in tile units. Chromosome i is at x = i; the Total column starts GAP after the last chromosome tile and is TW
## tiles wide. Sample k (from the top) is at y = ns - k + 1, all shifted up by TH + GAP when the Total row (TH tiles high,
## from y = 0.5) is drawn. TH > 1 only when the chromosome totals are turned 90 degrees to fit their tiles.
GAP <- 0.3; nch <- length(ch); ns <- length(sm)
D[, `:=`(x = match(chrom, ch), yb = ns - match(sample_id, sm) + 1)]
RT <- D[, .(n = sum(n_peaks)), by = .(sample_id, yb)]                                # per sample, over the chromosomes drawn
CT <- D[, .(n = sum(n_peaks)), by = .(chrom, x)]                                     # per chromosome, over the samples drawn
grand <- sum(D$n_peaks)
stopifnot(sum(RT$n) == grand, sum(CT$n) == grand)
if (tot_on) {
  message(sprintf("[heatmap] totals: %s peaks over %d samples x %d chromosomes; per sample %s to %s, per chromosome %s to %s",
                  format(grand, big.mark = ","), ns, nch, format(min(RT$n), big.mark = ","), format(max(RT$n), big.mark = ","),
                  format(min(CT$n), big.mark = ","), format(max(CT$n), big.mark = ",")))
  if (!is.null(opt("--chroms"))) message("[heatmap] note: sample totals count only the chromosomes drawn (--chroms)")
}

make_plot <- function(xa, cf, tf, tw, th, trot) {
  sh <- if (tot_on) th + GAP else 0
  xl <- nch + 0.5 + GAP; xt <- xl + tw / 2; xmax <- if (tot_on) xl + tw else nch + 0.5
  yt <- 0.5 + th / 2; ymax <- ns + 0.5 + sh
  DD <- copy(D)[, y := yb + sh]
  p <- ggplot() +
    geom_tile(data = DD, aes(x, y, fill = pmax(pmin(log2_density_ratio, lim), -lim)), width = 1, height = 1, colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0, limits = c(-lim, lim),
                         na.value = "grey85", name = "log2(share of peaks / share of genome)",
                         breaks = c(-lim, 0, lim), labels = function(x) format(x, drop0trailing = TRUE, trim = TRUE)) +
    scale_x_continuous(breaks = c(seq_len(nch), if (tot_on) xt), labels = c(sub("^chr", "", ch), if (tot_on) "Total"),
                       limits = c(0.5, xmax), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(ns - seq_len(ns) + 1 + sh, if (tot_on) yt), labels = c(unname(lab[sm]), if (tot_on) "Total"),
                       limits = c(0.5, ymax), expand = c(0, 0)) +
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
  if (cf > 0) p <- p + geom_text(data = DD, aes(x, y, label = ifelse(n_peaks > 0, n_peaks, "")), size = cf / .pt, family = fam)
  if (tot_on) {
    TC <- rbind(RT[, .(x = xt, y = yb + sh, n, w = tw, h = 1)], data.table(x = xt, y = yt, n = grand, w = tw, h = th))  # Total column
    TR <- CT[, .(x, y = yt, n, w = 1, h = th)]                                                                          # Total row
    p <- p + geom_tile(data = rbind(TC, TR), aes(x, y, width = w, height = h), fill = "#E6E6E6", colour = "white", linewidth = 0.3)
    if (tf > 0) p <- p + geom_text(data = TC, aes(x, y, label = n), size = tf / .pt, family = fam, fontface = "bold") +
      geom_text(data = TR, aes(x, y, label = n), size = tf / .pt, family = fam, fontface = "bold", angle = if (trot) 90 else 0)
  }
  p
}

## measure with the real font on a device of the final size: panel size from the laid-out plot, text widths from the font
use_cairo <- capabilities("cairo")
open_dev <- function() {
  mdev <<- tempfile(fileext = ".pdf")
  if (use_cairo) grDevices::cairo_pdf(mdev, width = W, height = H, family = fam) else grDevices::pdf(mdev, width = W, height = H)
  grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = fam, fontsize = 10)))
}
close_dev <- function() { grid::popViewport(); invisible(dev.off()); unlink(mdev) }
panel <- function(p) {
  g <- ggplotGrob(p); pc <- unique(g$layout[grepl("^panel", g$layout$name), c("l", "t")])
  aw <- sum(sapply(seq_along(g$widths)[-pc$l], function(i) grid::convertWidth(g$widths[i], "in", valueOnly = TRUE)))
  ah <- sum(sapply(seq_along(g$heights)[-pc$t], function(i) grid::convertHeight(g$heights[i], "in", valueOnly = TRUE)))
  c(w = W - aw, h = H - ah)
}
txt_w <- function(s, pt, bold = FALSE) {
  grid::pushViewport(grid::viewport(gp = grid::gpar(fontface = if (bold) 2 else 1)))
  on.exit(grid::popViewport())
  max(sapply(s, function(x) grid::convertWidth(grid::stringWidth(x), "in", valueOnly = TRUE) * pt / grid::get.gpar("fontsize")$fontsize))
}
widest <- function(v) as.character(v[which.max(nchar(v))])
open_dev()
pw <- unname(panel(make_plot(0, 0, 0, 1, 1, FALSE))["w"])             # the panel width does not depend on the Total column's width
ph0 <- unname(panel(make_plot(90, 0, 0, 1, 1, FALSE))["h"])           # first estimate of the panel height (names turned)
uy0 <- ph0 / (ns + (if (tot_on) 1 + GAP else 0))
cnt <- as.character(D$n_peaks[D$n_peaks > 0]); cnt_w1 <- txt_w(widest(cnt), 1)            # widest count at 1 pt
## Counts and totals share one size: the largest z (0.25 pt steps, up to --count-max) at which a count fits a chromosome
## tile once the Total column is wide enough for the widest sample total and the grand total at z (bold)
rt_w1 <- if (tot_on) txt_w(widest(c(RT$n, grand)), 1, TRUE) else 0
ct_w1 <- if (tot_on) txt_w(widest(CT$n), 1, TRUE) else 0
geom_x <- function(z) { col_in <- if (tot_on) z * rt_w1 / 0.85 else 0
  u <- (pw - col_in) / (nch + (if (tot_on) GAP else 0)); tw <- if (tot_on) max(1, col_in / u) else 0
  list(u = pw / (nch + (if (tot_on) GAP + tw else 0)), tw = tw) }
zs <- rev(seq(4, max(4, cmax), by = 0.25)); z_ok <- NA
for (z in zs) { g <- geom_x(z); if (g$u > 0 && z * cnt_w1 <= 0.85 * g$u && z <= 0.75 * uy0 * 72) { z_ok <- z; break } }
cf <- if (cf_opt == "auto") { if (is.na(z_ok)) 0 else z_ok } else as.numeric(cf_opt)
tf <- if (!tot_on) 0 else if (tf_opt == "auto") (if (is.na(z_ok)) 4 else z_ok) else as.numeric(tf_opt)
gx <- geom_x(max(tf, 4)); unit_x <- gx$u; tw <- gx$tw                 # the Total column is sized for the totals only
if (unit_x < 2 / 72) { close_dev()
  die(sprintf("the row labels and key leave %.1f pt per chromosome at --width %g; shorten the labels or widen the figure", unit_x * 72, W)) }
## chromosome totals: upright when they fit a tile, otherwise turned 90 degrees in a taller Total row
trot <- tot_on && tf * ct_w1 > 0.85 * unit_x
th <- if (trot) max(1, (tf * ct_w1 / 0.85) / uy0) else 1
if (trot && is.null(opt("--height"))) { H <- H + (th - 1) * uy0; close_dev(); open_dev() }   # grow the figure, not shrink the rows
chr_w <- txt_w(sub("^chr", "", ch), fs)
xa <- if (xa_opt == "auto") (if (chr_w <= 0.8 * unit_x) 0 else 90) else as.numeric(xa_opt)   # 0.8: keep a gap between neighbours
ph <- unname(panel(make_plot(xa, 0, 0, tw, th, trot))["h"])
unit_y <- ph / (ns + (if (tot_on) th + GAP else 0))
fit <- floor(min(0.85 * unit_x / cnt_w1, 0.75 * unit_y * 72) * 4) / 4                      # largest count size that fits
close_dev()
if (cf_opt == "auto" && is.na(z_ok)) message(sprintf("[heatmap] WARNING: counts would need %.2f pt to fit a tile; left out. Shorten the labels or widen the figure", fit))
if (cf_opt != "auto" && cf > fit) message(sprintf("[heatmap] WARNING: counts at %g pt are wider than a tile; %.2f pt would fit (or use --count-font auto)", cf, fit))
if (tot_on && trot && tf > 0.85 * unit_x * 72) message(sprintf("[heatmap] WARNING: chromosome totals at %g pt do not fit their tiles even turned; use a smaller --total-font or widen the figure", tf))
if (xa == 0 && chr_w > 0.8 * unit_x) message("[heatmap] WARNING: chromosome names are wider than a tile; use --x-angle 90 or auto")

p <- make_plot(xa, cf, tf, tw, th, trot)
dev <- switch(ext, png = "png", pdf = if (use_cairo) grDevices::cairo_pdf else "pdf", svg = if (use_cairo) grDevices::svg else "svg")
ggsave(out, p, width = W, height = H, dpi = as.numeric(opt("--dpi", 600)), bg = "white", device = dev)
message(sprintf("[heatmap] %d samples x %d chromosomes -> %s (%.2f x %.2f in; %s %g pt; tiles %.1f x %.1f pt; chromosome names %s; counts %s; totals %s)",
                ns, nch, out, W, H, fam, fs, unit_x * 72, unit_y * 72, if (xa == 0) "upright" else "turned 90",
                if (cf > 0) paste(cf, "pt") else "off", if (tot_on) sprintf("%g pt bold, Total column %.1f pt wide, chromosome totals %s", tf, tw * unit_x * 72, if (trot) "turned 90" else "upright") else "off"))
