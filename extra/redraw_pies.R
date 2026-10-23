## Redraw pie charts from the TSV that R/plot_pies.R or extra/compare_pies.R writes next to each figure (columns panel,
## contig, value; other columns are ignored), for a figure panel: pies in a grid, no legend, every slice labelled.
## Slices of at least --inside-min percent carry their chromosome name and percentage inside; smaller slices are
## labelled outside ("chr18 1.2%"), stacked down each side of the pie without overlapping, with a leader line to the slice.
## Usage: Rscript extra/redraw_pies.R pies.tsv [options]
##   --out FILE            .png (default), .pdf or .svg (default <tsv>_labelled.png)
##   --ncol N              pies per row (default: all in one row)
##   --width W --height H  figure size in inches (default 2.4 per pie wide, 1.4 per row high)
##   --font PT             label size in points (default 5.5); titles are 1.15 x this, bold
##   --inside-min PCT      smallest slice labelled inside (default 6)
##   --radius F            pie radius as a fraction of the cell width (default 0.22; capped by the cell height)
##   --dpi N               PNG resolution (default 600)
##   --panels a,b,...      which panels and in what order (default: as in the file)
##   --totals a=n,b=m      peak total per panel, added to its title (optional; a column peaks_total is used if present)
##   --rds FILE            also save the figure as an R object (patchwork), to place it in a larger figure without re-scaling
## Colours and slice order are those of R/plot_pies.R: chromosomes clockwise from 12 o'clock in genome order.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
a <- commandArgs(TRUE)
if (!length(a) || a[1] %in% c("-h", "--help")) {
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)); h <- readLines(f)
  cat(sub("^## ?", "", h[seq_len(which(!grepl("^##", h))[1] - 1)]), sep = "\n"); quit(status = 0)
}
opt <- function(n, d = NULL) { i <- match(n, a); if (is.na(i)) d else a[i + 1] }
die <- function(...) { message("ERROR: ", ...); quit(status = 1) }
known <- c("--out", "--ncol", "--width", "--height", "--font", "--inside-min", "--dpi", "--panels", "--totals", "--radius", "--rds")
bad <- setdiff(grep("^--", a, value = TRUE), known); if (length(bad)) die("unknown option(s): ", paste(bad, collapse = ", "))
D <- fread(a[1]); if (!all(c("panel", "contig", "value") %in% names(D))) die(a[1], " needs columns panel, contig, value")
D[, panel := as.character(panel)]
pn <- if (!is.null(opt("--panels"))) strsplit(opt("--panels"), ",")[[1]] else unique(D$panel)
miss <- setdiff(pn, D$panel); if (length(miss)) die("panels not in the file: ", paste(miss, collapse = ", "))
D <- D[panel %in% pn]
tot <- if ("peaks_total" %in% names(D)) D[, .(t = peaks_total[1]), by = panel][, setNames(t, panel)] else NULL
if (!is.null(opt("--totals"))) { kv <- strsplit(strsplit(opt("--totals"), ",")[[1]], "="); tot <- setNames(as.numeric(sapply(kv, `[`, 2)), sapply(kv, `[`, 1)) }
ncol <- as.integer(opt("--ncol", length(pn))); nrow <- ceiling(length(pn) / ncol)
fs <- as.numeric(opt("--font", 5.5)); inside_min <- as.numeric(opt("--inside-min", 6))
W <- as.numeric(opt("--width", 2.4 * ncol)); H <- as.numeric(opt("--height", 1.4 * nrow))
out <- opt("--out", sub("\\.tsv$", "_labelled.png", a[1]))

## genome order and the colours of R/plot_pies.R
rank_of <- function(x) { n <- suppressWarnings(as.numeric(sub("^chr", "", x))); ifelse(!is.na(n), n, ifelse(x == "chrX", 100, ifelse(x == "chrY", 101, 200))) }
ctg <- unique(D$contig); ctg <- ctg[order(rank_of(ctg), ctg)]
cols <- setNames(grDevices::hcl(h = seq(15, 375, length.out = length(ctg) + 1)[-1], c = 45, l = rep(c(62, 78), length.out = length(ctg))), ctg)

cw <- W / ncol; ch <- H / nrow                        # one cell, inches
lh <- fs * 1.25 / 72                                  # line height of a label, inches
fmt <- function(p) ifelse(p >= 10, sprintf("%.0f%%", p), ifelse(p >= 1, sprintf("%.1f%%", p), sprintf("%.2g%%", p)))
spread <- function(y, gap, lo, hi) {                   # keep order, at least `gap` apart, as close to y as possible
  n <- length(y); if (!n) return(y); o <- order(y); ys <- y[o]; off <- (seq_len(n) - 1) * gap
  blk <- list(); for (q in seq_len(n)) { blk[[length(blk) + 1]] <- c(q, q)
    repeat { nb <- length(blk); cur <- blk[[nb]]; ix <- cur[1]:cur[2]; pos <- mean(ys[ix] - off[ix]) + off[ix]
      if (nb == 1) break; pv <- blk[[nb - 1]]; px <- pv[1]:pv[2]; ppos <- mean(ys[px] - off[px]) + off[px]
      if (tail(ppos, 1) + gap <= pos[1] + 1e-9) break; blk[[nb - 1]] <- c(pv[1], cur[2]); blk[[nb]] <- NULL } }
  r <- numeric(n); for (b in blk) { ix <- b[1]:b[2]; r[ix] <- mean(ys[ix] - off[ix]) + off[ix] }
  if (r[1] < lo) r <- r + (lo - r[1]); if (r[n] > hi) r <- r - (r[n] - hi); out <- numeric(n); out[o] <- r; out
}
one_pie <- function(pnl) {
  d <- D[panel == pnl][, .(value = sum(value)), by = contig]
  d <- d[match(ctg, contig)][!is.na(contig)][value > 0]
  d[, share := 100 * value / sum(value)]
  d[, `:=`(a1 = 2 * pi * (cumsum(share) - share) / 100, a2 = 2 * pi * cumsum(share) / 100)]
  d[, am := (a1 + a2) / 2]                            # clockwise from 12 o'clock
  title_h <- 2.2 * lh * 1.15
  R <- min(cw * as.numeric(opt("--radius", 0.22)), (ch - title_h) * 0.40)   # pie radius, inches: room for labels on both sides
  poly <- d[, { t <- seq(a1, a2, length.out = max(3, ceiling((a2 - a1) / 0.02))); .(x = c(0, R * sin(t)), y = c(0, R * cos(t))) }, by = contig]
  p <- ggplot() + geom_polygon(data = poly, aes(x, y, group = contig, fill = contig), colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = cols, guide = "none")
  ins <- d[share >= inside_min]
  ins[, lab := paste0(contig, "\n", fmt(share))]   # labels computed here, so a saved figure (--rds) needs no helper functions
  if (nrow(ins)) p <- p + geom_text(data = ins, aes(x = 0.64 * R * sin(am), y = 0.64 * R * cos(am), label = lab),
                                    size = fs / .pt, lineheight = 0.9, colour = "grey10")
  outl <- d[share < inside_min][, side := ifelse(sin(am) >= 0, 1, -1)]
  if (nrow(outl)) {
    outl[, ty := 1.12 * R * cos(am)]
    lim <- (ch - title_h) / 2 - lh / 2
    cap <- floor(2 * lim / lh) + 1                     # labels that fit down one side of the cell
    if (nrow(outl) > 2 * cap) message("[pies] WARNING: ", pnl, ": ", nrow(outl), " outside labels but room for ", 2 * cap,
                                      "; labels will overlap (raise --height, lower --font or raise --inside-min)")
    ## a side with more labels than fit hands its lowest labels (the slices nearest 6 o'clock) to the other side
    outl[, moved := FALSE]
    while (sum(outl$side < 0) > cap && sum(outl$side > 0) < cap) { i <- which(outl$side < 0 & !outl$moved); outl[i[which.min(outl$am[i])], `:=`(side = 1, moved = TRUE)] }
    while (sum(outl$side > 0) > cap && sum(outl$side < 0) < cap) { i <- which(outl$side > 0 & !outl$moved); outl[i[which.max(outl$am[i])], `:=`(side = -1, moved = TRUE)] }
    ## stack each side in clockwise order (top to bottom), so handed-over labels go to the bottom of the stack;
    ## strictly decreasing targets keep that order through spread()
    outl[, ord := ifelse(side > 0, am, -am)]
    setorder(outl, side, ord)
    outl[, ty := cummin(ty) - seq_len(.N) * 1e-6, by = side]
    outl[, ly := spread(ty, lh, -lim, lim), by = side]
    ## leaders: a handed-over label's leader leaves its slice radially, follows its own arc round the pie (lower labels on
    ## wider arcs, so leaders never cross each other or the pie) and ends with a horizontal run to the label
    outl[, k := ifelse(moved, cumsum(moved), 0L), by = side]
    dr <- 0.03 * R
    rmax <- outl[, .(r = 1.08 * R + max(k) * dr), by = side]
    outl[rmax, lx := side * max(1.22 * R, r + 0.08 * R), on = "side", by = .EACHI]
    seg <- rbindlist(lapply(seq_len(nrow(outl)), function(j) {
      L <- outl[j]
      if (!L$moved) return(data.table(contig = L$contig, x = c(R * sin(L$am), 1.08 * R * sin(L$am), L$lx),
                                      y = c(R * cos(L$am), 1.08 * R * cos(L$am), L$ly)))
      r <- 1.08 * R + L$k * dr
      c0 <- acos(max(-1, min(1, L$ly / r)))                       # angle (from 12 o'clock) of the arc point at the label height
      th <- if (L$side > 0) c0 else 2 * pi - c0
      t <- seq(L$am, th, length.out = max(3, ceiling(abs(th - L$am) / 0.02)))
      data.table(contig = L$contig, x = c(R * sin(L$am), r * sin(t), L$lx), y = c(R * cos(L$am), r * cos(t), L$ly))
    }))
    outl[, lab := paste(contig, fmt(share))]
    p <- p + geom_path(data = seg, aes(x, y, group = contig), linewidth = 0.15, colour = "grey40") +
      geom_text(data = outl, aes(x = lx + side * 0.03, y = ly, label = lab, hjust = ifelse(side > 0, 0, 1)),
                size = fs / .pt, colour = "grey10")
  }
  ttl <- if (!is.null(tot) && !is.na(tot[pnl])) sprintf("%s\n%s peaks", pnl, format(tot[pnl], big.mark = ",")) else pnl
  p + annotate("text", x = 0, y = (ch - title_h) / 2 + title_h * 0.45, label = ttl, size = fs * 1.15 / .pt, fontface = "bold", lineheight = 0.95) +
    coord_fixed(xlim = c(-cw / 2, cw / 2), ylim = c(-(ch - title_h) / 2, (ch - title_h) / 2 + title_h), expand = FALSE, clip = "off") +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}
fig <- wrap_plots(lapply(pn, one_pie), ncol = ncol)
ggsave(out, fig, width = W, height = H, dpi = as.numeric(opt("--dpi", 600)), bg = "white")
if (!is.null(opt("--rds"))) { saveRDS(fig, opt("--rds")); message("[pies] figure object -> ", opt("--rds")) }
message("[pies] wrote ", out, sprintf(" (%.2f x %.2f in)", W, H))
